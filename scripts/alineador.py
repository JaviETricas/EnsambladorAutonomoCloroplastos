#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import platform
import re
import subprocess
import sys


def fasta_to_oneline_py(src: Path, dst: Path) -> None:
    """Convierte cualquier FASTA a formato 'una línea por secuencia' (cabecera + secuencia sin saltos)."""
    with src.open("r", encoding="utf-8", errors="replace") as fi, \
         dst.open("w", encoding="utf-8") as fo:
        seq_parts = []
        header = None
        for line in fi:
            if line.startswith(">"):
                if header is not None:
                    fo.write("".join(seq_parts) + "\n")
                header = line.rstrip("\n")
                fo.write(header + "\n")
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            fo.write("".join(seq_parts) + "\n")


def fasta_to_oneline_awk(src: Path, dst: Path) -> None:
    """Versión AWK (Unix)."""
    cmd = [
        "awk",
        r'/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {print seq}',
        str(src)
    ]
    with dst.open("w", encoding="utf-8") as fo:
        subprocess.run(cmd, stdout=fo, check=True)


def get_oneliner() -> callable:
    """Elige implementación oneline (Python en Windows, AWK en Unix si está disponible)."""
    if platform.system() == "Windows":
        return fasta_to_oneline_py
    # intenta AWK, cae a Python si falla
    try:
        subprocess.run(["awk", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return fasta_to_oneline_awk
    except Exception:
        return fasta_to_oneline_py


def normalize_basename(name: str) -> str:
    """
    Limpia sufijos comunes de nombres generados (Novowrap, .fastq.gz, etc.)
    para usarlo como 'hint' de cabecera.
    """
    n = name
    n = re.sub(r"Option_\d+_", "", n)
    n = re.sub(r"_Novowrap.*", "", n)
    n = re.sub(r"_pileup", "", n)
    n = re.sub(r"\.fastq(?:\.gz)?_?", "_", n)
    n = re.sub(r"\.(fasta|fa|fna|tsv|txt|fai)$", "", n, flags=re.IGNORECASE)
    return n


def load_mutations(tsv_path: Path) -> dict[int, str]:
    """
    Lee el TSV con cabecera:
    CHROM, POS, ., ,, A, C, G, T, N, +, -, ^, $, other
    Devuelve un dict pos(1-based) -> base (A/C/G/T) a corregir.
    """
    muts: dict[int, str] = {}
    with tsv_path.open("r", encoding="utf-8", errors="replace") as tsv:
        for line in tsv:
            if not line.strip() or line.startswith("#") or line.startswith("CHROM"):
                continue
            cols = line.rstrip().split("\t")
            try:
                pos   = int(cols[1])  # POS
                dot   = int(cols[2])  # '.'
                comma = int(cols[3])  # ','
                A = int(cols[4]); C = int(cols[5]); G = int(cols[6]); T = int(cols[7])
            except (ValueError, IndexError):
                continue

            total = dot + comma + A + C + G + T
            if total == 0:
                continue

            # si la referencia no domina, elegir mayoría A/C/G/T
            if (dot + comma == 0) or ((dot + comma) / total < 0.5):
                bases  = ("A", "C", "G", "T")
                counts = (A, C, G, T)
                muts[pos] = bases[max(range(4), key=lambda i: counts[i])]
    return muts


def align_with_reference_and_number(oneline_fasta: Path) -> Path | None:
    """
    Alinea oneline_fasta contra:
      1) El Fusion_Cloroplastos_<max>.fasta en ROOT/Resultados, o
      2) ROOT/temporalDocs/cloroplasto_referencia/Fusion_Cloroplastos_1.fasta (fallback),
    y escribe el resultado en ROOT/Resultados/Fusion_Cloroplastos_<next>.fasta.

    Implementación: MAFFT requiere un único archivo de entrada; se crea un FASTA temporal
    concatenando referencia + query y se pasa como único input a MAFFT.
    """
    ROOT_DIR = Path(__file__).resolve().parents[1]
    LIB_DIR  = ROOT_DIR / "libreris"
    RES_DIR  = ROOT_DIR / "Resultados"
    REF_DIR  = ROOT_DIR / "temporalDocs" / "cloroplasto_referencia"

    RES_DIR.mkdir(parents=True, exist_ok=True)

    # --- elegir referencia ---
    import re, os
    pattern = re.compile(r"^Fusion_Cloroplastos_(\d+)\.fasta$", re.IGNORECASE)
    existentes: list[tuple[int, Path]] = []
    for p in RES_DIR.iterdir():
        if p.is_file():
            m = pattern.match(p.name)
            if m:
                try:
                    existentes.append((int(m.group(1)), p))
                except ValueError:
                    pass

    if existentes:
        existentes.sort(key=lambda x: x[0])
        max_n, ref_fasta = existentes[-1]
        next_n = max_n + 1
    else:
        ref_fasta = REF_DIR / "Fusion_Cloroplastos_1.fasta"
        if not ref_fasta.exists():
            print(f"[WARN] No hay previos en {RES_DIR} y no existe referencia: {ref_fasta}")
            return None
        next_n = 1

    out_path = RES_DIR / f"Fusion_Cloroplastos_{next_n}.fasta"

    # --- comprobaciones previas ---
    for label, f in (("referencia", ref_fasta), ("query", oneline_fasta)):
        if not f.exists():
            print(f"[WARN] El FASTA de {label} no existe: {f}")
            return None
        if f.stat().st_size == 0:
            print(f"[WARN] El FASTA de {label} está vacío: {f}")
            return None

    # --- preparar input temporal para MAFFT (un solo archivo) ---
    tmp_in = RES_DIR / f".tmp_mafft_input_{next_n}.fasta"
    try:
        with tmp_in.open("w", encoding="utf-8") as fo:
            with ref_fasta.open("r", encoding="utf-8", errors="replace") as fr:
                fo.write(fr.read())
                if not fo.tell() or not str(fo).endswith("\n"):
                    fo.write("\n")
            with oneline_fasta.open("r", encoding="utf-8", errors="replace") as fq:
                fo.write(fq.read())
                fo.write("\n")
    except Exception as e:
        print(f"[WARN] No se pudo crear el input temporal para MAFFT: {e}")
        return None

    # --- resolver mafft binario ---
    def resolve_mafft() -> list[str]:
        cand = LIB_DIR / "mafft"
        if cand.exists() and cand.is_file() and os.access(cand, os.X_OK):
            return [str(cand)]
        return ["mafft"]

    def run_mafft(mafft_cmd: list[str]) -> tuple[int, str, str]:
        cmd = mafft_cmd + ["--auto", str(tmp_in)]
        print(f"[MAFFT] {' '.join(cmd)} -> {out_path}")
        res = subprocess.run(cmd, text=True, capture_output=True)
        return res.returncode, res.stdout, res.stderr

    # intento 1: libreris/mafft (si ejecutable) o PATH
    mafft_cmd = resolve_mafft()
    rc, out, err = run_mafft(mafft_cmd)

    # si falla y era libreris/mafft, reintentar con PATH
    if rc != 0 and mafft_cmd[0] != "mafft":
        print(f"[WARN] MAFFT en libreris falló (rc={rc}). stderr:\n{err}")
        rc, out, err = run_mafft(["mafft"])

    # limpieza del temporal
    try:
        tmp_in.unlink(missing_ok=True)
    except Exception:
        pass

    if rc != 0:
        print(f"[WARN] MAFFT falló (rc={rc}). stderr:\n{err}")
        return None

    try:
        with out_path.open("w", encoding="utf-8") as fo:
            fo.write(out)
        print(f"[OK] Alineado: {out_path}")
        return out_path
    except Exception as e:
        print(f"[WARN] No se pudo escribir el alineado: {e}")
        return None



def write_oneline_copy(src_fasta: Path, dst_oneline: Path) -> None:
    oneliner = get_oneliner()
    dst_oneline.parent.mkdir(parents=True, exist_ok=True)
    oneliner(src_fasta, dst_oneline)


def apply_mutations_single_contig(fasta_in: Path, fasta_out: Path, muts: dict[int, str]) -> int:
    """
    Aplica mutaciones sobre el PRIMER contig del FASTA (típico cloroplasto circular 1-contig).
    Devuelve la longitud final de la secuencia.
    """
    with fasta_in.open("r", encoding="utf-8", errors="replace") as fi, \
         fasta_out.open("w", encoding="utf-8") as fo:
        header = None
        seq_parts = []
        for line in fi:
            if line.startswith(">"):
                header = line.rstrip("\n")
                fo.write(header + "\n")
                break
        for line in fi:
            if line.startswith(">"):
                # ignoramos contigs adicionales en esta versión
                break
            seq_parts.append(line.strip())

        seq = list("".join(seq_parts))
        for pos, base in muts.items():
            idx = pos - 1
            if 0 <= idx < len(seq):
                seq[idx] = base

        fo.write("".join(seq) + "\n")
        return len(seq)


def rename_header_inplace(fasta_path: Path, name_hint: str, species: str = "Hordeum vulgare") -> None:
    """
    Renombra el encabezado del PRIMER contig como:
        >{name_hint}, {species}, {len}
    (len = nº de nucleótidos de la secuencia del primer contig)
    """
    tmp = fasta_path.with_suffix(fasta_path.suffix + ".tmp")
    with fasta_path.open("r", encoding="utf-8", errors="replace") as fi, \
         tmp.open("w", encoding="utf-8") as fo:
        header = None
        seq_parts = []
        for line in fi:
            if line.startswith(">"):
                header = line.rstrip("\n")
                # Dejamos la construcción del nuevo header más adelante (cuando sepamos len)
                break
        for line in fi:
            if line.startswith(">"):
                # sólo primer contig
                break
            seq_parts.append(line.strip())

        seq = "".join(seq_parts)
        new_header = f">{name_hint}, {species}, {len(seq)}"
        fo.write(new_header + "\n")
        fo.write(seq + "\n")

    tmp.replace(fasta_path)


def main():
    ap = argparse.ArgumentParser(
        description="Corrige un FASTA usando un TSV (mpileup), genera oneline y alinea contra referencia numerando la salida."
    )
    ap.add_argument("--fasta", required=True, type=Path, help="FASTA a corregir (1 contig).")
    ap.add_argument("--tsv", required=True, type=Path, help="TSV de recuentos generado por BAMtsv.py.")
    ap.add_argument("--out", type=Path, default=None, help="FASTA corregido (por defecto: *_corrected.fasta).")
    ap.add_argument("--species", type=str, default="Hordeum vulgare", help="Nombre de especie para el header.")
    ap.add_argument("--oneline", action="store_true", help="(Compat) Aceptado pero se genera oneline siempre en un fichero aparte.")
    args = ap.parse_args()

    fasta_in = args.fasta.resolve()
    tsv = args.tsv.resolve()
    fasta_out = (args.out.resolve() if args.out else fasta_in.with_name(fasta_in.stem + "_corrected.fasta"))

    # 1) Corregir
    muts = load_mutations(tsv)
    length = apply_mutations_single_contig(fasta_in, fasta_out, muts)

    # 2) Renombrar cabecera (base limpia + especie + longitud)
    name_hint = normalize_basename(fasta_in.stem)
    rename_header_inplace(fasta_out, name_hint=name_hint, species=args.species)
    print(f"[OK] FASTA corregido: {fasta_out} (len={length})")

    # 3) Generar SIEMPRE un oneline separado
    oneline_out = fasta_out.with_name(fasta_out.stem + "_oneline.fasta")
    write_oneline_copy(fasta_out, oneline_out)
    print(f"[OK] FASTA oneline: {oneline_out}")

    # 4) Alinear contra referencia y numerar salida
    align_with_reference_and_number(oneline_out)


if __name__ == "__main__":
    main()


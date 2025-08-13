#!/usr/bin/env python3

import argparse
import glob
import os
import platform
import re
import shutil
import subprocess
import sys
from typing import Dict, List, Tuple


def fasta_to_oneline_py(src: str, dst: str) -> None:
    """Convierte cualquier FASTA en formato «una línea» (cabecera + secuencia sin saltos)."""
    with open(src, encoding="utf-8", errors="replace") as fi, \
         open(dst, "w", encoding="utf-8") as fo:
        seq_parts: List[str] = []
        header: str | None = None
        for line in fi:
            if line.startswith(">"):
                # Vuelca la secuencia acumulada del contig previo
                if header is not None:
                    fo.write("".join(seq_parts) + "\n")
                # Escribe nueva cabecera y reinicia buffer
                header = line.rstrip("\n")
                fo.write(header + "\n")
                seq_parts.clear()
            else:
                seq_parts.append(line.strip())
        # Último contig
        if seq_parts:
            fo.write("".join(seq_parts) + "\n")

def fasta_to_oneline_awk(src: str, dst: str) -> None:
    AWK_SCRIPT = '/^>/ { if (seq) print seq; print; seq=""; next } { seq = seq $0 } END { print seq }'
    with open(dst, "w", encoding="utf-8") as fo:
        subprocess.run(["awk", AWK_SCRIPT, src], stdout=fo, check=True)

def pick_oneline_backend(force_awk: bool):
    has_awk = shutil.which("awk") is not None
    if force_awk and not has_awk:
        print("--use-awk pedido pero awk no está en el PATH → usando modo Python", file=sys.stderr)
        force_awk = False
    if platform.system() == "Windows":
        return fasta_to_oneline_py
    return fasta_to_oneline_awk if force_awk else fasta_to_oneline_py

def _normalize(name: str) -> str:
    n = name
    n = re.sub(r"Option_\d+_", "", n)
    n = re.sub(r"_Novowrap.*", "", n)
    n = re.sub(r"_pileup", "", n)
    n = re.sub(r"\.fastq(?:\.gz)?_?", "_", n)
    n = re.sub(r"\.(fasta|fa|tsv|txt|fai)$", "", n)
    return n

def pair_fasta_tsv(fasta_dir: str, tsv_dir: str, strict: bool) -> List[Tuple[str, str]]:
    fastas = glob.glob(os.path.join(fasta_dir, "*.fasta"))
    tsvs = glob.glob(os.path.join(tsv_dir, "*.tsv"))

    tsv_map: Dict[str, str] = { _normalize(os.path.basename(t)): t for t in tsvs }
    pairs: List[Tuple[str, str]] = []
    for fa in fastas:
        key = _normalize(os.path.basename(fa))
        if key in tsv_map:
            pairs.append((fa, tsv_map[key]))
            continue
        if strict:
            print(f"No se encontró TSV exacto para {fa}", file=sys.stderr)
            continue
        hits = [p for k, p in tsv_map.items() if k in key or key in k]
        if len(hits) == 1:
            pairs.append((fa, hits[0]))
        elif len(hits) == 0:
            print(f"No se encontró TSV para {fa}", file=sys.stderr)
        else:
            print(f"Múltiples TSV para {fa}: {hits} → ignoro", file=sys.stderr)
    return pairs


def load_mutations(tsv_path: str) -> Dict[int, str]: # Carga mutaciones de un TSV de pileup.

    muts: Dict[int, str] = {}

    with open(tsv_path, encoding="utf-8", errors="replace") as tsv:
        for line in tsv:
            if line.startswith("#") or not line.strip():
                continue  # comentarios / líneas en blanco
            cols = line.rstrip().split("	")

            try:
                pos = int(cols[0])
                offset = 0  # formato 1: POS es columna 0
            except ValueError:
                try:
                    pos = int(cols[1])
                    offset = 1  # formato 2: POS es columna 1
                except (ValueError, IndexError):
                    continue  # línea malformada

            try:
                dot   = int(cols[2 + offset])  # columnas de . y ,
                comma = int(cols[3 + offset])
                counts = list(map(int, cols[4 + offset : 8 + offset]))  # A,C,G,T (4 columnas)
            except (ValueError, IndexError):
                continue  # datos incompletos

            total = dot + comma + sum(counts)
            if total == 0:
                continue

            # Condición de reemplazo: .+, <50 % del total o 0
            if (dot + comma == 0) or ((dot + comma) / total < 0.5):
                bases = "ACGT"  # orden en el TSV típico (A,C,G,T)
                max_idx = max(range(4), key=lambda i: counts[i])
                muts[pos] = bases[max_idx]
    return muts

def apply_mutations(fasta_in: str, fasta_out: str, muts: Dict[int, str]) -> None:
    """Une todas las líneas de secuencia (si fuera necesario) y aplica las mutaciones."""
    with open(fasta_in, encoding="utf-8", errors="replace") as fi:
        header = fi.readline()  # primera línea → cabecera
        seq_parts: List[str] = []
        for line in fi:
            if line.startswith(">"):
                break  # Solo corregimos el primer contig
            seq_parts.append(line.strip())
    seq = list("".join(seq_parts))
    for pos, base in muts.items():
        idx = pos - 1
        if 0 <= idx < len(seq):
            seq[idx] = base
    with open(fasta_out, "w", encoding="utf-8") as fo:
        fo.write(header)
        fo.write("".join(seq) + "\n")



def process_pair(fa_path: str, tsv_path: str, oneline_fn, dirs, fusion_state, species: str):
    base = os.path.splitext(os.path.basename(fa_path))[0]

    name_hint = re.sub(r"^Option_\d+_", "", base)

    name_hint = name_hint.split(".", 1)[0]


    tmp_oneline = os.path.join(dirs['oneline'], f"{base}_oneline.fasta")


    print(f"[1/4] {base}: generando oneline…")
    oneline_fn(fa_path, tmp_oneline)

    muts = load_mutations(tsv_path)
    corrected = os.path.join(dirs['oneline'], f"{base}_corr.fasta")
    print(f"[2/4] {base}: aplicando {len(muts)} mutaciones…")
    apply_mutations(tmp_oneline, corrected, muts)

    # Renombrado de cabeceras
    renamed = os.path.join(dirs['oneline'], f"{base}_corr_renamed.fasta")
    print(f"[3/4] {base}: renombrando cabeceras…")
    rename_headers(corrected, renamed, species, name_hint)

    new_idx = fusion_state['index'] + 1
    fused_out = os.path.join(dirs['fused_out'], f"Fusion_Cloroplastos_{new_idx}.fasta")
    print(f"[4/4] {base}: fusionando → Fusion_{new_idx}.fasta")
    with open(fused_out, "wb") as fo:
        for frag in (fusion_state['file'], renamed):
            with open(frag, "rb") as fi:
                shutil.copyfileobj(fi, fo)
    fusion_state.update(index=new_idx, file=fused_out)

def rename_headers(fasta_in: str, fasta_out: str, species: str, name_hint: str) -> None:
    """Reescribe las cabeceras de un FASTA con el formato:
    >{species}, <longitud_secuencia>
    donde *longitud_secuencia* es el número de caracteres de la secuencia
    hasta la siguiente cabecera (o fin de archivo).
    """
    with open(fasta_in, encoding="utf-8", errors="replace") as fi, \
         open(fasta_out, "w", encoding="utf-8") as fo:
        seq_parts: List[str] = []
        for line in fi:
            if line.startswith(">"):
                if seq_parts:  # vuelca el registro previo
                    seq = "".join(seq_parts)
                    fo.write(f">{name_hint}, {species}, {len(seq)}\n")
                    fo.write(seq + "\n")
                    seq_parts.clear()
                # ignoramos la cabecera original; escribimos la nueva cuando tengamos la longitud
            else:
                seq_parts.append(line.strip())
        # último registro
        if seq_parts:
            seq = "".join(seq_parts)
            fo.write(f">{name_hint}, {species}, {len(seq)}\n")
            fo.write(seq + "\n")

def alignment_has_excess_gaps(aln_path: str, max_gaps_per_line: int = 1000) -> bool:
    """Devuelve True si alguna línea de secuencia tiene >max_gaps_per_line guiones."""
    with open(aln_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line[0] == '>':
                continue
            if line.count('-') > max_gaps_per_line:
                return True
    return False

def run_mafft(in_fasta: str, out_fasta: str):
    print("[MAFFT] Ejecutando MAFFT…")
    with open(out_fasta, "w", encoding="utf-8") as fo:
        subprocess.run(["mafft", in_fasta], stdout=fo, check=True)



def main():
    ap = argparse.ArgumentParser("Fusionador de cloroplastos + corrección TSV")
    ap.add_argument("--fasta", help="Archivo FASTA de cloroplasto")
    ap.add_argument("--tsv",   help="Archivo TSV con mutaciones")
    ap.add_argument("--use-awk", action="store_true", help="Usar awk para convertir FASTA a formato una línea")
    ap.add_argument("--strict-pairing", action="store_true", help="Emparejamiento estricto: solo parejas exactas FASTA↔TSV")
    # Nueva opción para el nombre de la especie
    ap.add_argument("--species", default="Hurdeum vulgare",
                    help="Nombre de la especie para las cabeceras FASTA (por defecto: 'Hurdeum vulgare')")

    args = ap.parse_args()

    #Cogemos de cloroplasto_referencia el fasta y lo copiamos en resultados.
    


    script_dir = os.path.dirname(os.path.realpath(__file__))
    root_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
    temp = os.path.join(root_dir, "temporalDocs")
    ref = os.path.join(root_dir, "cloroplasto_referencia")
    dirs = {
        'selection': os.path.join(temp, "Novowrapselection"),
        'tsv':        os.path.join(temp, "bam"),
        'oneline':    os.path.join(temp, "Oneline"),
        'fused_in':   os.path.join(ref ),
        'fused_out':  os.path.join(root_dir, "Resultados" ),
        'alignment':  os.path.join(temp, "Alineacion"),
        'results':    os.path.join(root_dir, "Resultados"),
    }
    for d in ('oneline', 'fused_out', 'alignment', 'results'):
        os.makedirs(dirs[d], exist_ok=True)

    if args.fasta and args.tsv:
        if args.fasta.endswith(".fai"):
            sys.exit("ERROR: especificaste un índice .fai; pasa el .fasta real")
        pairs = [(args.fasta, args.tsv)]
    else:
        pairs = pair_fasta_tsv(dirs['selection'], dirs['tsv'], args.strict_pairing)
        if not pairs:
            sys.exit("No se encontraron parejas FASTA↔TSV válidas")

    oneline_fn = pick_oneline_backend(args.use_awk)
    # ── comprobar qué fusión existe ya ───────────────────────────────
    patron = os.path.join(dirs['fused_out'], "Fusion_Cloroplastos_*.fasta")
    existentes = glob.glob(patron)
    nums = [int(m.group(1))                                      # números ya usados
            for f in existentes
            for m in [re.search(r"Fusion_Cloroplastos_(\d+)\.fasta$", f)]
            if m]
    last = max(nums) if nums else 1                              # 1 si no hay ninguna

    fusion_state = {
        'index': last,                                           # arrancamos donde lo dejamos
        'file': os.path.join(dirs['fused_out'],
                             f"Fusion_Cloroplastos_{last}.fasta")
    }
    if not os.path.exists(fusion_state['file']):
        sys.exit(f"Archivo de fusión inicial no encontrado: {fusion_state['file']}")

    for fa, tsv in pairs:
        process_pair(fa, tsv, oneline_fn, dirs, fusion_state, args.species)

    aln_temp = os.path.join(dirs['alignment'], f"Fusionado_Cloroplasto_{fusion_state['index']}.aln.fasta")
    run_mafft(fusion_state['file'], aln_temp)

    if alignment_has_excess_gaps(aln_temp, max_gaps_per_line=1000):
        rejected_path = os.path.splitext(aln_temp)[0] + ".rejected.fasta"
        i = 1
        base_rej, ext_rej = os.path.splitext(rejected_path)
        while os.path.exists(rejected_path):
            rejected_path = f"{base_rej}.{i}{ext_rej}"
            i += 1
        os.replace(aln_temp, rejected_path)
        print(f"[DESCARTADO] Alineamiento con >1000 '-' en alguna línea.\n"
              f"Guardado en: {rejected_path}\n"
              f"El archivo NO se moverá a Resultados.")
    else:
        final_aln = os.path.join(dirs['results'], os.path.basename(aln_temp))
        os.replace(aln_temp, final_aln)
        print(f"Resultado final: {final_aln}")

if __name__ == "__main__":
    main()


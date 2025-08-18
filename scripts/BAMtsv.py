#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
import re
from pathlib import Path

# Raíz del repo (.. .. / script.py)
ROOT_DIR = Path(__file__).resolve().parents[1]
LIB_DIR  = ROOT_DIR / "libreris"
MINIMAP2_BIN = LIB_DIR / "minimap2"
SAMTOOLS_BIN = LIB_DIR / "samtools"


def localizar_fasta(ref_path: Path, fq1: Path) -> Path:
    """
    ref_path puede ser un archivo FASTA o un directorio que contiene varios.
    Si es directorio: filtra por extensiones típicas y elige el que contenga el 'base' de fq1;
    si hay varios, usa el más reciente.
    """
    if ref_path.is_file():
        return ref_path.resolve()

    base = fq1.name[:-len("_1.fastq.gz")] if fq1.name.endswith("_1.fastq.gz") else fq1.stem
    exts = (".fa", ".fasta", ".fna")
    candidatos = [p for p in ref_path.iterdir()
                  if p.is_file() and p.suffix.lower() in exts and base in p.name]

    if not candidatos:
        # si no hay coincidencias por base, usa cualquier FASTA del dir (elige el más reciente)
        candidatos = [p for p in ref_path.iterdir() if p.is_file() and p.suffix.lower() in exts]
        if not candidatos:
            raise FileNotFoundError(f"No se encontró FASTA en {ref_path}")
    candidatos = sorted(candidatos, key=lambda p: p.stat().st_mtime, reverse=True)
    return candidatos[0].resolve()


def ensure_faidx(fasta: Path) -> None:
    """Asegura que existe el índice .fai del FASTA (samtools faidx)."""
    fai = fasta.with_suffix(fasta.suffix + ".fai")
    if not fai.exists():
        subprocess.run([str(SAMTOOLS_BIN), "faidx", str(fasta)], check=True)


def align_to_bam(fasta: Path, fq1: Path, fq2: Path, bam_out: Path) -> None:
    """minimap2 (SAM) → samtools sort (BAM) → samtools index (BAI)."""
    bam_out.parent.mkdir(parents=True, exist_ok=True)
    # SAM por stdout
    p1 = subprocess.Popen(
        [str(MINIMAP2_BIN), "-ax", "sr", str(fasta), str(fq1), str(fq2)],
        stdout=subprocess.PIPE
    )
    # sort BAM
    p2 = subprocess.Popen(
        [str(SAMTOOLS_BIN), "sort", "-o", str(bam_out)],
        stdin=p1.stdout
    )
    p1.stdout.close()
    p2.communicate()
    if p2.returncode != 0:
        raise RuntimeError("samtools sort falló")
    # index
    subprocess.check_call([str(SAMTOOLS_BIN), "index", str(bam_out)])


def parse_pileup_line(line: str) -> tuple[str, int, str, dict[str, int]]:
    """
    Parsea una línea de mpileup: CHROM POS REF DEPTH BASES QUALS
    Devuelve (chrom, pos_int, ref_base, counts_dict).
    """
    chrom, pos, refb, depth, bases, _quals = line.rstrip().split("\t", 5)

    counts = {k: 0 for k in ['.', ',', 'A', 'C', 'G', 'T', 'N', '+', '-', '^', '$', 'other']}

    i = 0
    while i < len(bases):
        c = bases[i]

        if c == '^':  # inicio de lectura, siguiente char es MAPQ
            counts['^'] += 1
            i += 2
            continue
        if c == '$':  # fin de lectura
            counts['$'] += 1
            i += 1
            continue
        if c == '.':  # match ref (forward)
            counts['.'] += 1
            i += 1
            continue
        if c == ',':  # match ref (reverse)
            counts[','] += 1
            i += 1
            continue
        if c in '+-':  # indel +<n><seq> o -<n><seq>
            counts[c] += 1
            i += 1
            m = re.match(r'(\d+)', bases[i:])
            if not m:
                counts['other'] += 1
                break
            n = int(m.group(1))
            i += len(m.group(1)) + n
            continue
        cu = c.upper()
        if cu in ('A', 'C', 'G', 'T', 'N'):
            counts[cu] += 1
            i += 1
            continue
        if c == '*':
            # placeholder de deleción en mpileup
            i += 1
            continue

        # IUPAC ambiguas: reparte 1 cuenta entre las bases componentes
        iupac = {
            'R': ('A', 'G'), 'Y': ('C', 'T'), 'S': ('G', 'C'),
            'W': ('A', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
            'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'),
            'H': ('A', 'C', 'T'), 'V': ('A', 'C', 'G')
        }
        if cu in iupac:
            for b in iupac[cu]:
                counts[b] += 1
        else:
            counts['other'] += 1
        i += 1

    return chrom, int(pos), refb, counts


def generar_tsv(bam: Path, fasta: Path, tsv_out: Path,
                min_alt_count: int = 3, min_alt_frac: float = 0.60,
                relax_count: int = 2, relax_frac: float = 0.45) -> None:

    tsv_out.parent.mkdir(parents=True, exist_ok=True)
    ensure_faidx(fasta)

    mp = subprocess.Popen(
        [
            str(SAMTOOLS_BIN), "mpileup",
            "-a",          # incluir todas las posiciones
            "-B",          # sin BAQ (útil en indeles)
            "-Q", "13",    # calidad mínima de base
            "-f", str(fasta),
            str(bam)
        ],
        stdout=subprocess.PIPE,
        text=True
    )
    with tsv_out.open("w", encoding="utf-8") as fo:
        header = ['CHROM', 'POS', '.', ',', 'A', 'C', 'G', 'T', 'N', '+', '-', '^', '$', 'other']
        fo.write("\t".join(header) + "\n")

        for line in mp.stdout:
            if not line.strip():
                continue
            chrom, pos, refb, counts = parse_pileup_line(line)
            alt_sum = counts['A'] + counts['C'] + counts['G'] + counts['T']
            if alt_sum == 0:
                continue

            max_alt = max(counts['A'], counts['C'], counts['G'], counts['T'])
            frac = max_alt / alt_sum

            is_ref_ambig = refb.upper() not in ('A', 'C', 'G', 'T')
            # Criterio:
            # - ref normal: umbral clásico
            # - ref ambigua (N/IUPAC): con que haya >=1 A/C/G/T, incluimos
            if (not is_ref_ambig and (max_alt >= min_alt_count and frac >= min_alt_frac)) \
               or (is_ref_ambig and alt_sum >= 1):
                row = [
                    chrom, str(pos),
                    str(counts['.']), str(counts[',']),
                    str(counts['A']), str(counts['C']), str(counts['G']), str(counts['T']),
                    str(counts['N']), str(counts['+']), str(counts['-']),
                    str(counts['^']), str(counts['$']), str(counts['other'])
                ]
                fo.write("\t".join(row) + "\n")

    mp.wait()
    if mp.returncode not in (0, None):
        raise RuntimeError("samtools mpileup terminó con error")
    

def dellbam(bam_out: Path) -> None:
    try:
        bam_index = bam_out.with_suffix(bam_out.suffix + ".bai")
        bam_out.unlink(missing_ok=True)
        bam_index.unlink(missing_ok=True)
        print(f"[INFO] BAM y BAI eliminados: {bam_out}, {bam_index}")
    except Exception as e:
        print(f"[WARN] No se pudo eliminar BAM o BAI: {e}")


def main():
    ap = argparse.ArgumentParser(
        description="Genera BAM + TSV de recuentos a partir de FASTQ pareados y un FASTA, y lanza alineador.py."
    )
    ap.add_argument("--fq1", required=True, type=Path, help="FASTQ R1 (.fastq.gz)")
    ap.add_argument("--fq2", required=True, type=Path, help="FASTQ R2 (.fastq.gz)")

    # Acepta --fasta y el alias legacy --ref
    ap.add_argument("--fasta", "--ref", dest="fasta", required=True, type=Path,
                    help="FASTA (o directorio con FASTA). --ref también válido.")

    ap.add_argument("--outdir", default=Path(ROOT_DIR / "temporalDocs" / "bam"),
                    type=Path, help="Directorio de salida")
    ap.add_argument("--bam-name", default=None, help="Nombre del BAM (opcional).")
    ap.add_argument("--tsv-name", default=None, help="Nombre del TSV (opcional).")
    ap.add_argument('--species', default="Hordeum vulgare",
                    help="Nombre de la especie usado por alineador.py")

    ap.add_argument("--dellbam", action="store_true", help="Eliminar BAM y BAI al finalizar (opcional).")

    ap.add_argument("--min-alt-count", type=int, default=3)
    ap.add_argument("--min-alt-frac", type=float, default=0.60)
    ap.add_argument("--relax-count", type=int, default=2)
    ap.add_argument("--relax-frac", type=float, default=0.45)

    # Flags para controlar el lanzamiento del alineador
    ap.add_argument("--no-aligner", action="store_true", help="No lanzar alineador.py al finalizar.")
    ap.add_argument("--aligner-path", type=Path, default=None, help="Ruta a alineador.py (por defecto, el vecino de este script).")
    ap.add_argument("--aligner-oneline", action="store_true", help="Añadir --oneline al alineador.")
    ap.add_argument("--aligner-out", type=Path, default=None, help="Ruta del FASTA corregido (por defecto: <ref>_corrected.fasta en --outdir).")

    args = ap.parse_args()

    fq1 = args.fq1.resolve()
    fq2 = args.fq2.resolve()

    # Seleccionar FASTA (acepta archivo o directorio)
    fasta_selected = localizar_fasta(args.fasta.resolve(), fq1)

    outdir = (args.outdir.resolve() if args.outdir else Path.cwd())
    outdir.mkdir(parents=True, exist_ok=True)

    # Nombres de salida por defecto
    bam_out = outdir / (args.bam_name if args.bam_name else (fasta_selected.stem + ".bam"))
    tsv_out = outdir / (args.tsv_name if args.tsv_name else (fasta_selected.stem + ".pileup.tsv"))

    # Alinear y generar BAM
    align_to_bam(fasta_selected, fq1, fq2, bam_out)

    # Generar TSV de pileup
    generar_tsv(
        bam=bam_out,
        fasta=fasta_selected,
        tsv_out=tsv_out,
        min_alt_count=args.min_alt_count,
        min_alt_frac=args.min_alt_frac,
        relax_count=args.relax_count,
        relax_frac=args.relax_frac
    )

    print(f"[OK] BAM: {bam_out}")
    print(f"[OK] TSV: {tsv_out}")

    if args.dellbam:
        dellbam(bam_out)

    # ---- Lanzar alineador.py (opcional) ----
    if not args.no_aligner:
        aligner_path = (args.aligner_path.resolve()
                        if args.aligner_path
                        else Path(__file__).resolve().with_name("alineador.py"))
        if not aligner_path.exists():
            raise FileNotFoundError(f"No encuentro alineador.py en: {aligner_path}")

        aligner_out = (args.aligner_out.resolve()
                       if args.aligner_out
                       else outdir / (fasta_selected.stem + "_corrected.fasta"))

        cmd = [
            sys.executable, str(aligner_path),
            "--fasta", str(fasta_selected),
            "--tsv", str(tsv_out),
            "--out", str(aligner_out),
        ]
        if args.species:
            cmd += ["--species", str(args.species)]
        if args.aligner_oneline:
            cmd += ["--oneline"]

        print(f"[ALIGNER] Ejecutando: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"[OK] FASTA corregido: {aligner_out}")


if __name__ == "__main__":
    main()


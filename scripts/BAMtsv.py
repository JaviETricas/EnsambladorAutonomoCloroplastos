#!/usr/bin/env python3

import argparse                 
import subprocess 
import time              
import sys                      
import re                       
from pathlib import Path        
from typing import List      

# Definimos rutas relativas para que sea portable.
ROOT_DIR = Path(__file__).resolve().parents[1]
MINIMAP2_BIN = ROOT_DIR / 'libreris' / 'minimap2'  
SAMTOOLS_BIN = ROOT_DIR / 'libreris' / 'samtools'

# Función que localiza el FASTA correcto dado un directorio o un archivo.
def localizar_fasta(ref_path: Path, fq1: Path) -> Path:
    if ref_path.is_file():
        return ref_path.resolve()

    base = fq1.name[:-len("_1.fastq.gz")] if fq1.name.endswith("_1.fastq.gz") else fq1.stem
    print(f"[DEBUG] base buscada → {base}", file=sys.stderr)

    candidatos = [p for p in ref_path.iterdir() if p.is_file() and base in p.name]

    if len(candidatos) > 1:                                
        candidatos = sorted(candidatos,                      
                            key=lambda p: p.stat().st_mtime,  
                            reverse=True)                    
        print(f"[INFO] Varios FASTA coinciden; uso el más reciente →"     
              f" {candidatos[0].name}", file=sys.stderr)
 
    fasta_path = candidatos[0].resolve()
    print(f"[DEBUG] FASTA elegido --> {fasta_path}", file=sys.stderr)
    return fasta_path


def alinear_y_sort(fasta: Path, fq1: Path, fq2: Path, bam_out: Path) -> None:
    """Alinea lecturas con minimap2 y genera BAM ordenado + BAI."""
    # minimap2 a SAM (stdout).
    p1 = subprocess.Popen(
        [str(MINIMAP2_BIN), '-ax', 'sr', str(fasta), str(fq1), str(fq2)],
        stdout=subprocess.PIPE
    )
    # samtools sort a BAM ordenado.
    p2 = subprocess.Popen(
        [str(SAMTOOLS_BIN), 'sort', '-o', str(bam_out)],
        stdin=p1.stdout
    )
    p1.stdout.close()
    p2.communicate()
    subprocess.check_call([str(SAMTOOLS_BIN), 'index', str(bam_out)])


# Función auxiliar que cuenta símbolos de la columna 5 de mpileup
def parse_pileup_line(line: str) -> tuple[str, str, dict[str, int]]:
    chrom, pos, refb, depth, bases, _ = line.rstrip().split('\t', 5)

    counts = {k: 0 for k in ['.', ',', 'A', 'C', 'G', 'T', 'N', '+', '-', '^', '$', 'other']}

    i = 0
    while i < len(bases):
        c = bases[i]

        if c == '^':                      # inicio de lectura (lleva 1 char extra = MAPQ)
            counts['^'] += 1
            i += 2
        elif c == '$':                    # fin de lectura
            counts['$'] += 1
            i += 1
        elif c == '.':                    # match a ref (forward)
            counts['.'] += 1
            i += 1
        elif c == ',':                    # match a ref (reverse)
            counts[','] += 1
            i += 1
        elif c in '+-':                   # indel: +<n><seq> o -<n><seq>
            counts[c] += 1
            i += 1
            m = re.match(r'(\d+)', bases[i:])
            if not m:
                counts['other'] += 1
                break
            n = int(m.group(1))
            i += len(m.group(1)) + n
        elif c.upper() in ('A','C','G','T','N'):
            counts[c.upper()] += 1
            i += 1
        else:
            iupac = {
                'R': ('A', 'G'), 'Y': ('C', 'T'), 'S': ('G', 'C'),
                'W': ('A', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
                'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'),
                'H': ('A', 'C', 'T'), 'V': ('A', 'C', 'G')
            }
            u = c.upper()
            if u in iupac:                # reparte ambigua en A/C/G/T
                for b in iupac[u]:
                    counts[b] += 1
            elif c == '*':                # deleción placeholder; ignoramos
                pass
            else:
                counts['other'] += 1
            i += 1

    return chrom, pos, counts


# Función que ejecuta mpileup y escribe un TSV con todos los conteos. 
def generar_tsv(bam: Path, fasta: Path, tsv_out: Path) -> None:
    print(f"[DEBUG] Generando TSV de pileup para {bam} usando {fasta}")

    mpileup = subprocess.Popen(
        [
            str(SAMTOOLS_BIN),
            'mpileup',
            '-Q', '13',     # descarta bases de calidad < 13
            '-f', str(fasta),
            str(bam)
        ],
        stdout=subprocess.PIPE,
        text=True
    )

    MIN_ALT_COUNT = 3       # mín. lecturas de la base mayoritaria A/C/G/T
    MIN_ALT_FRAC  = 0.60    # y que represente al menos el 60% de A/C/G/T en ese sitio

    with tsv_out.open('w') as fo:
        header = ['CHROM', 'POS', '.', ',', 'A', 'C', 'G', 'T', 'N', '+', '-', '^', '$', 'other']
        fo.write('\t'.join(header) + '\n')

        for line in mpileup.stdout:
            chrom, pos, counts = parse_pileup_line(line)

            alt_sum = counts['A'] + counts['C'] + counts['G'] + counts['T']
            if alt_sum == 0:
                # no hay evidencia de A/C/G/T; no hay nada que corregir aquí
                continue

            max_alt = max(counts['A'], counts['C'], counts['G'], counts['T'])
            frac = max_alt / alt_sum

            # Criterio: hay base mayoritaria “suficiente” → incluimos fila
            if (max_alt >= MIN_ALT_COUNT) and (frac >= MIN_ALT_FRAC):
                row = [
                    chrom, pos,
                    str(counts['.']), str(counts[',']),
                    str(counts['A']), str(counts['C']), str(counts['G']), str(counts['T']),
                    str(counts['N']), str(counts['+']), str(counts['-']),
                    str(counts['^']), str(counts['$']), str(counts['other'])
                ]
                fo.write('\t'.join(row) + '\n')

    mpileup.wait()

def load_mutations(tsv_path: str) -> dict[int, str]:
    muts: dict[int, str] = {}

    with open(tsv_path, encoding="utf-8", errors="replace") as tsv:
        for line in tsv:
            if not line.strip() or line.startswith("#") or line.startswith("CHROM"):
                continue

            cols = line.rstrip().split('\t')

            try:
                pos = int(cols[1])        # intenta POS en columna 1
                idx_dot, idx_comma = 2, 3
                idx_A, idx_C, idx_G, idx_T = 4, 5, 6, 7
            except (ValueError, IndexError):
                try:
                    pos = int(cols[0])    # o POS en columna 0
                except (ValueError, IndexError):
                    continue
                idx_dot, idx_comma = 1, 2
                idx_A, idx_C, idx_G, idx_T = 3, 4, 5, 6

            try:
                dot   = int(cols[idx_dot])
                comma = int(cols[idx_comma])
                A = int(cols[idx_A]); C = int(cols[idx_C]); G = int(cols[idx_G]); T = int(cols[idx_T])
            except (ValueError, IndexError):
                continue

            total = dot + comma + A + C + G + T
            if total == 0:
                continue

            # Si hay muy pocos matches a ref o la ref no domina → elegimos la mayoritaria A/C/G/T
            if (dot + comma == 0) or ((dot + comma) / total < 0.5):
                bases = ('A','C','G','T')
                counts = (A, C, G, T)
                muts[pos] = bases[max(range(4), key=lambda i: counts[i])]

    return muts

# Función principal que orchesta todo el proceso.                              
def main(argv: List[str] | None = None) -> None:
    start = time.time()
    parser = argparse.ArgumentParser(
        description="Genera BAM + TSV de pile-up a partir de FASTQ pareados y un FASTA"
    )

    parser.add_argument('--fq1', required=True, type=Path, help="FASTQ R1 (.fastq.gz)")
    parser.add_argument('--fq2', required=True, type=Path, help="FASTQ R2 (.fastq.gz)")
    parser.add_argument('--ref', default=Path(ROOT_DIR / "temporalDocs" / "Novowrapselection"), type=Path, help="Archivo FASTA o carpeta")
    parser.add_argument('--out', default=Path(ROOT_DIR / "temporalDocs" / "bam"), type=Path, help="Directorio de salida")
    parser.add_argument('--species', default="Hurdeum vulgare", help="Nombre de la especie usado por alineador.py")
    args = parser.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    fasta = localizar_fasta(args.ref, args.fq1)
    base = args.fq1.name.replace('_1.fastq.gz', '') if args.fq1.name.endswith('_1.fastq.gz') else args.fq1.stem

    bam_path = args.out / f"{base}.sorted.bam"
    tsv_path = args.out / f"{base}_pileup.tsv"

    alinear_y_sort(fasta, args.fq1, args.fq2, bam_path)
    generar_tsv(bam_path, fasta, tsv_path)

    print(f"Proceso completado.\n   BAM: {bam_path}\n   TSV: {tsv_path}")
    end = time.time()
    print(f"El tiempo transcurrido en este proceso a sido: {end-start:.1f}s\n")

    ensamblador_script = Path(__file__).with_name('alineador.py')
    if not ensamblador_script.exists():
        print(f"[WARN] alineador.py no encontrado en {ensamblador_script.parent}; omitiendo paso de alineador",
              file=sys.stderr)
    else:
        print("Lanzando alineador.py …")
        try:
            subprocess.run([
                sys.executable, str(ensamblador_script),
                '--fasta', str(fasta),
                '--tsv',   str(tsv_path),
                '--species', args.species 
                ], check=True)
        except subprocess.CalledProcessError as exc:
            print(f"[ERROR] alineador.py terminó con código {exc.returncode}", file=sys.stderr)
            sys.exit(exc.returncode)


if __name__ == "__main__":
    main()


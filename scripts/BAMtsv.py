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
MINIMAP2_BIN = ROOT_DIR / 'libreris' / 'minimap2' / 'minimap2'  
SAMTOOLS_BIN = ROOT_DIR / 'libreris' / 'samtools-1.22' / 'samtools'

# Función que localiza el FASTA correcto dado un directorio o un archivo.
def localizar_fasta(ref_path: Path, fq1: Path) -> Path:

    if ref_path.is_file():
        return ref_path.resolve()

    base = fq1.name[:-len("_1.fastq.gz")] if fq1.name.endswith("_1.fastq.gz") else fq1.stem
    print(f"[DEBUG] base buscada → {base}", file=sys.stderr)

    candidatos = [p for p in ref_path.iterdir() if p.is_file() and base in p.name]


    if not candidatos:
        raise FileNotFoundError(
            f"No se encontró un archivo FASTA que contenga ‘{base}’ en {ref_path}"
        )
    if len(candidatos) > 1:
        raise ValueError(
            f"Se encontraron múltiples FASTA que contienen ‘{base}’:\n  " +
            "\n  ".join(map(str, candidatos))
        )

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

    # Inicializa todos los contadores a cero
    counts = {k: 0 for k in ['.', ',', 'A', 'C', 'G', 'T', 'N',
                             '+', '-', '^', '$','other']}

    i = 0
    while i < len(bases):
        c = bases[i]

        if c == '^':                        
            counts['^'] += 1
            i += 2                             
        elif c == '$':                          
            counts['$'] += 1
            i += 1
        elif c in '+-':                          
            counts[c] += 1
            i += 1
            m = re.match(r'(\d+)', bases[i:])
            n = int(m.group(1))
            i += len(m.group(1)) + n             

        elif c.upper() in counts:                
            counts[c.upper()] += 1               
            i += 1
        else:                                    
            iupac = {
                'R': ('A', 'G'), 'Y': ('C', 'T'), 'S': ('G', 'C'),
                'W': ('A', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
                'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'),
                'H': ('A', 'C', 'T'), 'V': ('A', 'C', 'G')
            }
            if c.upper() in iupac:               # reparte en bases reales
                for b in iupac[c.upper()]:
                    counts[b] += 1
            elif c == '*':                       # placeholder de deleción
                pass
            else:
                print(f"[DEBUG] símbolo desconocido → {repr(c)}", file=sys.stderr)
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
        '-Q', '13',           # descarta bases de calidad < 13.
        '-f',  str(fasta),
        str(bam)
    ],
    stdout=subprocess.PIPE,
    text=True
    )

    # Abre el TSV de salida.
    with tsv_out.open('w') as fo:
        header = ['CHROM', 'POS', '.', ',', 'A', 'C', 'G', 'T', 'N',
                  '+', '-', '^', '$', 'other']
        fo.write('\t'.join(header) + '\n')

        # Recorre línea a línea el mpileup
        MIN_VAIANT_READS = 10 # Número mínimo de lecturas para considerar una variante.
        for line in mpileup.stdout:
            chrom, pos, counts = parse_pileup_line(line)

            # Suma de TODOS los símbolos que NO son "." ni "," ni "$" ni "^".
            variant_total = (counts['A'] + counts['C'] + counts['G'] + counts['T'] +
                             counts['N'] + counts['+'] + counts['-'] + counts['other'])

            # Si no hay variantes, se salta la escritura y pasa a la siguiente línea.
            if variant_total < MIN_VAIANT_READS:
                continue  

            # En caso de desajuste
            row = [
                chrom, pos,
                str(counts['.']), str(counts[',']), str(counts['A']), str(counts['C']),
                str(counts['G']), str(counts['T']), str(counts['N']), str(counts['+']),
                str(counts['-']), str(counts['^']), str(counts['$']), str(counts['other'])
            ]
            fo.write('\t'.join(row) + '\n')

    mpileup.wait()


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

    ensamblador_script = Path(__file__).with_name('ensamblador.py')
    if not ensamblador_script.exists():
        print(f"[WARN] ensamblador.py no encontrado en {ensamblador_script.parent}; omitiendo paso de ensamblaje",
              file=sys.stderr)
    else:
        print("Lanzando ensamblador.py …")
        try:
            subprocess.run([
                sys.executable, str(ensamblador_script),
                '--fasta', str(fasta),
                '--tsv',   str(tsv_path)
            ], check=True)
        except subprocess.CalledProcessError as exc:
            print(f"[ERROR] ensamblador.py terminó con código {exc.returncode}", file=sys.stderr)
            sys.exit(exc.returncode)


if __name__ == "__main__":
    main()


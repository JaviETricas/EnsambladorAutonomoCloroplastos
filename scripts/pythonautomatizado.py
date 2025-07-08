#!/usr/bin/env python3


from __future__ import annotations
import os
import subprocess
import time
import sys
import glob
from pathlib import Path


# 1. Generamos rutas relativas para que sea portable.
SCRIPT_DIR = Path(__file__).resolve().parent            # AutomatizerV01/Script
ROOT_DIR   = SCRIPT_DIR.parent                          # AutomatizerV01
LIB_DIR    = ROOT_DIR / "libreris"
TMP_DIR    = ROOT_DIR / "temporalDocs"
RES_DIR    = ROOT_DIR / "resultados"

PAIRS_FILE = SCRIPT_DIR / "parejas.txt"

# Ejecutables y recursos
TALLY_EXE            = LIB_DIR / "tally"
TRIMMOMATIC_DIR      = LIB_DIR / "Trimmomatic-0.39"
TRIMMOMATIC_JAR      = TRIMMOMATIC_DIR / "trimmomatic-0.39.jar"
TRIMMOMATIC_ADAPTERS = TRIMMOMATIC_DIR / "adapters" / "TruSeq3-PE-2.fa"
NOVOWRAP_EXE         = LIB_DIR / "novowrap_env" / "bin" / "novowrap"
CLORO_REF_GB         = ROOT_DIR / "cloroplasto_referencia" / "sequence_cloroplasto.gb"
BLAST_DB             = ROOT_DIR / "cloroplasto_referencia" / "mi_bd"

# Directorios de salida
TALLY_DIR = TMP_DIR / "tally"
TRIM_DIR  = TMP_DIR / "trimmomatic"
SEQ_DIR   = TMP_DIR / "seq_crumbs"  # Reservado por si se activa seqcrumbs
NOV_DIR   = RES_DIR / "novowrap"

for d in (TALLY_DIR, TRIM_DIR, SEQ_DIR, NOV_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Comprobamos las dependencias minimas
for p, desc in [
	(TALLY_EXE, "tally"),
	(TRIMMOMATIC_JAR, "Trimmomatic.jar"),
	(TRIMMOMATIC_ADAPTERS, "Trimmomatic adapters"),
	(NOVOWRAP_EXE, "NOVOwrap"),
	(CLORO_REF_GB, "GenBank referencia"),
]:
    if not p.exists():
        sys.exit(f"[ERROR] No se encontró {desc}: {p}")

#2. Cargar parejas.
if not os.path.isfile(PAIRS_FILE):
	raise FileNotFoundError(f"No existen archivos para cargar")

pairs = []
with open(PAIRS_FILE, 'r') as f:
	for line in f:
		line = line.strip()
		if not line:
			continue
		a, b = line.split(',', 1)
		pairs.append((a, b))
print(f"Cargado {len(pairs)} parejas desde {PAIRS_FILE}")

# Fase 3: Tally Eliminacion de duplicados.
start = time.time()
tally_results: list[Path] = []
for idx, (f1, f2) in enumerate(pairs, start=1):
	name = Path(f1).name  # "C2LWJACXX_5_1_1.fastq.gz"
	base = name[:-len(".fastq.gz")]  # "C2LWJACXX_5_1_1"
	out1 = TALLY_DIR / f"{base}_R1.fastq.gz"
	out2 = TALLY_DIR / f"{base}_R2.fastq.gz"

	print(f"[{idx}/{len(pairs)}] Tally: {base}")
	subprocess.run([
		str(TALLY_EXE),
		'-i', f1, '-j', f2,
		'-o', out1, '-p', out2,
		'--with-quality', '--no-tally', '--pair-by-offset'
	], check=False)
	if os.path.isfile(out1) and os.path.isfile(out2):
		print(f" Se han generado: {out1} {out2}")
		tally_results.extend([out1, out2])
	else:
		print(f"Error en la genracion de {base}")
end = time.time()
print(f"El tiempo transcurrido en este proceso a sido: {end-start:.1f}s\n")

# Fase 4: Trimmomatic Recorte de calidad y eliminar adaptadores.
start = time.time()
trim_paired: list[Path] = []
trim_unpaired: list[Path] = []

for i in range(0, len(tally_results), 2):
	t1, t2 = tally_results[i], tally_results[i+1]
	base = Path(t1).stem.replace("_R1", "" )

	pep1 = TRIM_DIR / f"{base}_pareado_1.fastq.gz"
	upe1 = TRIM_DIR / f"{base}_unpareado_1.fastq.gz"
	pep2 = TRIM_DIR / f"{base}_pareado_2.fastq.gz"
	upe2 = TRIM_DIR / f"{base}_unpareado_2.fastq.gz"
	print(f"Trimmomatic con: {base}")
	cmd = [
		'java', '-jar', str(TRIMMOMATIC_JAR), 'PE', '-threads', '12', '-phred33',
		t1, t2, pep1, upe1, pep2, upe2, f"ILLUMINACLIP:{TRIMMOMATIC_ADAPTERS}:2:30:10", 'LEADING:3', 'TRAILING:3'
	]
	subprocess.run(cmd, check=False)
	if os.path.isfile(pep1) and os.path.isfile(pep2):
		print(f"Documentos generados correctamente: {pep1}, {pep2}")
		trim_paired.extend([pep1, pep2])
	else:
		print(f"Error de pareamiento documentos {base}")
	if os.path.isfile(upe1) and os.path.isfile(upe2):
		print(f"Documentos Unpareados: {upe1}, {upe2}")
		trim_unpaired.extend([upe1, upe2])
	else:
		print(f"Error con los documentos unpareados")

end = time.time()
print(f"Tiempo del proceso Trimmomatic: {end-start:.1f}s\n")


# Fase 5. NOVOWrap Crear consenso frente al cloroplasto.
start = time.time()
novowrap_results: list[Path] = []
for seq_in in range(0, len(trim_paired), 2):
	in1, in2 = trim_paired[seq_in], trim_paired[seq_in + 1]
	base = Path(in1).stem.replace("_paired_R1", "")
	out_dir = NOV_DIR / f"{base}_Novowrap"
	NOV_DIR.mkdir(parents=True, exist_ok=True)

	print(f"[NOVOwrap] {base}")
	cmd = [
		str(NOVOWRAP_EXE),		
		'-input', in1, in2,
		'-ref', str(CLORO_REF_GB),
		'-out', out_dir,
		'-len_diff', '0.8',
		'-debug'
	]
	proc = subprocess.run(
		cmd,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE,
		text=True
		)
	if proc.returncode != 0:
		print(f"NOVOwrap fallo (returncode={proc.returncode})")
		print("stdout:", proc.stdout)
		print("stderr:", proc.stderr)
		continue


	if os.path.isdir(out_dir):
		fasta_files = glob.glob(os.path.join(out_dir, '**', '*.fasta'), recursive=True)
		if fasta_files:
			print(f" Ensamblado completado. Archivos encontrados:")
			for f in fasta_files:
				print (" ", f)
				novowrap_results.append(f)
				# Tras NOVOWrap, si generamos algún ensamblado, borramos los temporales
		else:
			print(f" No se encontro nada dentro de {out_dir}")
	else:
		print(f"No existe el directorio {out_dir}")
end = time.time()
print(f"El tiempo transcurrido en este proceso a sido: {end-start:.1f}s\n")


#Fase 6. Eliminamos ficheros intermedios para evitar la acumulacion de fastq de gran tamaño.
print("[INFO] Eliminando ficheros temporales…")
for tmp in tally_results + trim_unpaired:
    try:
        tmp.unlink(missing_ok=True)
    except Exception as e:
        print(f"   ⚠️  No se pudo borrar {tmp}: {e}")


#Fase 7. Resumen final.
def print_summary():
	print("\nResumen de resultados:")
	print(f" Tally ({len(tally_results)}):")
	for path in tally_results:
		print(f" {path}")
	print(f" Pareados Trimmomatic ({len(trim_paired)}):")
	for path in trim_paired:
		print(f" {path}")
	print(f" Unpareados Trimmomatic ({len(trim_unpaired)}):")
	for path in trim_unpaired:
		print(f" {path}")
	print(f" Novowrap({len(novowrap_results)}):")
	for path in novowrap_results:
		print(f" {path}")
	


# Fase 8. Selección de FASTAs generadas con NOVOwrap

selection_script = Path(__file__).parent / "SeleccionNovowrap.py"
# Nombre del .txt de salida que generará el script de selección
output_list = Path("lista_seleccion_novowrap.txt")

# Construimos el comando:
cmd_select = [
    "python3",
    str(selection_script),
    "-r", str(NOV_DIR),          # raíz de búsqueda: carpeta con todos los *_Novowrap
    "-o", str(output_list)       # fichero de salida
]

print(f"[SELECT] Ejecutando selección de FASTAs: {' '.join(cmd_select)}")
result = subprocess.run(
    cmd_select,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
)

# Comprobamos éxito o fallo
if result.returncode != 0:
    print(f"[SELECT] ERROR (returncode={result.returncode})")
    print(" stdout:", result.stdout)
    print(" stderr:", result.stderr)
else:
    print(f"[SELECT] Selección completada. FASTAs listados en: {output_list}")


print("Pipeline completado.")
print_summary()


#!/usr/bin/env python3.

import os
import subprocess
import time
import glob

PAIRS_FILE = "parejas.txt"
TALLY_DIR = "Dataset_OutPut/tally"
TRIM_DIR = "Dataset_OutPut/Trimmomatic"
SEQ_DIR = 'Dataset_OutPut/seq_crumbs'
NOV_DIR = 'Dataset_OutPut/novowrap'
TRIMMOMATIC_JAR = "/home/jesteban/Jav/DatasetClor/Trimmomatic-0.39/trimmomatic-0.39.jar"
TRIMMOMATIC_ADAPTERS = "/home/jesteban/Jav/DatasetClor/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
BLAST_DB = '/home/jesteban/Jav/DatasetClor/Automatizer/Cloroplasto_Comparacion/mi_bd'
here = os.path.dirname(os.path.abspath(__file__))  
tally_exe = os.path.normpath(os.path.join(here, '..', 'libreris', 'tally')

#Directorios de salida.
os.makedirs(TALLY_DIR, exist_ok=True)
os.makedirs(TRIM_DIR, exist_ok=True)
os.makedirs(NOV_DIR, exist_ok=True)
os.makedirs(SEQ_DIR, exist_ok=True)

#1. Cargar parejas.
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

# Fase 2: Tally Eliminacion de duplicados.
start = time.time()
tally_results = []
for idx, (f1, f2) in enumerate(pairs, start=1):
	base1 = os.path.basename(f1).removesuffix('.fastq.gz')
	base2 = os.path.basename(f2).removesuffix('.fastq.gz')
	out1 = os.path.join(TALLY_DIR, f"{base1}_R1.fastq.gz")
	out2 = os.path.join(TALLY_DIR, f"{base2}_R2.fastq.gz")
	print(f"[{idx}/{len(pairs)}] Tally: {base1}, {base2}")
	subprocess.run([
		'tally_exe',
		'-i', f1, '-j', f2,
		'-o', out1, '-p', out2,
		'--with-quality', '--no-tally', '--pair-by-offset'
	], check=False)
	if os.path.isfile(out1) and os.path.isfile(out2):
		print(f" Se han generado: {out1} {out2}")
		tally_results.extend([out1, out2])
	else:
		print(f"Error en la genracion de {base1}, y {base2}")
end = time.time()
print(f"El tiempo transcurrido en este proceso a sido: {end-start:.1f}s\n")

# Fase 3: Trimmomatic Recorte de calidad y eliminar adaptadores.
start = time.time()
trim_pareados = []
trim_unpareados = []

for i in range(0, len(tally_results), 2):
	t1, t2 = tally_results[i], tally_results[i+1]
	base1 = os.path.basename(t1).removesuffix('_R1.fastq.gz')
	base2 = os.path.basename(t2).removesuffix('_R2.fastq.gz')
	pep1 = os.path.join(TRIM_DIR, f"{base1}_pareado.fastq.gz")
	upe1 = os.path.join(TRIM_DIR, f"{base1}_unpareado.fastq.gz")
	pep2 = os.path.join(TRIM_DIR, f"{base2}_pareado.fastq.gz")
	upe2 = os.path.join(TRIM_DIR, f"{base2}_unpareado.fastq.gz")
	print(f"Trimmomatic con: {base1}, {base2}")
	cmd = [
		'java', '-jar', TRIMMOMATIC_JAR, 'PE', '-threads', '12', '-phred33',
		t1, t2, pep1, upe1, pep2, upe2, f"ILLUMINACLIP:{TRIMMOMATIC_ADAPTERS}:2:30:10", 'LEADING:3', 'TRAILING:3'
	]
	subprocess.run(cmd, check=False)
	if os.path.isfile(pep1) and os.path.isfile(pep2):
		print(f"Documentos generados correctamente: {pep1}, {pep2}")
		trim_pareados.extend([pep1, pep2])
	else:
		print(f"Error de pareamiento documentos {base1}, {base2}")
	if os.path.isfile(upe1) and os.path.isfile(upe2):
		print(f"Documentos Unpareados: {upe1}, {upe2}")
		trim_unpareados.extend([upe1, upe2])
	else:
		print(f"Error con los documentos unpareados")

end = time.time()
print(f"Tiempo del proceso Trimmomatic: {end-start:.1f}s\n")

# Fase 4. Filter_by_blast con seqcrumbs_env.
#Nos saltamos este paso por que NOVOwrap no lo usa.


# Fase 5. NOVOWrap Crear consenso frente al cloroplasto.
start = time.time()
novowrap_results = []
for seq_in in range(0, len(trim_pareados), 2):
	in1, in2 = trim_pareados[i], trim_pareados[i+1]
	base = os.path.basename(in1).removesuffix('_cloro_R1.fastq.gz')
	out_dir = os.path.join(NOV_DIR, f"{base}_Novowrap")
	print(f"\n=== NOVOWrap: {base} ===")
	cmd = [
		'novowrap',
		'-input', in1, in2,
		'-ref', '/home/jesteban/Jav/DatasetClor/Automatizer/Cloroplasto_Comparacion/sequence_cloroplasto.gb',
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
    			print("\nEliminando archivos temporales de Tally y Trimmomatic...")
    			# junta ambas listas
    			temporales = tally_results + trim_pareados + trim_unpareados
    			for tmp in temporales:
    			    try:
    		        	os.remove(tmp)
    			    except OSError:
    				    pass
				print("Archivos temporales eliminados.\n")
		else:
			print(f" No se encontro nada dentro de {out_dir}")
	else:
		print(f"No existe el directorio {out_dir}")
end = time.time()
print(f"El tiempo transcurrido en este proceso a sido: {end-start:.1f}s\n")



#Fase 6. Resumen final.
def print_summary():
	print("\nResumen de resultados:")
	print(f" Tally ({len(tally_results)}):")
	for path in tally_results:
		print(f" {path}")
	print(f" Pareados Trimmomatic ({len(trim_pareados)}):")
	for path in trim_pareados:
		print(f" {path}")
	print(f" Unpareados Trimmomatic ({len(trim_unpareados)}):")
	for path in trim_unpareados:
		print(f" {path}")
	print(f" Novowrap({len(novowrap_results)}):")
	for path in novowrap_results:
		print(f" {path}")
	
	
print("Pipeline completado.")
print_summary()





















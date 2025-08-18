#!/usr/bin/env python3

from __future__ import annotations
import os
import subprocess
import time
import sys
from instaladordependencias import _ensure_python38
from datetime import datetime 
from pathlib import Path

# 1. Generamos rutas relativas para que sea portable.
SCRIPT_DIR = Path(__file__).resolve().parent            # AutomatizerV01/Script
ROOT_DIR   = SCRIPT_DIR.parent                          # AutomatizerV01
LIB_DIR    = ROOT_DIR / "libreris"
TMP_DIR    = ROOT_DIR / "temporalDocs"
RES_DIR    = ROOT_DIR / "resultados"
PAIRS_FILE = SCRIPT_DIR / "parejas.txt"
NOV_SEL    = ROOT_DIR / "temporalDocs" / "Novowrapselection"

# Ejecutables y recursos
TALLY_EXE            = LIB_DIR / "./tally"
TRIMMOMATIC_DIR      = LIB_DIR / "Trimmomatic-0.39"
TRIMMOMATIC_JAR      = TRIMMOMATIC_DIR / "trimmomatic-0.39.jar"
TRIMMOMATIC_ADAPTERS = TRIMMOMATIC_DIR / "adapters" / "TruSeq3-PE-2.fa"
NOVOWRAP_EXE         = LIB_DIR / "novowrap"
CLORO_REF_GB         = ROOT_DIR / "cloroplasto_referencia" / "sequence_cloroplasto.gb"
BLAST_DB             = ROOT_DIR / "cloroplasto_referencia" / "mi_bd"

# Directorios de salida
TALLY_DIR = TMP_DIR / "tally"
TRIM_DIR  = TMP_DIR / "trimmomatic"
NOV_DIR   = TMP_DIR / "novowrap"

for d in (TALLY_DIR, TRIM_DIR, NOV_DIR):
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

# --- después de cargar pairs
tally_results   : list[Path] = []
trim_paired     : list[Path] = []
trim_unpaired   : list[Path] = []
novowrap_results: list[Path] = []

for idx, (f1, f2) in enumerate(pairs, start=1):
    print(f"\n=== PAREJA {idx}/{len(pairs)} ========================= \n")
    start_time = time.time()
    hora_actual = datetime.now()
    print(f"Hora de inicio: {hora_actual.hour:02d}:{hora_actual.minute:02d}")
    
    
    # ── TALLY ─────────────────────────────────────────────
    base = Path(f1).stem[:-2]      # recorta el "_1"
    out1 = TALLY_DIR / f"{base}_1R.fastq.gz"
    out2 = TALLY_DIR / f"{base}_2R.fastq.gz"

    print(f"[TALLY] {base}")
    subprocess.run([
        str(TALLY_EXE), '-i', f1, '-j', f2,
        '-o', out1, '-p', out2,
        '--with-quality', '--no-tally', '--pair-by-offset'
    ], check=False)
    if not (out1.exists() and out2.exists()):
        print("  ⨯ Error en Tally. Se omite esta pareja. \n")
        continue
    tally_results.extend([out1, out2])
    print(f"  ✓ Tally completado en {time.time() - start_time:.1f}s\n")

    
    # ── TRIMMOMATIC ──────────────────────────────────────
    pep1 = TRIM_DIR / f"{base}_pareado_1.fastq.gz"
    upe1 = TRIM_DIR / f"{base}_unpareado_1.fastq.gz"
    pep2 = TRIM_DIR / f"{base}_pareado_2.fastq.gz"
    upe2 = TRIM_DIR / f"{base}_unpareado_2.fastq.gz"

    print(f"[TRIM] {base}")
    cmd_trim = [
        'java', '-jar', str(TRIMMOMATIC_JAR), 'PE', '-threads', '16', '-phred33',
        out1, out2, pep1, upe1, pep2, upe2,
        f"ILLUMINACLIP:{TRIMMOMATIC_ADAPTERS}:2:30:10",
        'LEADING:3', 'TRAILING:3'
    ]
    subprocess.run(cmd_trim, check=False)
    if not (pep1.exists() and pep2.exists()):
        print("  ⨯ Error en Trimmomatic. Se omite esta pareja.")
        continue
    trim_paired.extend([pep1, pep2])
    trim_unpaired.extend([upe1, upe2])
    print(f"  ✓ Trimmomatic completado en {time.time() - start_time:.1f}s\n")


    # ── NOVOWRAP ─────────────────────────────────────────
    out_dir = NOV_DIR / f"{base}_Novowrap"
    print(f"[NOVOwrap] {base}")


    py38 = _ensure_python38()          # Path al intérprete 3.8 (o None si falla)
    if py38 is None:                   # si no tenemos Python 3.8 abandonamos la pareja
        print("  ⨯ No hay Python 3.8 disponible – se omite esta pareja.")
        continue

    cmd_novo = [
        str(py38), "-m", "novowrap",   # ejecuta novowrap con ese intérprete
        "-input", pep1, pep2,
        "-ref",   str(CLORO_REF_GB),
        "-out",   out_dir,
        "-len_diff", "0.8", "-debug"
    ]
    proc = subprocess.run(cmd_novo, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        print("  ⨯ NOVOwrap falló.")
        continue

    fasta_files = list(out_dir.rglob('*.fasta'))
    if fasta_files:
        print("  ✓ Ensamblado completado:")
        for f in fasta_files:
            print("    ", f)
            novowrap_results.append(f)
    else:
        print("  ⨯ No se encontraron FASTA; se omite esta pareja.")
        continue
    print(f"\nPipeline completado en {time.time()-start_time:.1f}s\n")


    # ── LIMPIEZA PARCIAL ─────────────────────────────────
    for tmp in (out1, out2, upe1, upe2):
        tmp.unlink(missing_ok=True)

    # Selección de FASTAs generadas con NOVOwrap
    selection_script = Path(__file__).parent / "SeleccionNovowrap.py"
    # Nombre del .txt de salida que generará el script de selección
    output_list = Path("Errores_de_novowrap.txt")

    # Construimos el comando:
    cmd_select = [
        "python3",
        str(selection_script),
        "-r", str(NOV_DIR),          #si usas NOV_DIR accede a la raíz de búsqueda: carpeta con todos los *_Novowrap
        "-o", str(output_list)       # fichero de salida
    ]

    print(f"[SELECT] Ejecutando selección de FASTAs: {' '.join(cmd_select)}")
    result = subprocess.run(
        cmd_select,
        cwd=str(Path(__file__).parent),  # directorio del script
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    #Fase 6. Eliminamos ficheros intermedios para evitar la acumulacion de fastq de gran tamaño.
    print("[INFO] Eliminando ficheros temporales…")
    for tmp in tally_results + trim_unpaired:
        try:
            tmp.unlink(missing_ok=True)
        except Exception as e:
            print(f" No se pudo borrar {tmp}: {e}")
        
# Comprobamos éxito o fallo
if result.returncode != 0:
    print(f"[SELECT] ERROR (returncode={result.returncode})")
    print(" stdout:", result.stdout)
    print(" stderr:", result.stderr)
else:
    print(f"[SELECT] Selección completada. FASTAs listados en: {output_list}")

print("Pipeline completado.")


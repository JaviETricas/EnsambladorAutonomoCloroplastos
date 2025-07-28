#!/usr/bin/env python3

import os
import csv
import argparse
import shutil
from pathlib import Path


# 1. Generamos rutas relativas para que sea portable.
SCRIPT_DIR = Path(__file__).resolve().parent            # AutomatizerV01/Script
ROOT_DIR   = SCRIPT_DIR.parent                          # AutomatizerV01
TMP_DIR    = ROOT_DIR / "temporalDocs"
NOV_SEL    = TMP_DIR / "Novowrapselection"              # Directorio de selección de NOVOwrap


NOV_DIR   = TMP_DIR / "novowrap"

NOV_SEL.mkdir(parents=True, exist_ok=True)

# Función que comprueba si un valor está dentro de +-tol*ref
def is_within(val, ref, tol=0.1):
    try:
        val = float(val)
        ref = float(ref)
    #Si falla la conversion devolver false.
    except (ValueError, TypeError):
        return False
    return abs(val - ref) <= tol * abs(ref)

# Procesa un CSV y copia el FASTA correspondiente o anota la carpeta si no hay coincidencias
def process_csv(csv_path, output_handle, dest_dir, dest_fail):
    dirpath = Path(csv_path).parent
    found = False

    print(f"Cargando documento .csv: {csv_path}")

    # Lee el CSV usando DictReader
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            # Revisa Length vs Ref_length
            if not is_within(row.get('Length'), row.get('Ref_length')):
                continue

            # Cuenta métricas adicionales dentro de tolerancia
            metrics = ['LSC', 'IRa', 'SSC', 'IRb']
            r_metrics = ['r_LSC', 'r_IRa', 'r_SSC', 'r_IRb']
            count_close = sum(
                1 for m, rn in zip(metrics, r_metrics)
                if is_within(row.get(m), row.get(rn))
            )
            if count_close < 3:
                continue

            # Nombre base del FASTA
            input_name = row.get('Input')
            if not input_name:
                continue

            # Busca archivo .fasta en el mismo directorio
            fasta_file = dirpath / f"{input_name}.fasta"
            if not fasta_file.is_file():
                for fname in dirpath.iterdir():
                    if fname.name.startswith(input_name) and fname.suffix.lower() == '.fasta':
                        fasta_file = fname
                        break

            # Si lo encontramos, copiamos y avisamos
            if fasta_file.is_file():
                dest_dir.mkdir(parents=True, exist_ok=True)
                dest_path = dest_dir / fasta_file.name
                shutil.copy2(fasta_file, dest_path)
                print(f"[COPY] Copiado {fasta_file} → {dest_path}")
                found = True
            break

    # Si no se encontró ningún registro válido, anotamos carpeta y mensaje
    if not found:
        print(f"La carpeta '{dirpath}' no tiene ningún documento válido.")
        os.makedirs(dest_fail, exist_ok=True)              # crea dest_fail si no existe
        try:                                               
            shutil.move(dirpath, dest_fail)                # «cortar y pegar» la carpeta
        except shutil.Error as e:                          # mismo nombre ya existe, etc.
            print(f"[WARN] No se pudo mover '{dirpath}': {e}")  
        output_handle.write(str(dirpath) + '\n')


# Función principal: recorre root_dir y procesa cada CSV
def main(root_dir='.', output_txt='failed_folders.txt'):
    # Directorio donde copiar FASTAs seleccionados
    print (f" root_dir: {root_dir}")
    print (f" output_txt: {output_txt}")
    dest_dir = Path('..') / 'temporalDocs' / 'Novowrapselection'
    dest_fail = Path('..') / 'temporalDocs' / 'NovowrapFail'
    
    # Archivo para anotar carpetas sin coincidencias
    with open(output_txt, 'w') as out_f:
        for dirpath, _, files in os.walk(root_dir):
            for fname in files:
                if fname.lower().endswith('.csv'):
                    process_csv(os.path.join(dirpath, fname), out_f, dest_dir)

    print(f"Procesamiento de selección completado. Carpetas sin coincidencias en: {output_txt}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select FASTA based on CSV criteria')
    parser.add_argument('-r', '--root', default= NOV_DIR, help='Directorio raíz de búsqueda')
    parser.add_argument('-o', '--output', default= SCRIPT_DIR / 'Errores_de_novowrap.txt', help='TXT de carpetas sin coincidencias')
    args = parser.parse_args()
    main(root_dir=args.root, output_txt=args.output)

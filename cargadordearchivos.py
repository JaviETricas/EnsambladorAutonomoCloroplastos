#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
AUTO_DIR = SCRIPT_DIR / "scripts"
PAIRS_FILE = AUTO_DIR / "parejas.txt"
AUTO_SCRIPT = AUTO_DIR / "pythonautomatizado.py"
SELEC_SCRIPT = AUTO_DIR / "SeleccionNovowrap.py"

#Define la funcion para introducir las rutas de los archivos manualmente.
def input_pairs_manual():
    pairs = []

    while True:
        f1 = input('Ruta archivo 1 (.fastq.gz): ').strip()
        f2 = input('Ruta archivo 2 (.fastq.gz): ').strip()
        pairs.append((f1, f2))
        ans = input('¿Añadir otra pareja? [y/n]: ').strip().lower()
        if not ans.startswith('y'):
            break
    return pairs

# Define la fucion para descubrir parejas de archivos en un directorio
def discover_pairs_in_dir():
    directory = input('Ruta de la carpeta con .fastq.gz: ').strip()
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"No existe el directorio: {directory}")

    # files son cada elemento dentro de directory que terminan en .fastq.gz
    files = [f for f in os.listdir(directory) if f.endswith('.fastq.gz')]
    base_map = {}
    for fname in files:
        if fname.endswith('_1.fastq.gz'):
            base = fname[:-len('_1.fastq.gz')]
            # Añade a base como parte 1
            base_map.setdefault(base, {})['1'] = os.path.join(directory, fname)
        elif fname.endswith('_2.fastq.gz'):
            base = fname[:-len('_2.fastq.gz')]
            base_map.setdefault(base, {})['2'] = os.path.join(directory, fname)

    pairs = []
    for base, ends in base_map.items():
        if '1' in ends and '2' in ends:
            pairs.append((ends['1'], ends['2']))
    return pairs

#Escribe un solo par en PAIRS_FILE, sobreescribiendo.
def write_pair_file(pair):
    f1, f2 = pair
    with open(PAIRS_FILE, 'w') as f:
        f.write(f"{f1},{f2}\n")

# Funcion principal que coordina la entrada procesamiento y ejecucion final
def main():
    ans = input('¿Quieres introducir manualmente dos archivos .fasta.gz? [y/n]: ').strip().lower()
    if ans.startswith('y'):
        pairs = input_pairs_manual()
    else:
        pairs = discover_pairs_in_dir()

    # Si no existen parejas: informa y sal de la funcion.
    if not pairs:
        print('No se encontraron parejas. Saliendo.')
        return

    #Escribe todas las parejas procesadas
    print(f"Total parejas a procesar: {len(pairs)}")
    for idx, (f1, f2) in enumerate(pairs, start=1):
        print(f"\nProcesando pareja {idx}/{len(pairs)}:")

        print(f" - {f1}\n - {f2}\n")
        write_pair_file((f1, f2))
        # Ejecuta el script de procesamiento con python3 en el directorio indicado
        ret = subprocess.run(['python3', str(AUTO_SCRIPT)], cwd=str(SCRIPT_DIR))
        if ret.returncode != 0:
            print(f"Error al ejecutar {AUTO_SCRIPT} para esta pareja. Continuando...")
        else:
            print('Procesamiento completado correctamente.')

    print('\nTodas las parejas han sido procesadas.')

    # Ejecuta el script de seleccion en las subcarpetas procesadas.
    print(f"Ejecutando '{SELEC_SCRIPT.name}' en {SELEC_SCRIPT.parent}...")
    sel_ret = subprocess.run(['python3',str(SELEC_SCRIPT)], cwd=str(SELEC_SCRIPT.parent))
    if sel_ret.returncode != 0:
        print(f"Error al ejecutar {SELEC_SCRIPT}. Código de salida: {sel_ret.returncode}")
    else:
        print('SeleccionNovowrap ejecutado correctamente.')

if __name__ == '__main__':
    main()


#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path

SCRIPT_DIR   = Path(__file__).resolve().parent
AUTO_DIR     = SCRIPT_DIR / "scripts"
PAIRS_FILE   = AUTO_DIR / "parejas.txt"
AUTO_SCRIPT  = AUTO_DIR / "pythonautomatizado.py"
SELEC_SCRIPT = AUTO_DIR / "SeleccionNovowrap.py"
NOVOWRAP_DIR = SCRIPT_DIR / 'temporalDocs' / 'Novowrapselection'
BAMTSV       = AUTO_DIR   / 'BAMtsv.py'

#Define la funcion para introducir las rutas de los archivos manualmente.
def input_pairs_manual():
    #Inicia la lista de parejas vacia
    pairs = []
    #Bucle para pedir rutas hasta que el usuario diga que no quiere mas parejas
    while True:
        # Pide la ruta del primer archivo y quita la terminacion fastq.gz
        f1 = input('Ruta archivo 1 (.fastq.gz): ').strip()
        # Pide la ruta del segundo archivo y quita la terminacion fastq.gz
        f2 = input('Ruta archivo 2 (.fastq.gz): ').strip()
        # Añade la pareja de rutas a la lista
        pairs.append((f1, f2))
        # Pregunta si desea añadir otra pareja
        ans = input('¿Añadir otra pareja? [y/n]: ').strip().lower()
        # Si no se presiona y rompe el bucle
        if not ans.startswith('y'):
            break
    # Retorna la lista pairs (parejas)
    return pairs

# Define la fucion para descubrir parejas de archivos en un directorio
def discover_pairs_in_dir():
    """
    Busca en un directorio archivos .fastq.gz que formen parejas
    basadas en sufijos '_1.fastq.gz' y '_2.fastq.gz'.
    """
    #Directory es la entrada que le das
    directory = input('Ruta de la carpeta con .fastq.gz: ').strip()
    # En caso de que directory no sea la ruta de un directorio
    if not os.path.isdir(directory):
        # Sube el error no existe el directorio.
        raise NotADirectoryError(f"No existe el directorio: {directory}")

    # files son cada elemento dentro de directory que terminan en .fastq.gz
    files = [f for f in os.listdir(directory) if f.endswith('.fastq.gz')]
    # Crea un diccionario
    base_map = {}
    # Recorre cada archivo encontrado denominandolo fname
    for fname in files:
        # Si este archivo termina en _1.fastq.gz
        if fname.endswith('_1.fastq.gz'):
            # Base es el nombre de fname menos _1.fastq.gz
            base = fname[:-len('_1.fastq.gz')]
            # Añade a base como parte 1
            base_map.setdefault(base, {})['1'] = os.path.join(directory, fname)
        elif fname.endswith('_2.fastq.gz'):
            base = fname[:-len('_2.fastq.gz')]
            base_map.setdefault(base, {})['2'] = os.path.join(directory, fname)

    # Lista de parejas
    pairs = []
    # Recorre el diccionario para emparejar los documentos.
    for base, ends in base_map.items():
        # Si existen ambos extremos añade la pareja a la lista
        if '1' in ends and '2' in ends:
            pairs.append((ends['1'], ends['2']))
    # Devuelve la lista de parejas.
    return pairs


def write_pair_file(pair):
    """
    Escribe un solo par en PAIRS_FILE, sobreescribiendo.
    """
    # Separa la tupla en 2
    f1, f2 = pair
    # abre el archivo de parejas y sobreescribe su contenido
    with open(PAIRS_FILE, 'w') as f:
        # Escribe las rutas f1 y f2 y salta de linea
        f.write(f"{f1},{f2}\n")

# Funcion principal que coordina la entrada procesamiento y ejecucion final
def main():
    # ans es el input de esa frase
    ans = input('¿Quieres introducir manualmente dos archivos .fasta.gz? [y/n]: ').strip().lower()
    # Si ans empieza por y
    if ans.startswith('y'):
        # Ejecuta la funcion manual
        pairs = input_pairs_manual()
    # Si no, la funcion directorios
    else:
        pairs = discover_pairs_in_dir()

    # Si no existen parejas: informa y sal de la funcion.
    if not pairs:
        print('No se encontraron parejas. Saliendo.')
        return

    #Escribe todas las parejas procesadas
    print(f"Total parejas a procesar: {len(pairs)}")

# lee todas las parejas del fichero 'parejas.txt'
    with open(PAIRS_FILE) as fh:
        pairs = [tuple(line.strip().split(',')) for line in fh if line.strip()]

    if not pairs:
        print("El fichero parejas.txt está vacío.")
        return

    print(f"Total parejas a procesar: {len(pairs)}")

    for idx, (f1, f2) in enumerate(pairs, start=1):
        print(f"\nProcesando pareja {idx}/{len(pairs)}")

        #   Lanza  NovoWrap
        ret = subprocess.run(['python3', str(AUTO_SCRIPT),
                          '--fq1', f1, '--fq2', f2],
                         cwd=str(SCRIPT_DIR))

        if ret.returncode != 0:
            print("FALLÓ pythonautomatizado.py — se omite esta pareja.")
            continue

        # Comprueba que NovoWrap produjo el FASTA
        base = Path(f1).name.replace('_1.fastq.gz', '')
        fasta_found = list(NOVOWRAP_DIR.glob(f'*{base}*.fasta'))
        if not fasta_found:
            print(f"FASTA de {base} no encontrado tras NOVOWRAP — salto pareja.")
            continue

        # Ejecuta BAMtsv
        subprocess.run(
            ['python3', str(BAMTSV),
             '--fq1', f1, '--fq2', f2,
             '--ref', str(NOVOWRAP_DIR)],
            cwd=str(SCRIPT_DIR),
            check=True
        )
        print("✓ pareja procesada con éxito.")

    print('\nTodas las parejas han sido procesadas.')


if __name__ == '__main__':
    main()


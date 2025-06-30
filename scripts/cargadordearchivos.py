#!/usr/bin/env python3


import os
import subprocess

PAIRS_FILE = 'parejas.txt'
AUTO_SCRIPT = 'pythonautomatizado.py'


def input_pairs_manual():
    """
    Pide manualmente rutas de dos archivos y construye la lista de parejas.
    """
    pairs = []
    while True:
        f1 = input('Ruta archivo 1 (.fastq.gz): ').strip()
        f2 = input('Ruta archivo 2 (.fastq.gz): ').strip()
        pairs.append((f1, f2))
        ans = input('¿Añadir otra pareja? [y/n]: ').strip().lower()
        if not ans.startswith('y'):
            break
    return pairs


def discover_pairs_in_dir():
    """
    Busca en un directorio archivos .fastq.gz que formen parejas
    basadas en sufijos '_1.fastq.gz' y '_2.fastq.gz'.
    """
    directory = input('Ruta de la carpeta con .fastq.gz: ').strip()
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"No existe el directorio: {directory}")

    files = [f for f in os.listdir(directory) if f.endswith('.fastq.gz')]
    base_map = {}
    for fname in files:
        if fname.endswith('_1.fastq.gz'):
            base = fname[:-len('_1.fastq.gz')]
            base_map.setdefault(base, {})['1'] = os.path.join(directory, fname)
        elif fname.endswith('_2.fastq.gz'):
            base = fname[:-len('_2.fastq.gz')]
            base_map.setdefault(base, {})['2'] = os.path.join(directory, fname)

    pairs = []
    for base, ends in base_map.items():
        if '1' in ends and '2' in ends:
            pairs.append((ends['1'], ends['2']))
    return pairs


def write_pair_file(pair):
    """
    Escribe un solo par en PAIRS_FILE, sobreescribiendo.
    """
    f1, f2 = pair
    with open(PAIRS_FILE, 'w') as f:
        f.write(f"{f1},{f2}\n")


def main():
    ans = input('¿Quieres introducir manualmente dos archivos .fasta.gz? [y/n]: ').strip().lower()
    if ans.startswith('y'):
        pairs = input_pairs_manual()
    else:
        pairs = discover_pairs_in_dir()

    if not pairs:
        print('No se encontraron parejas. Saliendo.')
        return

    print(f"Total parejas a procesar: {len(pairs)}")
    for idx, (f1, f2) in enumerate(pairs, start=1):
        print(f"\nProcesando pareja {idx}/{len(pairs)}:")
        print(f" - {f1}\n - {f2}\n")
        write_pair_file((f1, f2))
        ret = subprocess.run(['python3', AUTO_SCRIPT])
        if ret.returncode != 0:
            print(f"Error al ejecutar {AUTO_SCRIPT} para esta pareja. Continuando...")
        else:
            print('Procesamiento completado correctamente.')

    print('\nTodas las parejas han sido procesadas.')


if __name__ == '__main__':
    main()

#!/usr/bin/env python 3


import os
import csv
import argparse
import shutil
from pathlib import Path

# Definimos la función is_within con parámetros (val, ref, tol=0.1)
# que comprueba si val está dentro de un tolerancia tol de ref
def is_within(val, ref, tol=0.1):
    try:
        val = float(val)
        ref = float(ref)

    except (ValueError, TypeError):
        return False

    return abs(val - ref) <= tol * abs(ref)

# Definimos la función process_csv con parámetros (csv_path, output_handle)
# que procesa un archivo CSV y escribe resultados en output_handle
def process_csv(csv_path, output_handle):
    dirpath = os.path.dirname(csv_path)
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            if not is_within(row.get('Length'), row.get('Ref_legth')):
                continue

            metrics = ['LSC', 'IRa', 'SSC', 'IRb']
            r_metrics = ['r_LSC', 'r_IRa', 'r_SSC', 'r_IRb']
            count_close = 0
            for m, rn in zip(metrics, r_metrics):
                if is_within(row.get(m), row.get(rm)):
                    #Si no supera esa diferencia se suma 1 a count_close
                    count_close += 1
            if count_close < 3:
                continue

            input_name = row.get('Input')
            if not input_name:
                continue

            #Nos aseguramos que sea fasta.
            fasta_file = os.path.join(dirpath, f"{input_name}.fasta")
            if not os.path.isfile(fasta_file):

                for fname in os.listdir(dirpath):
                    if fname.startswith(input_name) and fname.lower().endswith('.fasta'):
                        fasta_file = os.path.join(dirpath, fname)
                        break


            # Una vez encontremos el .fasta lo añadimos al output
            if os.path.isfile(fasta_file):
                output_handle.write(fasta_file + '\n')

            else:
                output_handle.write(f"# FASTA not found {input_name} in {dirpath}\n")
            found = True
            
            break

            # Si no se encontro ningun reistro valido en el CSV, anotamos la carpeta fallida
    if not found:
        output_handle.write(f" Cloroplastos sin coincidencias tras novowrap: {dirpath}\n")
            

# Definimos la función principal main con parámetros (root_dir, output_txt)
# que recorre el directorio raíz y llama a process_csv para cada CSV
def main(root_dir='.', output_txt='selected_fasta_list.txt'):
    with open(output_txt, 'w') as out_f:
            for dirpath, _, files in os.walk(root_dir):
                for fname in files:
                    # Si es un CSV
                    if fname.lower().endswith('.csv'):
                        csv_path = os.path.join(dirpath, fname)
                        process_csv(csv_path, out_f)

    print(f"Seleccion completada. Resultados en {output_txt}")

    dest_dir = Path("..") / "temporalDocs" / "Novowrapselection"
    dest_dir.mkdir(parents=True, exist_ok=True)

    # abre el fichero de resultados con la lista y copia los seleccionados en la carpeta
    with open(output_txt, 'r') as list_handle:
        for line in list_handle:
            fasta_path = line.strip()
            if not fasta_path or fasta_path.startswith('#'):
                continue
            src = Path(fasta_path)
            if src.is_file():
                dest = dest_dir / src.name
                shutil.copy2(src, dest)
                print(f"[COPY] Copiado {src} → {dest}")
            else:
                print(f"[COPY] ¡No encontrado! {src}")



# Si este script se ejecuta directamente
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Select FASTA files based on CSV criteria')
    parser.add_argument('-r', '--root', default='.', help='Root directory to search')
    parser.add_argument('-o', '--output', default='selected_fasta_list.txt', help='Output TXT file')
    args = parser.parse_args()
    main(root_dir=args.root, output_txt=args.output)




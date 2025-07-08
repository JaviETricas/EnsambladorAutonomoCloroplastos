#!/usr/bin/env python3

import os
import glob
import subprocess
import sys



def main():
    base_dir = '/home/jesteban/Jav/DatasetClor/Automatizer/EnsambladorAutonomo'
    oneline_dir = os.path.join(base_dir, 'Oneline')
    fused_input_dir = os.path.join(base_dir, 'Fusionados_Cloroplastos')
    fused_output_dir = os.path.join(base_dir, 'Cloroplastos_Fusionados')
    alignment_dir = os.path.join(base_dir, 'Alineacion')

    for d in [oneline_dir, fused_output_dir, alignment_dir]:
        os.makedirs(d, exist_ok=True)

    pattern = os.path.join(base_dir, 'Cloroplastos_seleccionados', '*.fasta')
    fasta_files = sorted(glob.glob(pattern))
    if not fasta_files:
        print(f"No se han encontrado archivos .fasta en {base_dir}", file=sys.stderr)
        sys.exit(1)

    #Eliminar saltos de linea.
    awk_script = '/^>/ { if (seq) print seq; print; seq=""; next } { seq = seq $0 } END { print seq }'

    fusion_index = 1
    fusion_file = os.path.join(fused_input_dir, f'Fusion_Cloroplastos_{fusion_index}.fasta')
    if not os.path.exists(fusion_file):
        print(f"Archivo de fusion inicial no encontrado: {fusion_file}", file=sys.stderr)
        sys.exit(1)

    for fasta_path in fasta_files:
        basename = os.path.basename(fasta_path)
        name, ext = os.path.splitext(basename)
        oneline_path = os.path.join(oneline_dir, f"{name}_oneline{ext}")

        print(f"Procesado {basename} -> {os.path.basename(oneline_path)}")

        #Ejecuta AWK
        with open(oneline_path, 'w') as out:
            subprocess.run(['awk', awk_script, fasta_path], stdout=out, check=True)
        
        # Fusionar con el archivo actual
        new_index = fusion_index + 1
        fused_out_path = os.path.join(fused_output_dir, f'Fusion_Cloroplastos_{new_index}.fasta')
        print(f"Fusionando {os.path.basename(fusion_file)} + {os.path.basename(oneline_path)} -> {os.path.basename(fused_out_path)}")
        with open(fused_out_path, 'wb') as out:
            for infile in [fusion_file, oneline_path]:
                with open(infile, 'rb') as f:
                    out.write(f.read())

        fusion_file = fused_out_path
        fusion_index = new_index

    aln_path = os.path.join(alignment_dir, f'Fusionado_Cloroplasto{fusion_index}.aln.fasta')
    print(f"Ejecutando MAFFT: {os.path.basename(fusion_file)} -> {os.path.basename(aln_path)}")
    with open(aln_path, 'w') as out:
        subprocess.run(['mafft', fusion_file], stdout=out, check=True)

    #Abrir con seaview
    print(f"Abrir alineacion en Seaview: {aln_path}")
    subprocess.run(['seaview', aln_path])

if __name__ == '__main__':
    main()

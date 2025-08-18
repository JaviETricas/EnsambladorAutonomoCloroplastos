#!/usr/bin/env python3
import os
import subprocess
import sys
import argparse
import urllib.request
from pathlib import Path

SCRIPT_DIR   = Path(__file__).resolve().parent
AUTO_DIR     = SCRIPT_DIR / "scripts"
PAIRS_FILE   = AUTO_DIR / "parejas.txt"
AUTO_SCRIPT  = AUTO_DIR / "ensambladorcloroplasto.py"
SELEC_SCRIPT = AUTO_DIR / "SeleccionNovowrap.py"
NOVOWRAP_DIR = SCRIPT_DIR / 'temporalDocs' / 'Novowrapselection'
BAMTSV       = AUTO_DIR / 'BAMtsv.py'
TALLY        = AUTO_DIR / 'tally.py'
INTALL_DIR   = AUTO_DIR / 'instaladordependencias.py'

# ---------------- utilidades ----------------

def empty_trash():
    """Intenta vaciar la papelera en Linux/macOS/Windows."""
    try:
        # Linux/FreeDesktop (gio)
        if shutil.which("gio"):
            subprocess.run(["gio", "trash", "--empty"], check=False)
        # KDE (kioclient5)
        elif shutil.which("kioclient5"):
            subprocess.run(["kioclient5", "emptytrash"], check=False)
        # macOS
        elif sys.platform == "darwin":
            subprocess.run(["osascript", "-e", 'tell app "Finder" to empty trash'], check=False)
        # Windows
        elif os.name == "nt":
            import ctypes
            SHERB_NOCONFIRMATION = 0x000001
            SHERB_NOPROGRESSUI   = 0x000002
            SHERB_NOSOUND        = 0x000004
            ctypes.windll.shell32.SHEmptyRecycleBinW(None, None,
                SHERB_NOCONFIRMATION | SHERB_NOPROGRESSUI | SHERB_NOSOUND)
        else:
            # Fallback Linux
            trash = Path.home() / ".local/share/Trash"
            for d in (trash / "files", trash / "info"):
                if d.exists():
                    for p in d.iterdir():
                        try:
                            if p.is_dir():
                                import shutil as _sh
                                _sh.rmtree(p, ignore_errors=True)
                            else:
                                p.unlink(missing_ok=True)
                        except Exception:
                            pass
    except Exception:
        pass

def write_pair_file(pair):
    f1, f2 = pair
    PAIRS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(PAIRS_FILE, 'w') as f:
        f.write(f"{f1},{f2}\n")

# ---------------- entrada interactiva ----------------

def input_pairs_manual():
    """Pide PAREJAS (fq1,fq2) manualmente, con bucle Y/N."""
    pairs = []
    while True:
        f1 = input('Ruta archivo 1 (.fastq.gz): ').strip()
        f2 = input('Ruta archivo 2 (.fastq.gz): ').strip()
        if f1 and f2:
            pairs.append((f1, f2))
        else:
            print("  (Se requiere 1 y 2; no se añadió la pareja)")
        ans = input('¿Añadir otra pareja? [y/n]: ').strip().lower()
        if not ans.startswith('y'):
            break
    return pairs

def discover_pairs_in_dir_path(directory: str):
    """Devuelve parejas (fq1,fq2) encontradas en 'directory'."""
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

def input_dirs_multiple():
    """
    Pide UNA O VARIAS carpetas; tras cada una pregunta si quieres añadir otra.
    Devuelve la lista total de parejas encontradas en todas las carpetas.
    """
    all_pairs = []
    while True:
        d = input('Ruta de carpeta con .fastq.gz: ').strip()
        if not d:
            print("  (Ruta vacía; intenta de nuevo)")
        else:
            try:
                pares = discover_pairs_in_dir_path(d)
                print(f"  Encontradas {len(pares)} parejas en: {d}")
                all_pairs.extend(pares)
            except NotADirectoryError as e:
                print(f"  {e}")
        ans = input('¿Añadir otra carpeta? [y/n]: ').strip().lower()
        if not ans.startswith('y'):
            break
    # eliminar duplicados conservando orden
    seen = set()
    uniq = []
    for p in all_pairs:
        if p not in seen:
            seen.add(p)
            uniq.append(p)
    return uniq


# ---------------- main ----------------
def main():
    import shutil  # usado dentro de empty_trash para el fallback
    parser = argparse.ArgumentParser(description="Pipeline Novowrap+BAMtsv")
    parser.add_argument('--test', action='store_true',
                        help='Descarga un dataset de prueba en ./test y lo procesa')
    parser.add_argument('--species', default="Hurdeum vulgare",
                        help="Nombre de la especie usado por alineador.py")
    parser.add_argument('--dellbam', action='store_true', help='Burra los archivos bam de gran peso al final')
    parser.add_argument('--dell', action='store_true',
                        help='Vaciar la papelera al empezar cada pareja')
    parser.add_argument('--force', action='store_true',
                        help='Pasar --force al instalador de dependencias')
    args = parser.parse_args()

    pairs = []

    # --- modo test (opcional) ---
    if args.test:
        test_dir = SCRIPT_DIR / 'test'
        test_dir.mkdir(exist_ok=True)
        urls = [
            'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR127/030/ERR12745630/ERR12745630_1.fastq.gz',
            'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR127/030/ERR12745630/ERR12745630_2.fastq.gz'
        ]
        local_files = []
        for url in urls:
            dest = test_dir / Path(url).name
            if not dest.is_file():
                print(f"Descargando {dest.name}…")
                urllib.request.urlretrieve(url, dest)
            else:
                print(f"{dest.name} ya existe, se omite la descarga.")
            local_files.append(str(dest))
        pairs = [tuple(local_files)]

    # --- si no hay pares aún, pedir entrada interactiva ---
    if not pairs:
        ans = input('¿Quieres introducir manualmente PAREJAS .fastq.gz? [y/n]: ').strip().lower()
        if ans.startswith('y'):
            pairs = input_pairs_manual()
        else:
            print("Modo: introducir UNA O VARIAS carpetas.")
            pairs = input_dirs_multiple()

    if not pairs:
        print('No se encontraron parejas. Saliendo.')
        return

    print(f"Total parejas encontradas: {len(pairs)}")

    # --- preparar herramientas (tally + instalador) ---
    subprocess.run([sys.executable, str(TALLY)], check=True)

    installer_cmd = [sys.executable, str(INTALL_DIR)]
    if args.force:
        installer_cmd.append('--force')
    subprocess.run(installer_cmd, check=True)

    # --- procesado por parejas ---
    for idx, (f1, f2) in enumerate(pairs, start=1):
        numero = 1
        print(f"\nProcesando pareja {idx}/{len(pairs)}")
        print(f"Pareja {numero}, de {len(pairs)}")
        numero = numero + 1

        if args.dell:
            print("  [dell] Vaciando papelera…")
            empty_trash()

        write_pair_file((f1, f2))

        # 1) NOVOwrap + trimming (ensamblador)
        ret = subprocess.run(
            [sys.executable, str(AUTO_SCRIPT), '--fq1', f1, '--fq2', f2],
            cwd=str(SCRIPT_DIR)
        )
        if ret.returncode != 0:
            print("FALLÓ ensambladorcloroplasto.py — se omite esta pareja.")
            continue

        # 2) comprobar FASTA de NOVOwrap
        base = Path(f1).name.replace('_1.fastq.gz', '')
        fasta_found = list(NOVOWRAP_DIR.glob(f'*{base}*.fasta'))
        if not fasta_found:
            print(f"FASTA de {base} no encontrado tras NOVOWRAP — salto pareja.")
            continue

        # 3) BAM + TSV + alineador
        cmd = [sys.executable, str(BAMTSV),
           '--fq1', f1, '--fq2', f2,
           '--ref', str(NOVOWRAP_DIR),
           '--species', args.species]
        if args.dellbam:
            cmd.append('--dellbam')

        subprocess.run(cmd, cwd=str(SCRIPT_DIR), check=True)
        print("✓ pareja procesada con éxito.")

    print('\nTodas las parejas han sido procesadas.')

if __name__ == '__main__':
    import shutil  # requerido por empty_trash()
    main()

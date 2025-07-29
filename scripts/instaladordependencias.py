#!/usr/bin/env python3
"""
instaladordependencias.py – EnsambladorAutonomoCloroplastos
-----------------------------------------
Instala, enlaza y verifica las herramientas externas necesarias
(minimap2, **Trimmomatic 0.39 (ZIP)**, samtools, novowrap y mafft).

**Rutas portables**: todo vive dentro del repositorio (`libreris/`,
`temporalDocs/`). No exige privilegios de administrador.
**Sin credenciales GitHub**: se descarga por HTTP anónimo o conda.
"""

from __future__ import annotations

# ====================== 1. IMPORTAR LIBRERÍAS =========================
import argparse
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple

# RUTAS RELATIVAS PORTABLES 
SCRIPT_DIR = Path(__file__).resolve().parent            # AutomatizerV01/scripts
ROOT_DIR   = SCRIPT_DIR.parent                          # AutomatizerV01
LIB_DIR    = ROOT_DIR / "libreris"                    # Binarios + wrappers
TMP_DIR    = ROOT_DIR / "libreris"                # Descargas / builds

# Rutas específicas (las necesitan los demás programas)
TRIMMOMATIC_DIR     = LIB_DIR / "Trimmomatic-0.39"
TRIMMOMATIC_JAR     = TRIMMOMATIC_DIR / "trimmomatic-0.39.jar"
TRIMMOMATIC_WRAPPER = LIB_DIR / "trimmomatic"          # lanzador shell

# ====================== 3. CONSTANTES Y CONFIGURACIÓN ================
CONDA_PREFIX = "conda://"   # Marca repos conda
BUILD_PREFIX = "build://"   # Marca instalaciones personalizadas

REPOS: Dict[str, str] = {
    "Trimmomatic": f"{BUILD_PREFIX}trimmomatic", # ZIP oficial
    "samtools":    "https://github.com/samtools/samtools.git",
    "novowrap":    "https://github.com/wpwupingwp/novowrap.git",
    "mafft":       f"{CONDA_PREFIX}mafft",       # Bioconda
    "minimap2":    f"{CONDA_PREFIX}bioconda/label/cf201901::minimap2" 
}

VERIFY_CMDS: Dict[str, List[str]] = {
    "samtools":  ["samtools", "--version"],
    "novowrap":  ["novowrap", "-h"],
    "minimap2":  ["minimap2", "--version"]
}

# Detectar si el entorno permite subprocesos (algunos notebooks web no lo hacen)
PROC_SUPPORTED: bool = sys.platform != "emscripten"


# FUNCIONES AUXILIARES
def run_cmd(cmd: List[str], cwd: Path | None = None, *, accept_nonzero: bool = False) -> Tuple[bool, str]:
    """Ejecuta *cmd* y devuelve (ok, salida)."""
    if not PROC_SUPPORTED:
        return False, "Subprocesos no soportados."
    try:
        res = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, encoding="utf-8", errors="replace")
        salida = (res.stdout + res.stderr).strip()
        ok = res.returncode == 0 or (accept_nonzero and salida)
        return ok, salida or "(sin salida)"
    except (FileNotFoundError, OSError) as e:
        return False, str(e)

# Instalar paquetes conda genéricos (p.ej. mafft)
def conda_install(pkg: str, *, force: bool = False) -> bool:
    if not PROC_SUPPORTED:
        print(f"[WARN] Entorno sin subprocesos: no se puede instalar {pkg}.")
        return False
    if shutil.which(pkg.split("::")[-1]) and not force:
        return True
    print(f"[COND] Instalando {pkg} desde bioconda …")
    ok, out = run_cmd(["conda", "install", "-y", "-c", "bioconda", pkg])
    if not ok:
        print(f"[ERROR] conda no pudo instalar {pkg}:", out)
    return ok


# Instalar Trimmomatic 0.39 desde ZIP oficial
TRIMMOMATIC_ZIP_URL = "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"

def install_trimmomatic_zip(force: bool = False) -> bool:
    if TRIMMOMATIC_WRAPPER.exists() and not force:
        return True
    if not PROC_SUPPORTED:
        print("[WARN] Entorno sin subprocesos: no se puede instalar Trimmomatic.")
        return False

    # Descargar ZIP
    zip_path = TMP_DIR / "Trimmomatic-0.39.zip"
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    print("[DL] Descargando Trimmomatic-0.39.zip …")
    ok, out = run_cmd(["curl", "-L", "-o", str(zip_path), TRIMMOMATIC_ZIP_URL])
    if not ok:
        print("[ERROR] No se pudo descargar Trimmomatic:", out)
        return False

    # Extraer
    extract_dir = TMP_DIR / "trimmomatic_extract"
    shutil.rmtree(extract_dir, ignore_errors=True)
    extract_dir.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(extract_dir)
    except Exception as e:
        print("[ERROR] No se pudo extraer Trimmomatic:", e)
        return False

    src_dir = extract_dir / "Trimmomatic-0.39"
    if not src_dir.exists():
        print("[ERROR] Directorio Trimmomatic-0.39 no encontrado en el ZIP.")
        return False

    # Copiar completo a libreris/
    shutil.rmtree(TRIMMOMATIC_DIR, ignore_errors=True)
    shutil.copytree(src_dir, TRIMMOMATIC_DIR)

    # Permisos ejecución
    if TRIMMOMATIC_JAR.exists():
        TRIMMOMATIC_JAR.chmod(0o755)

    # Crear wrapper
    if TRIMMOMATIC_WRAPPER.exists() or TRIMMOMATIC_WRAPPER.is_symlink():
        TRIMMOMATIC_WRAPPER.unlink()
    with open(TRIMMOMATIC_WRAPPER, "w", encoding="utf-8") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write("java -jar \"$(dirname $0)/Trimmomatic-0.39/trimmomatic-0.39.jar\" \"$@\"\n")
    TRIMMOMATIC_WRAPPER.chmod(0o755)
    print("[OK] Trimmomatic preparado (wrapper + adapters + jar).")
    return True

# Clonar repos Git (samtools, novowrap)
def git_clone(name: str, url: str, *, force: bool = False) -> Path | None:
    if not PROC_SUPPORTED:
        print(f"[WARN] Entorno sin subprocesos: no se puede clonar {name}.")
        return None
    repo_root = TMP_DIR / "repos"
    repo_root.mkdir(parents=True, exist_ok=True)
    dest = repo_root / name
    if dest.exists():
        if force:
            shutil.rmtree(dest)
        else:
            return dest
    print(f"[GIT] Clonando {name} …")
    ok, out = run_cmd(["git", "clone", "--depth", "1", url, str(dest)])
    if not ok:
        print(f"[ERROR] No se pudo clonar {name}:", out)
        return None
    return dest

# Crear enlaces simbólicos en libreris para asegurar que los programas funcionan correctamente
def symlink_if_exists(binary: Path) -> bool:
    if not binary.exists():
        print(f"[WARN] No se encontró binario para {binary.name}; enlace omitido.")
        return False
    target = LIB_DIR / binary.name
    if target.exists() or target.is_symlink():
        target.unlink()
    target.symlink_to(binary.resolve())
    print(f"[OK] Enlace creado para {binary.name}")
    return True

# Verificación de herramientas
def verify(name: str) -> Tuple[bool, str]:
    if not PROC_SUPPORTED:
        return False, "Subprocesos no disponibles."

    # Trimmomatic
    if name == "Trimmomatic":
        if not TRIMMOMATIC_JAR.exists():
            return False, "JAR no encontrado."
        return run_cmd(["java", "-jar", str(TRIMMOMATIC_JAR), "-version"], accept_nonzero=True)

    # mafft
    if name == "mafft":
        exe = shutil.which("mafft")
        if not exe:
            return False, "Ejecutable mafft no en PATH."
        ok, out = run_cmd([exe, "--help"], accept_nonzero=True)
        ok = ok and "MAFFT" in out
        return ok, out if ok else "mafft no respondió."

    # samtools / novowrap
    cmd = VERIFY_CMDS.get(name)
    if cmd is None:
        return False, "Comando de verificación no registrado."
    if not shutil.which(cmd[0]):
        return False, f"{cmd[0]} no en PATH."
    return run_cmd(cmd, accept_nonzero=True)

def main() -> None:
    parser = argparse.ArgumentParser(description="Instala y verifica dependencias externas.")
    parser.add_argument("--force", action="store_true", help="Forzar reinstalación completa")
    args = parser.parse_args()

    installed_now: List[str] = []
    verified_ok: List[str] = []
    errores: List[str] = []

    print("=== Verificando/instalando herramientas requeridas ===")
    for name, locator in REPOS.items():
        print(f"\n→ Procesando {name} …")
        success = False

        if locator.startswith(BUILD_PREFIX):
            if name == "Trimmomatic":
                success = install_trimmomatic_zip(force=args.force)
        elif locator.startswith(CONDA_PREFIX):
            pkg = locator[len(CONDA_PREFIX):]
            success = conda_install(pkg, force=args.force)
        else:
            repo_path = git_clone(name, locator, force=args.force)
            success = repo_path is not None

        if success:
            installed_now.append(name)
        else:
            errores.append(name)
            continue

        if name == "Trimmomatic":
            symlink_if_exists(TRIMMOMATIC_WRAPPER)
        else:
            exe = shutil.which(name.lower())
            if exe:
                symlink_if_exists(Path(exe))

    print("\n=== Verificando funcionamiento de las herramientas ===")
    for name in REPOS:
        ok, msg = verify(name)
        tag = "✔" if ok else "✖"
        print(f"[{tag}] {name}: {msg.splitlines()[0] if msg else msg}")
        (verified_ok if ok else errores).append(name)

    # Deduplicar listas
    verified_ok = sorted(set(verified_ok))
    errores = sorted({e for e in errores if e not in verified_ok})

    print("\n=== Resumen ===")
    if installed_now:
        print("Herramientas instaladas ahora:", ", ".join(installed_now))
    if verified_ok:
        print("Verificadas correctamente:", ", ".join(verified_ok))
    if errores:
        print("Errores/faltantes:", ", ".join(errores))
    if not errores:
        print("Todo instalado y verificado correctamente. ¡Listo! \n")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrumpido por el usuario.")


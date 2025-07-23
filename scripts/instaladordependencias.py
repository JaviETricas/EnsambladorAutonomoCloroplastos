#!/usr/bin/env python3
"""
instaladordependencias.py – AutomatizerV01
-----------------------------------------
Instala, enlaza y verifica las herramientas externas necesarias
(tally, **Trimmomatic 0.39 (ZIP)**, samtools, novowrap y mafft).

Principales reglas de diseño
===========================
1. **Rutas portables**: todo vive dentro del repositorio (`libreris/`,
   `temporalDocs/`). No exige privilegios de administrador.
2. **Sin credenciales GitHub**: se descarga por HTTP anónimo o conda.
3. **Mínimos cambios**: cada revisión solo toca lo imprescindible.

Cambios de esta versión
-----------------------
* **Trimmomatic** se instala desde el ZIP oficial; tras copiar su JAR se le
  da `chmod u+x` y la verificación usa `java -jar … -version`.
* El resto del flujo (tally, samtools, novowrap, mafft) queda intacto.
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

# ====================== 2. RUTAS RELATIVAS PORTABLES ==================
SCRIPT_DIR = Path(__file__).resolve().parent            # AutomatizerV01/scripts
ROOT_DIR   = SCRIPT_DIR.parent                          # AutomatizerV01
LIB_DIR    = ROOT_DIR / "libreris"                     # Binarios + enlaces
TMP_DIR    = ROOT_DIR / "temporalDocs"                 # Descargas / builds
RES_DIR    = ROOT_DIR / "resultados"                   # Reservado para otros módulos

# Rutas específicas que otros módulos consumen
TALLY_EXE           = LIB_DIR / "tally"
TRIMMOMATIC_DIR     = LIB_DIR / "Trimmomatic-0.39"
TRIMMOMATIC_JAR     = TRIMMOMATIC_DIR / "trimmomatic-0.39.jar"
TRIMMOMATIC_WRAPPER = LIB_DIR / "trimmomatic"          # lanzador shell
NOVOWRAP_EXE        = LIB_DIR / "novowrap_env" / "bin" / "novowrap"

# ====================== 3. CONSTANTES Y CONFIGURACIÓN ================
CONDA_PREFIX = "conda://"   # Marca para repos conda
BUILD_PREFIX = "build://"   # Marca para instalación a medida

REPOS: Dict[str, str] = {
    "tally":       f"{CONDA_PREFIX}reaper",              # Compilación de reaper
    "Trimmomatic": f"{BUILD_PREFIX}trimmomatic_zip",   # ZIP oficial
    "samtools":    "https://github.com/samtools/samtools.git",
    "novowrap":    "https://github.com/wpwupingwp/novowrap.git",
    "mafft":       f"{CONDA_PREFIX}mafft",              # Bioconda
}

# Comandos de verificación rápidos
VERIFY_CMDS: Dict[str, List[str]] = {
    "samtools":  ["samtools", "--version"],
    "novowrap":  ["novowrap", "-h"],
}

# Detectar entornos sin subprocesos (ej. emscripten)
PROC_SUPPORTED: bool = sys.platform != "emscripten"

# =====================================================================
# 4. FUNCIONES AUXILIARES
# =====================================================================

def run_cmd(cmd: List[str], cwd: Path | None = None, *, accept_nonzero: bool = False) -> Tuple[bool, str]:
    """Ejecuta *cmd* y devuelve `(ok, salida)`."""
    if not PROC_SUPPORTED:
        return False, "Subprocesos no soportados en este entorno."
    try:
        res = subprocess.run(
            cmd,
            cwd=cwd,
            text=True,
            capture_output=True,
            encoding="utf-8",
            errors="replace",
        )
        salida = (res.stdout + res.stderr).strip()
        ok = res.returncode == 0 or (accept_nonzero and salida)
        return ok, salida or "(sin salida)"
    except (FileNotFoundError, OSError) as e:
        return False, str(e)
    except OSError as e:
        return False, str(e)


# ---------------------------------------------------------------------
# 4.1  tally  → instalar vía conda (reaper) y copiar binario
# ---------------------------------------------------------------------

def install_tally(force: bool = False) -> bool:
    """Instala `tally` mediante el paquete *reaper* de Bioconda y lo copia a
    `libreris/tally` para que el resto de programas lo encuentren allí.
    """
    if TALLY_EXE.exists() and not force:
        return True  # ✅ Ya disponible

    if not PROC_SUPPORTED:
        print("[WARN] Entorno sin subprocesos: no se puede instalar tally.")
        return False

    # Instalar reaper (incluye tally) con conda
    print("[COND] Instalando reaper (contiene tally) desde bioconda …")
    ok, out = run_cmd(["conda", "install", "-y", "-c", "bioconda", "reaper"])
    if not ok:
        print("[ERROR] conda no pudo instalar reaper:", out)
        return False

    # Localizar binario tally instalado por conda
    tally_path = shutil.which("tally")
    if tally_path is None:
        print("[ERROR] No se encontró el ejecutable tally en PATH tras la instalación.")
        return False

    # Copiar a libreris/
    LIB_DIR.mkdir(parents=True, exist_ok=True)
    if TALLY_EXE.exists() or TALLY_EXE.is_symlink():
        TALLY_EXE.unlink()
    shutil.copy2(tally_path, TALLY_EXE)
    TALLY_EXE.chmod(0o755)
    print("[OK] tally instalado y copiado a libreris/.")
    return True

# ------------------------------------------------------------------
# 4.2  Instalar mafft vía conda (sin cambios)
# ------------------------------------------------------------------

def conda_install(pkg: str, *, force: bool = False) -> bool:
    if not PROC_SUPPORTED:
        print(f"[WARN] Entorno sin subprocesos: no se puede instalar {pkg}.")
        return False
    if shutil.which(pkg.replace("::", "")) and not force:
        return True
    print(f"[COND] Instalando {pkg} desde bioconda …")
    ok, out = run_cmd(["conda", "install", "-y", "-c", "bioconda", pkg])
    if not ok:
        print(f"[ERROR] conda no pudo instalar {pkg}:", out)
    return ok


# ---------------------------------------------------------------------
# 4.3  Instalar Trimmomatic desde ZIP oficial
# ---------------------------------------------------------------------

TRIMMOMATIC_ZIP_URL = "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"

def install_trimmomatic_zip(force: bool = False) -> bool:
    if TRIMMOMATIC_WRAPPER.exists() and not force:
        return True
    if not PROC_SUPPORTED:
        print("[WARN] Entorno sin subprocesos: no se puede instalar Trimmomatic.")
        return False

    zip_path = TMP_DIR / "Trimmomatic-0.39.zip"
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    print("[DL] Descargando Trimmomatic-0.39.zip …")
    ok, out = run_cmd(["curl", "-L", "-o", str(zip_path), TRIMMOMATIC_ZIP_URL])
    if not ok:
        print("[ERROR] No se pudo descargar Trimmomatic:", out)
        return False

    extract_dir = TMP_DIR / "trimmomatic_extract"
    shutil.rmtree(extract_dir, ignore_errors=True)
    extract_dir.mkdir(parents=True, exist_ok=True)
    ok, out = run_cmd(["unzip", "-q", str(zip_path), "-d", str(extract_dir)])
    if not ok:
        print("[ERROR] No se pudo extraer Trimmomatic:", out)
        return False

    src_dir = extract_dir / "Trimmomatic-0.39"
    if not src_dir.exists():
        print("[ERROR] Directorio Trimmomatic-0.39 no encontrado en el ZIP.")
        return False

    shutil.rmtree(TRIMMOMATIC_DIR, ignore_errors=True)
    shutil.copytree(src_dir, TRIMMOMATIC_DIR)

    if TRIMMOMATIC_JAR.exists():
        TRIMMOMATIC_JAR.chmod(0o755)  # permiso ejecución

        
# ------------------------------------------------------------------
# 4.4  Clonar repos git (samtools, novowrap)
# ------------------------------------------------------------------

def git_clone(name: str, url: str, *, force: bool = False) -> Path | None:
    """Clona un repo en `TMP_DIR/repos/<name>` y devuelve la ruta local."""
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

# ------------------------------------------------------------------
# 4.5  Enlaces simbólicos en libreris/
# ------------------------------------------------------------------

def symlink_if_exists(binary: Path) -> bool:
    """Crea/actualiza un enlace en LIB_DIR apuntando a *binary*."""
    if not binary.exists():
        print(f"[WARN] No se encontró binario para {binary.name}; enlace omitido.")
        return False
    target = LIB_DIR / binary.name
    if target.exists() or target.is_symlink():
        target.unlink()
    target.symlink_to(binary.resolve())
    print(f"[OK] Enlace creado para {binary.name}")
    return True

# ------------------------------------------------------------------
# 4.6  Verificación de cada herramienta
# ------------------------------------------------------------------

def verify(name: str) -> Tuple[bool, str]:
    """Ejecuta el comando de ayuda/versión y devuelve `(ok, output)`."""
    if not PROC_SUPPORTED:
        return False, "Subprocesos no disponibles."

    # --- tally ---
    if name == "tally":
        if not TALLY_EXE.exists():
            return False, "tally no presente."
        return run_cmd([str(TALLY_EXE), "-h"], accept_nonzero=True)

    # --- Trimmomatic ---
    if name == "Trimmomatic":
        if not TRIMMOMATIC_JAR.exists():
            return False, "JAR no encontrado."
        return run_cmd(["java", "-jar", str(TRIMMOMATIC_JAR), "-version"], accept_nonzero=True)

    # --- mafft ---
    if name == "mafft":
        exe = shutil.which("mafft")
        if not exe:
            return False, "Ejecutable mafft no en PATH."
        ok, out = run_cmd([exe, "--help"], accept_nonzero=True)
        ok = ok and "MAFFT" in out
        return ok, out if ok else "mafft no respondió correctamente."

    # --- samtools / novowrap ---
    cmd = VERIFY_CMDS.get(name)
    if not cmd:
        return False, "Sin comando de verificación registrado."
    if not shutil.which(cmd[0]):
        return False, f"{cmd[0]} no en PATH."
    return run_cmd(cmd, accept_nonzero=True)

# =====================================================================
# 5. FLUJO PRINCIPAL
# =====================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Instala y verifica herramientas externas.")
    parser.add_argument("--force", action="store_true", help="Forzar reinstalación completa")
    args = parser.parse_args()

    installed_now: List[str] = []
    verified_ok: List[str] = []
    errors: List[str] = []

    print("=== Verificando/instalando herramientas requeridas ===")
    for name, locator in REPOS.items():
        print(f"\n→ Procesando {name} …")
        success = False

        # ------------------------------ tally ----------------------
        if locator == f"{BUILD_PREFIX}tally":
            success = install_tally(force=args.force)

        # ------------------------------ Trimmomatic ZIP -----------
        elif locator == f"{BUILD_PREFIX}trimmomatic_zip":
            success = install_trimmomatic_zip(force=args.force)

        # ------------------------------ conda packages ------------
        elif locator.startswith(CONDA_PREFIX):
            pkg = locator[len(CONDA_PREFIX):]
            success = conda_install(pkg, force=args.force)

        # ------------------------------ git repos -----------------
        else:
            repo_path = git_clone(name, locator, force=args.force)
            success = repo_path is not None

        if success:
            installed_now.append(name)
        else:
            errors.append(name)
            continue  # No intentar enlazar ni verificar

        # --------------- crear enlaces simbólicos -----------------
        if name == "tally":
            # el binario ya está en libreris; opcional crear enlace duplicado
            pass
        elif name == "Trimmomatic":
            symlink_if_exists(TRIMMOMATIC_WRAPPER)
        elif name == "samtools":
            exe = shutil.which("samtools")
            if exe:
                symlink_if_exists(Path(exe))
        elif name == "novowrap":
            exe = shutil.which("novowrap")
            if exe:
                symlink_if_exists(Path(exe))
        elif name == "mafft":
            exe = shutil.which("mafft")
            if exe:
                symlink_if_exists(Path(exe))

    # =================================================================
    # Verificación final
    # =================================================================
    print("\n=== Verificando funcionamiento de las herramientas ===")
    for name in REPOS:
        ok, msg = verify(name)
        tag = "✔" if ok else "✖"
        print(f"[{tag}] {name}: {msg.splitlines()[0]}")
        if ok:
            verified_ok.append(name)
        else:
            errors.append(name)

    # =================================================================
    # Resumen
    # =================================================================
    print("\n=== Resumen ===")
    if installed_now:
        print("Herramientas instaladas ahora:", ", ".join(installed_now))
    if verified_ok:
        print("Verificadas correctamente:", ", ".join(sorted(set(verified_ok))))
    if errors:
        print("Errores/faltantes:", ", ".join(sorted(set(errors))))
    else:
        print("Todo instalado y verificado correctamente. ¡Listo! ✅")

# =====================================================================
# 6. PUNTO DE ENTRADA
# =====================================================================

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrumpido por el usuario.")

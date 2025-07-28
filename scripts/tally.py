#!/usr/bin/env python3
"""
install_tally.py
Descarga/instala Reaper (tally) con Conda-Bioconda, localiza el binario y lo
copia a ../libreris/tally para que tu proyecto lo use de forma portable.
"""

import os
import subprocess
import shutil
import sys


# Configuración
CONDA_ENV  = "bioenv"    # Nombre del entorno conda donde instalar reaper
PACKAGE    = "reaper"    # El paquete que incluye tally
TARGET_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "libreris"))
TARGET_BIN = os.path.join(TARGET_DIR, "tally")

# Funciones auxiliares
def run(cmd, **kw):
    """Ejecuta un comando y muestra stdout/stderr en tiempo real."""
    print(f"[CMD] {' '.join(cmd)}")
    subprocess.run(cmd, check=True, **kw)

# funcion para verificar si un entorno conda existe
def conda_env_exists(env):
    result = subprocess.run(["conda", "env", "list", "--json"], check=True, capture_output=True, text=True)
    import json
    envs = json.loads(result.stdout)["envs"]
    return any(env.endswith("/" + env_name) or env == env_name for env_name in envs)

# funcion para instalar un paquete en un entorno conda
def conda_install(env, package):
    run(["conda", "install", "-y", "-n", env, "-c", "bioconda", "-c", "conda-forge", package])

# funcion para crear un entorno conda y instalar un paquete
def conda_create(env, package):
    run(["conda", "create", "-y", "-n", env, "-c", "bioconda", "-c", "conda-forge", package])

# funcion para ejecutar un comando dentro de un entorno conda
def conda_run(env, *cmd): 
    run(["conda", "run", "-n", env, *cmd]) 

# funcion para localizar un binario dentro de un entorno conda
def which_in_env(env, binary):
    result = subprocess.run(
        ["conda", "run", "-n", env, "bash", "-c", f"command -v {binary}"],
        check=True, capture_output=True, text=True
    )
    path = result.stdout.strip()
    if not path:
        raise RuntimeError(f"No se encontró {binary} en el entorno {env}")
    return path

# ---------------------------------------------------------------------------
# 1. Crear o actualizar el entorno con reaper
# ---------------------------------------------------------------------------
print(f"[INFO] Verificando entorno '{CONDA_ENV}' …")
if conda_env_exists(CONDA_ENV):
    print(f"[INFO] Entorno '{CONDA_ENV}' ya existe - instalando/actualizando {PACKAGE}")
    conda_install(CONDA_ENV, PACKAGE)
else:
    print(f"[INFO] Entorno '{CONDA_ENV}' no existe - creándolo con {PACKAGE}")
    conda_create(CONDA_ENV, PACKAGE)

# ---------------------------------------------------------------------------
# 2. Localizar tally dentro del entorno
# ---------------------------------------------------------------------------
print("[INFO] Localizando 'tally' en el entorno …")
tally_path = which_in_env(CONDA_ENV, "tally")
print(f"[INFO] tally encontrado en: {tally_path}")

# ---------------------------------------------------------------------------
# 3. Copiar a ../libreris/tally
# ---------------------------------------------------------------------------
os.makedirs(TARGET_DIR, exist_ok=True)
if os.path.exists(TARGET_BIN):
    os.remove(TARGET_BIN)           # Evita bucles de symlinks
shutil.copy2(tally_path, TARGET_BIN)
os.chmod(TARGET_BIN, 0o755)
print(f"[INFO] Copiado tally a {TARGET_BIN}")




#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
import argparse, json, os, shutil, subprocess, sys, zipfile
from pathlib import Path
from typing import Dict, List, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
LIB_DIR    = ROOT_DIR / "libreris"        
TMP_DIR    = ROOT_DIR / "temporalDocs"     

TRIMMOMATIC_DIR     = LIB_DIR / "Trimmomatic-0.39"
TRIMMOMATIC_JAR     = TRIMMOMATIC_DIR / "trimmomatic-0.39.jar"
TRIMMOMATIC_WRAPPER = LIB_DIR / "trimmomatic"

CONDA_PREFIX = "conda://"
BUILD_PREFIX = "build://"

REPOS: Dict[str, str] = {
    "Trimmomatic": f"{BUILD_PREFIX}trimmomatic",
    "novowrap":    f"{BUILD_PREFIX}novowrap",       # pip + Py3.8
    "minimap2":    f"{BUILD_PREFIX}minimap2",       # git + make
    "samtools":    f"{BUILD_PREFIX}samtools",       # **cambia: env dedicado**
    "mafft":       f"{CONDA_PREFIX}bioconda::mafft",
}

# Comandos de verificación rápidos
VERIFY_CMDS: Dict[str, List[str]] = {
    "samtools": ["samtools", "--version"],
    "mafft":    ["mafft", "--help"],
    "minimap2":    ["minimap2", "--version"]
}

PROC_SUPPORTED = sys.platform != "emscripten"


def run_cmd(cmd: List[str], *, cwd: Path | None = None,
            accept_nonzero: bool = False) -> Tuple[bool, str]:
    """Ejecuta cmd y devuelve (ok, salida)."""
    if not PROC_SUPPORTED:
        return False, "Subprocesos no soportados."
    try:
        res = subprocess.run(
            cmd, cwd=cwd, text=True, capture_output=True,
            encoding="utf-8", errors="replace"
        )
        out = (res.stdout + res.stderr).strip()
        ok  = res.returncode == 0 or (accept_nonzero and out)
        return ok, out or "(sin salida)"
    except (FileNotFoundError, OSError) as e:
        return False, str(e)


#  Instalador Trimmomatic (ZIP oficial)
TRIM_ZIP_URL = ("http://www.usadellab.org/cms/uploads/supplementary/"
                "Trimmomatic/Trimmomatic-0.39.zip")

def install_trimmomatic_zip(force: bool = False) -> bool:
    if TRIMMOMATIC_WRAPPER.exists() and not force:
        return True
    if not PROC_SUPPORTED:
        print("[WARN] Sin subprocesos: Trimmomatic omitido.")
        return False
    try:
        import urllib.request
        TMP_DIR.mkdir(parents=True, exist_ok=True)
        zip_path = TMP_DIR / "Trimmomatic-0.39.zip"
        print("[DL] Descargando Trimmomatic-0.39.zip …")
        urllib.request.urlretrieve(TRIM_ZIP_URL, zip_path)
    except Exception as e:
        print("[ERROR] Descarga fallida:", e)
        return False

    extract_dir = TMP_DIR / "trim_extract"
    shutil.rmtree(extract_dir, ignore_errors=True)
    extract_dir.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(extract_dir)
    except Exception as e:
        print("[ERROR] No se pudo extraer:", e); return False

    src = extract_dir / "Trimmomatic-0.39"
    if not src.exists():
        print("[ERROR] Estructura ZIP inesperada."); return False
    shutil.rmtree(TRIMMOMATIC_DIR, ignore_errors=True)
    shutil.copytree(src, TRIMMOMATIC_DIR)

    if TRIMMOMATIC_JAR.exists():
        TRIMMOMATIC_JAR.chmod(0o755)

    if TRIMMOMATIC_WRAPPER.exists() or TRIMMOMATIC_WRAPPER.is_symlink():
        TRIMMOMATIC_WRAPPER.unlink()

    with open(TRIMMOMATIC_WRAPPER, "w", encoding="utf-8") as fh:
        fh.write("#!/usr/bin/env bash")
        fh.write("java -jar \"$(dirname $0)/Trimmomatic-0.39/trimmomatic-0.39.jar\" \"$@\"")
    TRIMMOMATIC_WRAPPER.chmod(0o755)

    print("[OK] Trimmomatic preparado.")
    return True


# Instalar minimap2
def install_minimap2_source(force: bool = False) -> bool:
    dest_src = LIB_DIR / "minimap2_src"
    bin_path = LIB_DIR / "minimap2"
    if bin_path.exists() and not force:
        return True
    if not PROC_SUPPORTED:
        print("[WARN] Sin subprocesos: minimap2 omitido."); return False

    if force and dest_src.exists():
        shutil.rmtree(dest_src)

    print("[GIT] Clonando minimap2 …")
    ok, out = run_cmd(["git", "clone", "--depth", "1",
                       "https://github.com/lh3/minimap2", str(dest_src)])
    if not ok:
        print("[ERROR] clon minimap2:", out); return False

    print("[MAKE] Compilando minimap2 …")
    ok, out = run_cmd(["make"], cwd=dest_src)
    if not ok:
        print("[ERROR] make minimap2:", out); return False

    built_bin = dest_src / "minimap2"
    if not built_bin.exists():
        print("[ERROR] Binario minimap2 no generado."); return False
    try:
        if bin_path.is_symlink() or bin_path.exists():
            bin_path.unlink() 
    except OSError:
        bin_path.rename(bin_path.with_suffix(".old"))

    shutil.copy2(built_bin, bin_path)
    bin_path.chmod(0o755)
    print("[OK] minimap2 compilado y disponible.")
    return True                                           


# Instalar mafft con conda
def conda_install(pkg: str, *, force: bool = False) -> bool:
    if not PROC_SUPPORTED:
        print(f"[WARN] Sin subprocesos: {pkg} omitido."); return False
    base_name = pkg.split("::")[-1]
    if shutil.which(base_name) and not force:
        return True

    cmd = ["conda", "install", "-y"]
    if "::" in pkg:
        canal, nombre = pkg.split("::", 1)
        cmd.extend(["-c", canal, pkg])
    else:
        canal, nombre = "bioconda", pkg
        cmd.extend(["-c", canal])

    if nombre.startswith("samtools") and "python=3.11" not in cmd:
        cmd.append("python=3.11")

    cmd.append(pkg) 

    print(f"[COND] Instalando {pkg} …")
    ok, out = run_cmd(cmd)
    if not ok:
        print("[ERROR] conda install:", out)
    return ok


# Instalar samtools
def install_samtools_env(force: bool = False) -> bool:
    env_name = "samtools_env"
    conda_prefix = Path(os.environ.get("CONDA_PREFIX", Path.home() / "miniconda3"))
    env_dir = conda_prefix / "envs" / env_name
    samtools_exe = env_dir / "bin" / "samtools"

    if samtools_exe.exists() and not force:
        symlink_if_exists("samtools")
        return True

    if force and env_dir.exists():
        print(f"[INFO] Eliminando entorno antiguo '{env_name}' …")
        shutil.rmtree(env_dir)

    if not env_dir.exists():
        print(f"[INFO] Creando entorno '{env_name}' con python=3.11 …")
        ok, out = run_cmd([
            "conda", "create", "-y", "-n", env_name,
            "python=3.11", "-c", "conda-forge"
        ])
        if not ok:
            print("[ERROR] No se pudo crear el entorno:", out)
            return False

    print(f"[COND] Instalando samtools en '{env_name}' …")
    ok, out = run_cmd([
        "conda", "install", "-y", "-n", env_name,
        "-c", "bioconda", "samtools"
    ])
    if not ok:
        print("[ERROR] conda install samtools:", out)
        return False
    if not samtools_exe.exists():
        print("[ERROR] samtools no apareció en el entorno.")
        return False

    symlink_if_exists("samtools")
    print("[OK] samtools instalado y enlazado.")
    return True


# Entorno python3.8
def _ensure_python38() -> Path | None:
    """Devuelve ruta a python3.8, instalando entorno conda si es necesario."""
    py38 = shutil.which("python3.8")
    if py38:
        return Path(py38)

    env_name = "novowrap38"
    print(f"[INFO] Creando entorno Conda '{env_name}' con python=3.8 …")
    ok, out = run_cmd(["conda", "create", "-y", "-n", env_name, "python=3.8"])
    if not ok:
        print("[ERROR] No se pudo crear el entorno python3.8:", out)
        return None
    conda_prefix = os.environ.get("CONDA_PREFIX", Path.home()/ "miniconda3")
    py38_path = Path(conda_prefix) / "envs" / env_name / "bin" / "python"
    return py38_path if py38_path.exists() else None


# Instalar novowrap
def install_novowrap_pip(force: bool = False) -> bool:
    link = LIB_DIR / "novowrap"
    if link.exists() and not force:
        return True

    py38 = _ensure_python38()
    if not py38:
        return False

    env_bin = py38.parent          
    pip_cmd = [str(py38), "-m", "pip"]

    run_cmd(pip_cmd + ["install", "--upgrade", "pip", "--user"], accept_nonzero=True)
    ok, out = run_cmd(pip_cmd + ["install", "--user", "novowrap"])
    if not ok:
        print("[ERROR] pip no pudo instalar novowrap:", out); return False

    exe = env_bin / "novowrap"
    if not exe.exists():
        exe = Path.home() / ".local" / "bin" / "novowrap"
        if not exe.exists():
            print("[ERROR] Executable novowrap no encontrado tras pip."); return False

    # Enlazar en libreris/
    LIB_DIR.mkdir(parents=True, exist_ok=True)
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(exe.resolve())
    print("[OK] novowrap instalado y enlazado.")
    return True


#  Enlaces simbólicos
def _conda_env_bins() -> List[Path]:
    ok, out = run_cmd(["conda", "env", "list", "--json"])
    if not ok:
        return []
    try:
        data = json.loads(out)
        envs = [Path(p) / "bin" for p in data.get("envs", [])]
        return [e for e in envs if e.exists()]
    except json.JSONDecodeError:
        return []


def symlink_if_exists(binary_name: str) -> bool:

    path: str | None = shutil.which(binary_name)
    if not path and (pre := os.environ.get("CONDA_PREFIX")):
        p = Path(pre) / "bin" / binary_name
        if p.exists():
            path = str(p)

    if not path:
        for b in _conda_env_bins():
            p = b / binary_name
            if p.exists():
                path = str(p)
                break

    if not path:
        print(f"[WARN] {binary_name} no encontrado para enlace.")
        return False

    target = LIB_DIR / binary_name
    real   = Path(path).resolve()

    if real == target.resolve():
        print(f"[OK] {binary_name} ya presente en libreris/ (sin enlace).")
        return True

    if target.exists() or target.is_symlink():
        target.unlink()
    target.symlink_to(real)
    print(f"[OK] Enlace creado para {binary_name}")
    return True


#  Verificación
def verify(name: str) -> Tuple[bool, str]:
    if not PROC_SUPPORTED:
        return False, "Sin subprocesos."

    if name == "Trimmomatic":
        if not TRIMMOMATIC_JAR.exists():
            return False, "JAR no encontrado."
        return run_cmd(["java", "-jar", str(TRIMMOMATIC_JAR), "-version"],
                       accept_nonzero=True)

    if name == "novowrap":
        py38 = _ensure_python38()
        if not py38:
            return False, "No hay python3.8"
        return run_cmd([str(py38), "-m", "novowrap", "-h"], accept_nonzero=True)

    cmd = VERIFY_CMDS.get(name)
    if not cmd:
        return False, "Sin comando verificación."
    return run_cmd(cmd, accept_nonzero=True)


# FLUJO PRINCIPAL
def main() -> None:
    parser = argparse.ArgumentParser(description="Instala/verifica dependencias.")
    parser.add_argument("--force", action="store_true",
                        help="Forzar reinstalación completa")
    args = parser.parse_args()

    if str(LIB_DIR) not in os.environ.get("PATH", ""):
        os.environ["PATH"] = f"{LIB_DIR}{os.pathsep}" + os.environ["PATH"]

    installed, verified, errores = [], [], []

    print("=== Verificando/instalando herramientas ===")
    for name, locator in REPOS.items():
        print(f"\n→ {name} …")
        ok = False

        # ----- build ---------------------------------------------------------
        if locator.startswith(BUILD_PREFIX):
            if name == "Trimmomatic":
                ok = install_trimmomatic_zip(force=args.force)
            elif name == "novowrap":
                ok = install_novowrap_pip(force=args.force)
            elif name == "minimap2":
                ok = install_minimap2_source(force=args.force)
            elif name == "samtools":                        
                ok = install_samtools_env(force=args.force)   

        # ----- conda ---------------------------------------------------------
        elif locator.startswith(CONDA_PREFIX):
            ok = conda_install(locator[len(CONDA_PREFIX):], force=args.force)

        if ok:
            installed.append(name)
        else:
            errores.append(name); continue

        # symlink -------------------------------------------------------------
        if name == "Trimmomatic":
            symlink_if_exists("trimmomatic")
        elif name == "novowrap":
            symlink_if_exists("novowrap")
        elif name == "minimap2":
            symlink_if_exists("minimap2")
        else:
            symlink_if_exists(name)

    # ----------------- verificación -----------------
    print("\n=== Verificación final ===")
    for name in REPOS:
        ok, msg = verify(name)
        print(f"[{'✔' if ok else '✖'}] {name}: {msg.splitlines()[0]}")
        (verified if ok else errores).append(name)

    verified = sorted(set(verified))
    errores  = sorted({e for e in errores if e not in verified})

    print("\n=== Resumen ===")
    if installed:
        print("Instaladas ahora:", ", ".join(installed))
    if verified:
        print("Verificadas:", ", ".join(verified))
    if errores:
        print("Errores/faltantes:", ", ".join(errores))
    if not errores:
        print("Todo correcto. ✅")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrumpido por el usuario.")


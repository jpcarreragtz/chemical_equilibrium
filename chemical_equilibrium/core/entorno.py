"""core.entorno — asegura que el script corra con un Python que tenga numpy,
scipy y matplotlib.

Uso, en cada script, justo después de poner la raíz del proyecto en sys.path
y ANTES de importar numpy, scipy o matplotlib:

    RAIZ = next(p for p in pathlib.Path(__file__).resolve().parents if (p / "core").is_dir())
    sys.path.insert(0, str(RAIZ))
    import core.entorno  # noqa: F401

La raíz se busca hacia arriba por la carpeta que contiene core/, así que el
script funciona a cualquier profundidad (ejercicios/sNN_ejN_tema/, examenes/...).

Al importarse:
  * Si numpy, scipy y matplotlib se pueden importar, no hace nada.
  * Si falta alguno y existe chemical_equilibrium/.venv/bin/python, vuelve a
    lanzar el mismo script con ese intérprete (os.execv, mismos argumentos).
    La variable de entorno CHEMEQ_ENTORNO_RELANZADO evita un ciclo.
  * Si no hay .venv, o el .venv tampoco los tiene, o no se está corriendo un
    archivo (consola interactiva, kernel de Jupyter, `python -c`), termina con
    un mensaje que dice qué falta y el comando para correrlo.

Solo usa la biblioteca estándar y funciona desde Python 3.6.
"""

import importlib
import os
import sys
from pathlib import Path

REQUERIDOS = ("numpy", "scipy", "matplotlib")
GUARDA = "CHEMEQ_ENTORNO_RELANZADO"
RAIZ = Path(__file__).resolve().parents[1]            # chemical_equilibrium/
VENV_PYTHON = RAIZ / ".venv" / "bin" / "python"


def faltantes():
    """Paquetes de REQUERIDOS que este intérprete no puede importar."""
    faltan = []
    for nombre in REQUERIDOS:
        try:
            importlib.import_module(nombre)
        except Exception:                 # ausente o instalación rota
            faltan.append(nombre)
    return faltan


def _es_archivo_en_ejecucion():
    """True si se corre `python archivo.py` (no consola, kernel ni -c/-m)."""
    if hasattr(sys, "ps1") or "ipykernel" in sys.modules:
        return False
    principal = getattr(sys.modules.get("__main__"), "__file__", None)
    return bool(principal and sys.argv and sys.argv[0]
                and os.path.isfile(sys.argv[0])
                and os.path.samefile(sys.argv[0], principal))


def _mensaje(faltan, motivo):
    script = sys.argv[0] if sys.argv and sys.argv[0] else "archivo.py"
    lineas = [
        "",
        "ERROR: este Python no tiene " + ", ".join(faltan) + ".",
        "  Intérprete: " + sys.executable + " (Python "
        + ".".join(map(str, sys.version_info[:3])) + ")",
        "  " + motivo,
    ]
    if VENV_PYTHON.exists():
        lineas += [
            "  Corre el script con el entorno del proyecto:",
            '    "' + str(VENV_PYTHON) + '" "' + script + '"',
            "  Si al entorno le faltan paquetes:",
            '    "' + str(VENV_PYTHON) + '" -m pip install numpy scipy matplotlib',
        ]
    else:
        lineas += [
            "  Crea el entorno del proyecto (una sola vez) y corre el script:",
            '    python3 -m venv "' + str(RAIZ / ".venv") + '"',
            '    "' + str(VENV_PYTHON) + '" -m pip install numpy scipy matplotlib',
            '    "' + str(VENV_PYTHON) + '" "' + script + '"',
        ]
    lineas.append("  En VS Code: Python: Select Interpreter -> .venv de "
                  "chemical_equilibrium.")
    return "\n".join(lineas)


def asegurar():
    """Relanza con el .venv si hace falta; si no se puede, termina con mensaje."""
    faltan = faltantes()
    if not faltan:
        return
    if not VENV_PYTHON.exists():
        motivo = "No existe " + str(VENV_PYTHON) + "."
    elif os.environ.get(GUARDA):
        motivo = "Ya se relanzó con el .venv del proyecto y siguen faltando."
    elif not _es_archivo_en_ejecucion():
        motivo = ("No se puede relanzar solo desde una consola interactiva o "
                  "un kernel de Jupyter: cambia de intérprete o de kernel.")
    else:
        os.environ[GUARDA] = "1"
        sys.stdout.flush()
        sys.stderr.flush()
        venv = str(VENV_PYTHON)           # sin resolver el enlace: conserva el venv
        os.execv(venv, [venv] + sys.argv)
    raise SystemExit(_mensaje(faltan, motivo))


asegurar()

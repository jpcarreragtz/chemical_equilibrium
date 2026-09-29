"""core — motor del solver de equilibrio químico (no se toca salvo bugs).

Expone las rutas del proyecto para que módulos y scripts no dependan del
directorio de trabajo:

    from core import ROOT, DATA, SALIDAS
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]   # chemical_equilibrium/
DATA = ROOT / "data"                         # sva_tables.json, raw/
SALIDAS = ROOT / "salidas"                   # CSV/PNG regenerables

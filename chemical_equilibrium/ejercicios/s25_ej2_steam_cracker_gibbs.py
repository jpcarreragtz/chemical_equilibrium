"""Sesión 25, Ejercicio 2 — Minimización de la energía de Gibbs.

Steam cracker a 1000 K y 1 bar, con 1 mol C2H6 + 4 mol H2O y 9 especies.
ΔG°f de JANAF; el C2H6 no está en JANAF y sale de la ec. 13.18 de SVA.
La misma solución en GAMS está en s25_ej2_gibbs.gms.
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto

from core.gibbs import cargar_dgf, imprimir_resultado, minimizar_gibbs, validar_gibbs

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
ESPECIES = ["CH4", "C2H4", "C2H2", "CO2", "CO", "O2", "H2", "H2O", "C2H6"]
ALIMENTACION = {"C2H6": 1.0, "H2O": 4.0}  # mol
T = 1000.0                                # K
P = 1.0                                   # bar

# %% ΔG°f y su fuente
dGf, fuentes = cargar_dgf(ESPECIES, T)
print(f"ΔG°f a {T:g} K (J/mol)")
for sp in ESPECIES:
    print(f"  {sp:<6}{dGf[sp]:>12.0f}   {fuentes[sp]}")

# %% Equilibrio
res = minimizar_gibbs(ESPECIES, dGf, ALIMENTACION, T, P)
imprimir_resultado(res)
validar_gibbs(res)

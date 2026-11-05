"""Sesión 25, Ejercicio 1 — Minimización de la energía de Gibbs.

Reformado de metano a 1000 K y 1 bar, con 2 mol CH4 + 3 mol H2O. Es el
Ejemplo 13.14 de Smith–Van Ness. La diapositiva remite a JANAF para los ΔG°f,
así que aquí se resuelve con JANAF y se contrasta con los valores del libro.
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto

from core.gibbs import cargar_dgf, imprimir_resultado, minimizar_gibbs, validar_gibbs

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
ESPECIES = ["CH4", "H2O", "CO", "CO2", "H2"]
ALIMENTACION = {"CH4": 2.0, "H2O": 3.0}   # mol
T = 1000.0                                # K
P = 1.0                                   # bar
DGF_LIBRO = {"CH4": 19720.0, "H2O": -192420.0, "CO": -200240.0,
             "CO2": -395790.0, "H2": 0.0}  # J/mol, SVA Ej. 13.14

# %% ΔG°f y su fuente
dGf, fuentes = cargar_dgf(ESPECIES, T)
print(f"ΔG°f a {T:g} K (J/mol)")
for sp in ESPECIES:
    print(f"  {sp:<6}{dGf[sp]:>12.0f}   {fuentes[sp]}")

# %% Equilibrio con JANAF
print("\nCon JANAF")
res = minimizar_gibbs(ESPECIES, dGf, ALIMENTACION, T, P)
imprimir_resultado(res)
validar_gibbs(res)

# %% Contraste con los valores del libro
print("\nCon los ΔG°f del Ejemplo 13.14 de SVA")
libro = minimizar_gibbs(ESPECIES, DGF_LIBRO, ALIMENTACION, T, P)
imprimir_resultado(libro)

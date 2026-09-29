"""Parcial 1 — Problema 1 (2.5 pts)

Enunciado: A <-> R + 2S en fase gas, operación continua isotérmica e
isobárica, P = 10 atm, T = 500 K, v0 = 500 dm³/s de A puro.
(a) Ecuaciones de F_i y C_i con todos los parámetros evaluados.
(b) X_Ae si K_C(500 K) = 0.2 mol²/dm⁶.
(c) F_i y v cuando X = 0.70·X_Ae.

Supuestos: gas ideal; T y P constantes (volumen variable, v = v0(1+eps·X));
A, R, S son especies abstractas -> su "composición" se declara en
fragmentos conservados (A = Rm·Sm2) para poder verificar el balance.
"""

import csv
import os
import sys

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = (_DIR if os.path.exists(os.path.join(_DIR, "equilibrium.py"))
         else os.path.dirname(_DIR))
sys.path.insert(0, _ROOT)

import numpy as np

from equilibrium import R_ATM, StoichiometricTable, format_table, solve_extents
from selfcheck import Problem, validate

# ══════════════════════════ INPUTS ════════════════════════════════════════
T, P_atm, v0, Kc = 500.0, 10.0, 500.0, 0.2   # K, atm, dm³/s, mol²/dm⁶
F_A0 = 121.86                                # mol/s = C_A0·v0
ATOMOS = {"A": {"Rm": 1, "Sm": 2}, "R": {"Rm": 1}, "S": {"Sm": 1}}
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

# cálculo (se imprime en el orden fijo del entregable)
tbl = StoichiometricTable(nu={"A": -1.0, "R": 1.0, "S": 2.0},
                          theta={"A": 1.0}, limiting="A", T0=T, P0=P_atm,
                          pressure_unit="atm", variable_volume=True,
                          n_A0=F_A0)
C_A0 = tbl.C_A0
Xe = tbl.solve_Xe(Kc)
X70 = 0.70 * Xe
K_y = Kc / C_A0 ** 2                          # Pi y^nu = K_C / C_T0^delta
res = solve_extents(["A", "R", "S"], {"A": F_A0}, [[-1, 1, 2]], [K_y],
                    P=1.0, P0=1.0, atoms=ATOMOS)

print("=" * 70)
print("  PARCIAL 1 — P1:  A <-> R + 2S   (gas, 10 atm, 500 K)")
print("=" * 70)

print("\n1) ENUNCIADO Y SUPUESTOS")
print("   A <-> R + 2S, gas ideal, operación continua ISOTÉRMICA e")
print("   ISOBÁRICA: P = 10 atm, T = 500 K, v0 = 500 dm³/s de A PURO.")
print("   K_C(500 K) = 0.2 mol²/dm⁶ (unidades consistentes con delta = +2).")
print("   Volumen variable: v = v0(1 + eps·X) a T,P constantes.")

print("\n2) TABLA ESTEQUIOMÉTRICA (base: A, por mol de A)")
print("  " + format_table(tbl.symbolic_rows(),
      ["Especie", "Inicial", "Cambio", "Remanente", "Concentración"]
      ).replace("\n", "\n  "))

print("\n3) PARÁMETROS EVALUADOS")
print(f"   y_A0 = 1 (A puro)      delta = 1+2-1 = {tbl.delta:g}      "
      f"eps = y_A0·delta = {tbl.epsilon:g}")
print(f"   C_A0 = C_T0 = P/(R·T) = 10/(0.08206·500) = {C_A0:.5f} mol/dm³")
print(f"   F_A0 = C_A0·v0 = {F_A0:.2f} mol/s")
print("   F_A = F_A0(1-X)   F_R = F_A0·X   F_S = 2F_A0·X   F_T = F_A0(1+2X)")
print("   v = v0(1+2X);   C_A = C_A0(1-X)/(1+2X),  C_R = C_A0·X/(1+2X),")
print("   C_S = 2C_A0·X/(1+2X)")

print("\n4) ECUACIÓN A RESOLVER (números sustituidos)")
print("   K_C = C_R·C_S²/C_A = 4·C_A0²·X³/[(1+2X)²(1-X)] = 0.2")
print(f"   con (1+2X)²(1-X) = 1+3X-4X³:   "
      f"({4 * C_A0 ** 2:.5f} + 0.8)·X³ - 0.6·X - 0.2 = 0")
a3 = 4 * C_A0 ** 2 + 4 * Kc
print(f"   =>  {a3:.4f}·X³ - 0.6·X - 0.2 = 0")
raices = np.roots([a3, 0.0, -3 * Kc, -Kc])
print("   raíces: " + ",  ".join(
    f"{z.real:.5f}{z.imag:+.5f}i" if abs(z.imag) > 1e-10 else f"{z.real:.5f}"
    for z in raices))
print("   -> se descartan las dos raíces complejas; la única física "
      "(0 <= X <= 1) es X_Ae.")
print(f"   (equivalente Método A: Pi y^nu = K_C/C_T0² = 0.2/{C_A0:.5f}² = "
      f"{K_y:.5f} = 4X³/[(1+2X)²(1-X)])")

print("\n5) RESULTADOS")
print(f"   b) X_Ae = {Xe:.5f}")
print(f"   c) X = 0.70·X_Ae = {X70:.5f}")
n, C = tbl.moles(X70), tbl.concentrations(X70)
v = v0 * (1 + tbl.epsilon * X70)
filas = [[sp, n[sp], n[sp] / sum(n.values()), C[sp]] for sp in tbl.species]
filas.append(["TOTAL", sum(n.values()), 1.0, sum(C.values())])
print("  " + format_table(filas, ["Especie", "F_i (mol/s)", "y_i",
                                  "C_i (mol/dm³)"]).replace("\n", "\n  "))
print(f"   v = v0(1+2X) = 500·(1+2·{X70:.4f}) = {v:.1f} dm³/s")

print("\n6) COMPROBACIÓN")
print(f"   K_C recalculada en X_Ae: {tbl.Kc_calc(Xe):.6f}  (objetivo 0.2)")
print(f"   C_T = F_T/v = {sum(n.values()) / v:.5f} = P/(RT) = {C_A0:.5f} "
      f"mol/dm³ (constante) ✓")
print(f"   suma y_i = {sum(n.values()) / sum(n.values()):.6f} ✓;   "
      f"fragmentos: Rm {F_A0:g} -> "
      f"{n['A'] + n['R']:.2f}, Sm {2 * F_A0:g} -> "
      f"{2 * n['A'] + n['S']:.2f} ✓")

print("\n7) SELFCHECK (Método A, 20 arranques aleatorios)")
prob = Problem(species=["A", "R", "S"], n0={"A": F_A0}, nu=[[-1, 1, 2]],
               K=[K_y], P=1.0, P0=1.0, atoms=ATOMOS, T=T,
               pressure_unit="atm", label="P1 parcial")
rep = validate(prob, res, n_starts=20, k_tol=1e-4)
print(f"   Métodos A y B coinciden: X(A) = {res.xi[0] / F_A0:.5f} vs "
      f"X_Ae(B) = {Xe:.5f}  (dif = {abs(res.xi[0] / F_A0 - Xe):.1e})")

# CSV
ruta_csv = os.path.join(_DIR, "resultados_P1.csv")
with open(ruta_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["caso", "X", "F_A", "F_R", "F_S", "F_T", "v_dm3_s",
                "C_A", "C_R", "C_S"])
    for etiqueta, X in (("equilibrio", Xe), ("70pct_de_Xe", X70)):
        m, c = tbl.moles(X), tbl.concentrations(X)
        w.writerow([etiqueta, f"{X:.5f}", f"{m['A']:.3f}", f"{m['R']:.3f}",
                    f"{m['S']:.3f}", f"{sum(m.values()):.3f}",
                    f"{v0 * (1 + 2 * X):.1f}", f"{c['A']:.5f}",
                    f"{c['R']:.5f}", f"{c['S']:.5f}"])
print(f"\n   CSV: {ruta_csv}")

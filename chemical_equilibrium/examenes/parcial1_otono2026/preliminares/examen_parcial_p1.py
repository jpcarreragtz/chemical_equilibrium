"""Parcial 1 — Problema 1 (2.5 pts): A <-> R + 2S, gas, continuo,
isotérmico e isobárico. P = 10 atm, T = 500 K, v0 = 500 dm³/s de A puro.

a) Ecuaciones de flujo molar y concentración de todas las especies +
   parámetros (C_A0, F_A0, delta, eps, ...).
b) X_Ae con K_C = 0.2 mol²/dm⁶.
c) F_i y v cuando X = 70 % de X_Ae.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))  # raíz del proyecto

from core.equilibrium import (R_ATM, StoichiometricTable, format_table,
                         solve_extents)
from core.selfcheck import Problem, validate

T, P_atm, v0, Kc = 500.0, 10.0, 500.0, 0.2      # K, atm, dm3/s, mol2/dm6

# Tabla de Fogler por mol de A (A puro, y_A0 = 1); n_A0 = F_A0 en mol/s
C_A0 = P_atm / (R_ATM * T)                       # = C_T0 (A puro)
F_A0 = C_A0 * v0
tbl = StoichiometricTable(nu={"A": -1.0, "R": 1.0, "S": 2.0},
                          theta={"A": 1.0}, limiting="A",
                          T0=T, P0=P_atm, pressure_unit="atm",
                          variable_volume=True, n_A0=F_A0)

print("=" * 70)
print("  P1: A <-> R + 2S   (10 atm, 500 K, v0 = 500 dm3/s, A puro)")
print("=" * 70)
print("\n  a) Parámetros y ecuaciones (base: A, y_A0 = 1):")
print(f"     C_A0 = C_T0 = P/(R·T) = 10/(0.08206·500) = {C_A0:.5f} mol/dm3")
print(f"     F_A0 = C_A0·v0 = {F_A0:.3f} mol/s")
print(f"     delta = 1 + 2 - 1 = {tbl.delta:g}     "
      f"eps = y_A0·delta = {tbl.epsilon:g}")
print("     Flujos:   F_A = F_A0(1-X)   F_R = F_A0·X   F_S = 2·F_A0·X")
print("               F_T = F_A0(1+2X)          v = v0(1+2X)")
print("     Concentraciones (T, P constantes):  C_i = F_i/v:")
print("       C_A = C_A0(1-X)/(1+2X)   C_R = C_A0·X/(1+2X)   "
      "C_S = 2·C_A0·X/(1+2X)")
print("\n  Tabla estequiométrica simbólica:")
print("  " + format_table(tbl.symbolic_rows(),
      ["Especie", "Inicial", "Cambio", "Remanente", "Concentración"]
      ).replace("\n", "\n  "))

# ---------------- b) conversión de equilibrio ----------------
Xe = tbl.solve_Xe(Kc)
print("\n  b) Equilibrio:  K_C = C_R·C_S²/C_A = "
      "4·C_A0²·X³/[(1+2X)²(1-X)] = 0.2 mol²/dm⁶")
print(f"     4·C_A0² = {4 * C_A0 ** 2:.5f}  ->  "
      f"{4 * C_A0 ** 2:.5f}·X³ = 0.2·(1+2X)²·(1-X)")
print(f"     X_Ae = {Xe:.5f}")
print(f"     verificación: K_C(X_Ae) = {tbl.Kc_calc(Xe):.6f} "
      f"(objetivo 0.2)")

# Contraste independiente por Método A (avance de reacción) + selfcheck.
# K_y = Pi y^nu = K_C / C_T^delta  (C_T = P/RT = C_A0 porque A es puro).
# Átomos: A se trata como el "compuesto" R·S2 para que el balance cierre.
K_y = Kc / C_A0 ** 2
atomos_moieties = {"A": {"Rm": 1, "Sm": 2}, "R": {"Rm": 1}, "S": {"Sm": 1}}
res = solve_extents(["A", "R", "S"], {"A": 1.0}, [[-1, 1, 2]], [K_y],
                    P=1.0, P0=1.0, atoms=atomos_moieties)
print(f"\n  Contraste Método A:  Pi y^nu = K_C/C_T0^2 = {K_y:.5f}  "
      f"->  X = xi = {res.xi[0]:.5f}")
prob = Problem(species=["A", "R", "S"], n0={"A": 1.0}, nu=[[-1, 1, 2]],
               K=[K_y], P=1.0, P0=1.0, atoms=atomos_moieties, T=T,
               pressure_unit="atm", label="P1 parcial")
validate(prob, res, n_starts=20, k_tol=1e-4)

# ---------------- c) al 70 % de la conversión de equilibrio ----------------
X70 = 0.70 * Xe
n, C = tbl.moles(X70), tbl.concentrations(X70)
v = v0 * (1 + tbl.epsilon * X70)
print(f"\n  c) X = 0.70·X_Ae = {X70:.5f}")
filas = [[sp, n[sp], C[sp]] for sp in tbl.species]
filas.append(["TOTAL", sum(n.values()), sum(C.values())])
print("  " + format_table(filas, ["Especie", "F_i (mol/s)", "C_i (mol/dm3)"]
                          ).replace("\n", "\n  "))
print(f"     v = v0(1 + eps·X) = 500·(1 + 2·{X70:.5f}) = {v:.2f} dm3/s")
print(f"     check: C_T = F_T/v = {sum(n.values()) / v:.5f} = P/(RT) = "
      f"{C_A0:.5f} mol/dm3")

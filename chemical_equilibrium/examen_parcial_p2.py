"""Parcial 1 — Problema 2 (3.0 pts): oxidación de etileno a acetaldehído.

Balanceadas (en el examen aparecen esqueléticas):
    R1: C2H4 + 1/2 O2 = C2H4O        (acetaldehído, CH3CHO)
    R2: C2H4 + 3 O2   = 2 CO2 + 2 H2O

Alimentación: 1.5 mol C2H4 + 3 mol O2, 1 bar. A la salida se MIDEN
0.78 mol de acetaldehído y 0.35 mol de H2O. NO es problema de
equilibrio: los avances salen de las mediciones.
"""

from equilibrium import format_table
from reactions import build_system
from selfcheck import Problem, Report

reacciones = ["C2H4 + 1/2 O2 = C2H4O", "C2H4 + 3 O2 = 2 CO2 + 2 H2O"]
alim = {"C2H4": 1.5, "O2": 3.0}
especies, nu, atomos, _ = build_system(reacciones, alim)   # verifica balance

print("=" * 70)
print("  P2: acetaldehído por oxidación de etileno (1.5 C2H4 + 3 O2)")
print("=" * 70)
for j, rx in enumerate(reacciones):
    print(f"  R{j + 1}: {rx}   (balanceada; en el examen viene esquelética)")

# ---------------- a) avances desde las mediciones ----------------
# El acetaldehído SOLO se forma en R1 (coef +1):  n_C2H4O = xi1 = 0.78
# El agua SOLO se forma en R2 (coef +2):          n_H2O = 2·xi2 = 0.35
xi1 = 0.78 / 1.0
xi2 = 0.35 / 2.0
print(f"\n  a) xi_1 = n_C2H4O/1 = 0.78 mol      "
      f"xi_2 = n_H2O/2 = 0.35/2 = {xi2:.3f} mol")

n = {sp: alim.get(sp, 0.0) + nu[0][i] * xi1 + nu[1][i] * xi2
     for i, sp in enumerate(especies)}
n_T = sum(n.values())
delta1 = sum(nu[0])
print(f"     n_T = 4.5 + ({delta1:g})·xi_1 + (0)·xi_2 = {n_T:.4f} mol")
expr = {"C2H4": "1.5 - xi1 - xi2", "O2": "3 - 0.5·xi1 - 3·xi2",
        "C2H4O": "xi1", "CO2": "2·xi2", "H2O": "2·xi2"}
print("  " + format_table(
    [[sp, expr[sp], n[sp], n[sp] / n_T] for sp in especies],
    ["Especie", "balance", "n_i (mol)", "y_i"]).replace("\n", "\n  "))

# checks de captura (átomos entrada vs salida)
from equilibrium import atom_totals
ent, sal = atom_totals(alim, atomos), atom_totals(n, atomos)
print("     check átomos in -> out: " + ", ".join(
    f"{e}: {ent.get(e, 0):g} -> {sal.get(e, 0):g}"
    for e in sorted(set(ent) | set(sal))))

# ---------------- b) a e) ----------------
X = (alim["C2H4"] - n["C2H4"]) / alim["C2H4"]
y_O2 = n["O2"] / n_T
Y = n["C2H4O"] / (alim["C2H4"] * 1.0)     # rendimiento: n_des/(n_A0·coef)
productos = n["C2H4O"] + n["CO2"] + n["H2O"]
S = n["C2H4O"] / productos                # respecto a TODOS los productos
print(f"\n  b) X_C2H4 = (1.5 - {n['C2H4']:.3f})/1.5 = {X:.5f}")
print(f"  c) y_O2 = {n['O2']:.3f}/{n_T:.2f} = {y_O2:.5f}")
print(f"  d) Rendimiento a acetaldehído = n_C2H4O/(n_C2H4,0·1) "
      f"= 0.78/1.5 = {Y:.5f}")
print(f"  e) Selectividad respecto a TODOS los productos = "
      f"n_C2H4O/(n_C2H4O + n_CO2 + n_H2O) = 0.78/{productos:.2f} = {S:.5f}")
print(f"     (alternativa deseado/no-deseado: 0.78/(0.35+0.35) = "
      f"{0.78 / 0.70:.4f} — usa la definición vista en clase)")

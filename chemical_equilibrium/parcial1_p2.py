"""Parcial 1 — Problema 2 (3.0 pts)

Enunciado: oxidación de etileno; rutas paralelas al acetaldehído (deseado)
y a combustión completa. Balanceadas (en el examen vienen esqueléticas):
    R1: C2H4 + 1/2 O2 = C2H4O  (CH3CHO)          delta_1 = -1/2
    R2: C2H4 + 3 O2   = 2 CO2 + 2 H2O            delta_2 = 0
Alimentación 1.5 mol C2H4 + 3 mol O2 a 1 bar. A la salida se MIDEN
0.78 mol de acetaldehído y 0.35 mol de H2O.
(a) avances y moles de salida; (b) X_C2H4; (c) y_O2; (d) rendimiento al
deseado; (e) selectividad del acetaldehído.

Supuestos: NO es problema de equilibrio — los avances se despejan de las
mediciones (sistema lineal). Fase gas (solo afecta y_i, no los balances).
"""

import csv
import os
import sys

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = (_DIR if os.path.exists(os.path.join(_DIR, "equilibrium.py"))
         else os.path.dirname(_DIR))
sys.path.insert(0, _ROOT)

from equilibrium import atom_totals, format_table
from reactions import build_system
from selfcheck import Report

# ══════════════════════════ INPUTS ════════════════════════════════════════
reacciones = ["C2H4 + 1/2 O2 = C2H4O", "C2H4 + 3 O2 = 2 CO2 + 2 H2O"]
alim = {"C2H4": 1.5, "O2": 3.0}
medidos = {"C2H4O": 0.78, "H2O": 0.35}      # mol a la salida
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

especies, nu, atomos, _ = build_system(reacciones, alim)  # verifica balance
xi1 = medidos["C2H4O"] / 1.0                 # solo R1 produce acetaldehído
xi2 = medidos["H2O"] / 2.0                   # solo R2 produce agua
n = {sp: alim.get(sp, 0.0) + nu[0][i] * xi1 + nu[1][i] * xi2
     for i, sp in enumerate(especies)}
n_T = sum(n.values())

print("=" * 70)
print("  PARCIAL 1 — P2: acetaldehído por oxidación de etileno")
print("=" * 70)

print("\n1) ENUNCIADO Y SUPUESTOS")
print("   Fase gas, 1 bar. Alimentación: 1.5 mol C2H4 + 3 mol O2.")
print("   Salida MEDIDA: 0.78 mol CH3CHO y 0.35 mol H2O -> NO es")
print("   equilibrio: los avances se despejan de las mediciones.")
for j, rx in enumerate(reacciones):
    print(f"   R{j + 1}: {rx}   (balanceada; en el examen viene esquelética)")

print("\n2) TABLA DE AVANCES (n_i = n_i0 + sum nu_ij·xi_j)")
expr = {"C2H4": "1.5 - xi1 - xi2", "O2": "3 - 0.5·xi1 - 3·xi2",
        "C2H4O": "0 + xi1", "CO2": "0 + 2·xi2", "H2O": "0 + 2·xi2"}
print("  " + format_table([[sp, alim.get(sp, 0.0), expr[sp]]
                           for sp in especies],
      ["Especie", "n_i0 (mol)", "n_i(xi1, xi2)"]).replace("\n", "\n  "))
print("   n_T = 4.5 - 0.5·xi1   (delta_1 = -1/2, delta_2 = 0)")

print("\n3) PARÁMETROS")
print("   nu (filas R1, R2; columnas " + ", ".join(especies) + "):")
for j, fila in enumerate(nu):
    print(f"     R{j + 1}: {fila}")
print(f"   n_i0 = {alim} (resto 0)")

print("\n4) SISTEMA A RESOLVER (números sustituidos, lineal)")
print("   n_C2H4O = xi1        = 0.78  ->  xi1 = 0.78 mol")
print("   n_H2O   = 2·xi2      = 0.35  ->  xi2 = 0.35/2 = 0.175 mol")

print("\n5) RESULTADOS")
print(f"   xi1 = {xi1:.3f} mol,  xi2 = {xi2:.3f} mol,  "
      f"n_T = 4.5 - 0.5·{xi1:.2f} = {n_T:.4f} mol")
print("  " + format_table([[sp, n[sp], n[sp] / n_T] for sp in especies],
      ["Especie", "n_i (mol)", "y_i"]).replace("\n", "\n  "))
X = (alim["C2H4"] - n["C2H4"]) / alim["C2H4"]
y_O2 = n["O2"] / n_T
rend = n["C2H4O"] / (alim["C2H4"] * 1.0)
S_glob = n["C2H4O"] / (n["CO2"] + n["H2O"])
S_CO2 = n["C2H4O"] / n["CO2"]
S_H2O = n["C2H4O"] / n["H2O"]
print(f"   b) X_C2H4 = (1.5 - {n['C2H4']:.3f})/1.5 = {X:.4f}")
print(f"   c) y_O2 = {n['O2']:.3f}/{n_T:.2f} = {y_O2:.4f}")
print(f"   d) Rendimiento al deseado = n_C2H4O/(n_C2H4,0·1) = 0.78/1.5 "
      f"= {rend:.4f}")
print(f"      (definición del curso: mol deseado / mol que se formaría si "
      f"TODO el limitante fuera a R1)")
print(f"   e) Selectividad global = deseado/no-deseado = "
      f"0.78/({n['CO2']:.2f}+{n['H2O']:.2f}) = {S_glob:.4f}")
print(f"      vs CO2: 0.78/0.35 = {S_CO2:.4f}      "
      f"vs H2O: 0.78/0.35 = {S_H2O:.4f}")
print("      (otra convención vista en clase: respecto a TODOS los "
      f"productos = 0.78/1.48 = {n['C2H4O'] / 1.48:.4f} — indica cuál usas)")

print("\n6) COMPROBACIÓN")
ent, sal = atom_totals(alim, atomos), atom_totals(n, atomos)
print("   átomos in -> out:  " + ",  ".join(
    f"{e}: {ent.get(e, 0):.2f} -> {sal.get(e, 0):.2f}"
    for e in sorted(set(ent) | set(sal))))
print(f"   suma y_i = {sum(n.values()) / n_T:.6f};   mín n_i = "
      f"{min(n.values()):.3f} mol >= 0")

print("\n7) SELFCHECK (adaptado: sin K, no es equilibrio)")
rep = Report()
rep.add("pass", 1, "reacciones balanceadas por el parser (C, H, O)")
peor = max(abs(ent.get(e, 0.0) - sal.get(e, 0.0))
           for e in set(ent) | set(sal))
rep.add("pass" if peor < 1e-9 else "fail", 7,
        f"átomos in vs out: máx |Δ| = {peor:.1e} mol")
rep.add("pass" if min(n.values()) >= 0 else "fail", 7,
        f"n_i >= 0 (mín = {min(n.values()):.3f} mol)")
rep.add("na", 6, "residuos de K: no aplica (avances medidos, no equilibrio)")
rep.verdict()

# CSV
ruta_csv = os.path.join(_DIR, "resultados_P2.csv")
with open(ruta_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["xi1", "xi2"] + [f"n_{sp}" for sp in especies]
               + [f"y_{sp}" for sp in especies]
               + ["n_T", "X_C2H4", "y_O2", "rendimiento", "S_global",
                  "S_vs_CO2", "S_vs_H2O"])
    w.writerow([xi1, xi2] + [f"{n[sp]:.4f}" for sp in especies]
               + [f"{n[sp] / n_T:.5f}" for sp in especies]
               + [f"{n_T:.4f}", f"{X:.5f}", f"{y_O2:.5f}", f"{rend:.4f}",
                  f"{S_glob:.4f}", f"{S_CO2:.4f}", f"{S_H2O:.4f}"])
print(f"\n   CSV: {ruta_csv}")

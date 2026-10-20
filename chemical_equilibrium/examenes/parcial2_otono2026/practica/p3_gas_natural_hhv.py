"""Ejercicios adicionales 2do parcial — PROBLEMA 3 (calor de combustión de
gases naturales CH4/C2H6/C3H8/N2, agua LÍQUIDA como producto → HHV)

Planteamiento (a mano)
----------------------
Reacciones de combustión (1 mol de cada combustible, H2O(l)):
    CH4  + 2   O2 →   CO2 + 2 H2O(l)
    C2H6 + 7/2 O2 → 2 CO2 + 3 H2O(l)
    C3H8 + 5   O2 → 3 CO2 + 4 H2O(l)
    N2: inerte (no aporta)
Ley de Hess con ΔH°f (Tabla C.4 SVA; H2O(l) = −285 830 J/mol):
    ΔH°c,i = Σ ν ΔH°f(productos) − ΔH°f(combustible)
Expresión pedida (base 1 mol de gas natural):
    ΔH°c(y) = y_CH4·ΔH°c,CH4 + y_C2H6·ΔH°c,C2H6 + y_C3H8·ΔH°c,C3H8
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz

from core.energia import calor_combustion, combustion_mezcla, HF_H2O_LIQ
from core.equilibrium import format_table

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
GASES = {
    "a": {"CH4": 0.95, "C2H6": 0.02, "C3H8": 0.02, "N2": 0.01},
    "b": {"CH4": 0.90, "C2H6": 0.05, "C3H8": 0.03, "N2": 0.02},
    "c": {"CH4": 0.85, "C2H6": 0.07, "C3H8": 0.03, "N2": 0.05},
}
CLAVE = {"a": -921.71, "b": -946.16, "c": -923.87}     # kJ/mol (maestra)
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

print("=" * 72)
print("  PROBLEMA 3 — ΔH°c de gases naturales (H2O líquida → HHV)")
print("=" * 72)

print("\n1) ΔH°c de cada combustible (ley de Hess, ΔH°f de la Tabla C.4)")
filas = []
c = {}
for sp in ("CH4", "C2H6", "C3H8"):
    r = calor_combustion(sp)
    c[sp] = r
    rx = " + ".join(f"{-v:g} {k}" for k, v in r["nu"].items() if v < 0)
    px = " + ".join(f"{v:g} {k}" for k, v in r["nu"].items() if v > 0)
    filas.append([f"{rx} → {px}", r["dHc_l"] / 1000, r["dHc_g"] / 1000])
print("  " + format_table(filas, ["reacción", "ΔH°c H2O(l) kJ/mol",
                                  "ΔH°c H2O(g) kJ/mol"]).replace("\n", "\n  "))
print(f"   (ΔH°f H2O(l) = {HF_H2O_LIQ:.0f} J/mol)")

e = {sp: c[sp]["dHc_l"] / 1000 for sp in c}
print("\n2) EXPRESIÓN GENERAL (kJ/mol de gas natural):")
print(f"   ΔH°c = {e['CH4']:.3f}·y_CH4 + ({e['C2H6']:.3f})·y_C2H6 + "
      f"({e['C3H8']:.3f})·y_C3H8")

print("\n3) RESULTADOS")
filas, res = [], {}
for k, y in GASES.items():
    m = combustion_mezcla(y)
    res[k] = m["dHc_l"] / 1000
    filas.append([k, y["CH4"], y["C2H6"], y["C3H8"], y["N2"], res[k],
                  CLAVE[k], res[k] - CLAVE[k]])
print("  " + format_table(filas, ["inciso", "y_CH4", "y_C2H6", "y_C3H8",
                                  "y_N2", "ΔH°c kJ/mol", "maestra",
                                  "diferencia"]).replace("\n", "\n  "))
mayor = min(res, key=res.get)
print(f"\n   → Mayor calor de combustión (|ΔH°c| máximo): inciso {mayor}) "
      f"con {res[mayor]:.2f} kJ/mol")

# ¿qué composición reproduce el −923.87 de la clave en c)?
yc = GASES["c"]
yCH4_imp = (-CLAVE["c"] + (e["C2H6"] * yc["C2H6"] + e["C3H8"] * yc["C3H8"])
            ) / -e["CH4"]
print("\n   NOTA inciso c) (DISCREPANCIA CON LA CLAVE)")
print(f"   0.85(890.649) + 0.07(1560.688) + 0.03(2219.167) = "
      f"{-res['c']:.2f} kJ/mol")
print(f"   La clave −923.87 es −932.87 con los dígitos '3' y '2' "
      f"intercambiados (error de")
print(f"   transcripción; para obtenerla de verdad haría falta y_CH4 ≈ "
      f"{yCH4_imp:.4f}).")
print("   La conclusión no cambia: el gas b) tiene el mayor calor de combustión.")

print("\nAUTOVALIDACIÓN")
ok = True
for k in ("a", "b"):
    b = abs(res[k] - CLAVE[k]) < 0.05
    ok &= b
    print(f"  [{'PASS' if b else 'FAIL'}] {k}) {res[k]:.2f} vs clave "
          f"{CLAVE[k]} (tol 0.05 kJ/mol)")
for sp, r in c.items():
    b = abs(r["check_HHV_LHV"]) < 1e-6
    ok &= b
    print(f"  [{'PASS' if b else 'FAIL'}] {sp}: HHV − LHV = n_H2O·ΔH_lv "
          f"({r['n_H2O']:g}·44 012 J)")
for k, y in GASES.items():
    b = abs(sum(y.values()) - 1) < 1e-12
    ok &= b
    print(f"  [{'PASS' if b else 'FAIL'}] {k}) Σ y_i = {sum(y.values()):.4f}")
print(f"  [info] c) {res['c']:.2f} vs clave {CLAVE['c']}: ver NOTA")
print(f"\nVEREDICTO: {'OK' if ok else 'REVISAR'}")

"""Ejercicios adicionales 2do parcial — PROBLEMA 4 (reformado de metano,
dos reacciones, proceso NO adiabático)

    R1: CH4 + H2O -> CO + 3 H2        R2: CO + H2O -> CO2 + H2
Alimentación CH4 + H2O a 1 atm y 500 °C; productos a 850 °C con
    y_CO2 = 0.0275, y_CO = 0.1725, y_H2O = 0.1725, y_H2 = 0.6275
¿Calor agregado por mol de gas producido?

Planteamiento (a mano)  — base: 1 mol de gas producido
------------------------------------------------------
1) Balance de materia (átomos):
     C:  n_CH4,0 = y_CO2 + y_CO                         = 0.2000
     H:  4 n_CH4,0 + 2 n_H2O,0 = 2 y_H2O + 2 y_H2        → n_H2O,0 = 0.4000
     O:  n_H2O,0 = 2 y_CO2 + y_CO + y_H2O = 0.4000        (comprobación)
   Avances: ξ1 = CH4 consumido = 0.2000;  ξ2 = CO2 formado = 0.0275
2) Balance de energía (Q = ΔH; W_s = 0; trayectoria por 298.15 K):
     Q = ΔH°_R + ξ1·ΔH°298,1 + ξ2·ΔH°298,2 + ΔH°_P
     ΔH°_R = −Σ n_i,0 ∫_{298.15}^{773.15} Cp_i dT   (reactivos 500 °C → 25 °C)
     ΔH°_P = +Σ n_i  ∫_{298.15}^{1123.15} Cp_i dT   (productos 25 °C → 850 °C)
     ΔH°298,1 = −110 525 − (−74 520) − (−241 818) = 205 813 J
     ΔH°298,2 = −393 509 − (−110 525) − (−241 818) = −41 166 J
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz

from core.energia import (C_a_K, alimentacion_por_atomos, avances_desde_flujos,
                          balance_energia)
from core.equilibrium import format_table

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
REACCIONES = ["CH4 + H2O = CO + 3 H2", "CO + H2O = CO2 + H2"]
SALIDA = {"CO2": 0.0275, "CO": 0.1725, "H2O": 0.1725, "H2": 0.6275}  # base 1 mol
ESPECIES_ALIM = ["CH4", "H2O"]
T_IN, T_OUT = C_a_K(500.0), C_a_K(850.0)
CLAVE_Q = 54881.0                                                    # J
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

alim = alimentacion_por_atomos(SALIDA, ESPECIES_ALIM)
xi, res_xi = avances_desde_flujos(alim, SALIDA, REACCIONES)
b = balance_energia([(alim, T_IN)], SALIDA, T_OUT, REACCIONES, xi)

print("=" * 72)
print("  PROBLEMA 4 — Reformado de CH4 (2 reacciones), Q por mol de producto")
print("=" * 72)
print("\n1) BALANCE DE MATERIA (base 1 mol de producto)")
print(f"   Alimentación: CH4 = {alim['CH4']:.4f} mol, H2O = {alim['H2O']:.4f}"
      f" mol  (relación H2O/CH4 = {alim['H2O'] / alim['CH4']:.2f})")
print(f"   ξ1 = {xi[0]:.4f}   ξ2 = {xi[1]:.4f}   (residual {res_xi:.1e})")
print(f"   Balance atómico (sal − ent): {b['balance_atomos']}")

print("\n2) BALANCE DE ENERGÍA   Q = ΔH°_R + Σξ·ΔH°298 + ΔH°_P")
print(f"   ΔH°298,1 = {b['dH298_j'][0]:.0f} J   ΔH°298,2 = "
      f"{b['dH298_j'][1]:.0f} J")
filas = [["ΔH°_R (reactivos 773.15 → 298.15 K)", b["dHR"]],
         ["ξ1·ΔH°298,1 + ξ2·ΔH°298,2", b["dH298"]],
         ["ΔH°_P (productos 298.15 → 1123.15 K)", b["dHP"]],
         ["Q = ΔH", b["Q"]],
         ["Q por entalpías de formación (comprobación)", b["Q_formacion"]]]
print("  " + format_table(filas, ["término", "J/mol de producto"]
                          ).replace("\n", "\n  "))
print(f"   ⟨Cp⟩_H productos = {b['Cp_medio_P']:.3f} J/K")
print(f"\n   → Q = {b['Q']:.0f} J por mol de gas producido   "
      f"(maestra {CLAVE_Q:.0f} J)")

print("\nAUTOVALIDACIÓN")
ok = True
for nom, cond in (
        ("balance atómico cierra", max(abs(v) for v in
                                     b["balance_atomos"].values()) < 1e-12),
        ("ξ reproduce la salida (n_sal = n_ent + νξ)", b["err_materia"] < 1e-12),
        ("Q trayectoria = Q formación", abs(b["Q"] - b["Q_formacion"]) < 1e-6),
        ("Q vs clave (rtol 1e-4)", abs(b["Q"] - CLAVE_Q) / CLAVE_Q < 1e-4),
        ("Σ y_salida = 1", abs(sum(SALIDA.values()) - 1) < 1e-12)):
    ok &= cond
    print(f"  [{'PASS' if cond else 'FAIL'}] {nom}")
print(f"\nVEREDICTO: {'OK' if ok else 'REVISAR'}")

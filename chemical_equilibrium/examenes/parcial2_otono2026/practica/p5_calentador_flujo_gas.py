"""Ejercicios adicionales 2do parcial — PROBLEMA 5 (calentador de agua,
combustión INCOMPLETA, flujo volumétrico de combustible)

Combustible: 94 % CH4, 6 % N2; 35 % de exceso de aire; combustible y aire
entran secos a 77 (°F). Agua: 75 lbm/s de 77 a 203 (°F). Gases salen a
410 (°F). Del CH4: 70 % → CO2, 30 % → CO. Sin pérdidas. ¿Flujo volumétrico
de combustible?

Planteamiento (a mano) — base 1 mol de combustible
---------------------------------------------------
  R1: CH4 + 2 O2 → CO2 + 2 H2O(g)        ξ1 = 0.70·0.94 = 0.658
  R2: CH4 + 3/2 O2 → CO + 2 H2O(g)       ξ2 = 0.30·0.94 = 0.282
  O2 teórico = 2·0.94 = 1.88 (el exceso se define sobre combustión COMPLETA)
  O2 alimentado = 1.35·1.88 = 2.538;  N2 del aire = 2.538·79/21 = 9.5477
  Salida: CO2 0.658, CO 0.282, H2O 1.88, O2 = 2.538 − 2ξ1 − 1.5ξ2 = 0.799,
          N2 = 9.5477 + 0.06 = 9.6077
  Entradas a 77 °F = 298.15 K → ΔH°_R = 0
  ΔH (gases) = ξ1·ΔH°298,1 + ξ2·ΔH°298,2 + ΔH°_P(298.15 → 483.15 K)  < 0
  Calor al agua = −ΔH(gases) por mol de combustible   (sin pérdidas)
  Q_agua = ṁ·(h_203°F − h_77°F)   (tablas de vapor, líquido saturado)
  ṅ_comb = Q_agua / (−ΔH) ;   V̇ = ṅ R T / P   a 77 (°F) y 1 atm
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz

from core.energia import (ATM_PA, F_a_K, LBM_A_KG, balance_energia,
                          h_liquido_sat, salida_desde_avances)
from core.equilibrium import format_table
from core.thermo import R_J

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
COMB = {"CH4": 0.94, "N2": 0.06}       # base 1 mol de combustible
EXCESO = 0.35
FRAC_CO2 = 0.70                         # resto → CO
T_GAS_IN = F_a_K(77.0)                  # combustible y aire
T_GAS_OUT = F_a_K(410.0)
M_AGUA = 75.0 * LBM_A_KG                # kg/s
T_AGUA = (77.0, 203.0)                  # °F
P_REF, T_REF_V = ATM_PA, F_a_K(77.0)    # condiciones del flujo volumétrico
CLAVE_V = 0.407                         # m³/s
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

REACCIONES = ["CH4 + 2 O2 = CO2 + 2 H2O", "CH4 + 3/2 O2 = CO + 2 H2O"]
nCH4 = COMB["CH4"]
xi = [FRAC_CO2 * nCH4, (1 - FRAC_CO2) * nCH4]
O2 = (1 + EXCESO) * 2 * nCH4
aire = {"O2": O2, "N2": O2 * 79 / 21}
ent_total = {"CH4": nCH4, "O2": O2, "N2": COMB["N2"] + aire["N2"]}
salida = salida_desde_avances(ent_total, REACCIONES, xi)
salida = {k: v for k, v in salida.items() if abs(v) > 1e-14}

b = balance_energia([(COMB, T_GAS_IN), (aire, T_GAS_IN)], salida, T_GAS_OUT,
                    REACCIONES, xi)
q_mol = -b["Q"]                                     # J que recibe el agua por mol
tC = [(t - 32) / 1.8 for t in T_AGUA]
dh = h_liquido_sat(tC[1]) - h_liquido_sat(tC[0])    # kJ/kg
Q_agua = M_AGUA * dh * 1000.0                       # W
n_comb = Q_agua / q_mol                             # mol/s
V = n_comb * R_J * T_REF_V / P_REF                  # m³/s
# sensibilidad: Cp del agua constante 4.184 kJ/(kg·K)
V_cp = M_AGUA * 4.184e3 * (tC[1] - tC[0]) / q_mol * R_J * T_REF_V / P_REF

print("=" * 72)
print("  PROBLEMA 5 — Calentador de agua: flujo volumétrico de combustible")
print("=" * 72)
print("\n1) BALANCE DE MATERIA (base 1 mol de combustible)")
print(f"   ξ1 (→CO2) = {xi[0]:.4f}   ξ2 (→CO) = {xi[1]:.4f}")
print(f"   O2 alimentado = 1.35·2·0.94 = {O2:.4f}   N2 aire = {aire['N2']:.4f}")
print("  " + format_table([[sp, v] for sp, v in salida.items()],
                          ["salida", "mol"]).replace("\n", "\n  "))
print(f"   Balance atómico: {b['balance_atomos']}")

print("\n2) BALANCE DE ENERGÍA DEL LADO DE LOS GASES (por mol de combustible)")
filas = [["ΔH°_R (entradas a 298.15 K)", b["dHR"]],
         ["ΔH°298,1 = " + f"{b['dH298_j'][0]:.0f} J", xi[0] * b["dH298_j"][0]],
         ["ΔH°298,2 = " + f"{b['dH298_j'][1]:.0f} J", xi[1] * b["dH298_j"][1]],
         [f"ΔH°_P (298.15 → {T_GAS_OUT:.2f} K)", b["dHP"]],
         ["ΔH gases = Q", b["Q"]],
         ["comprobación (entalpías de formación)", b["Q_formacion"]]]
print("  " + format_table(filas, ["término", "J/mol comb."]).replace("\n", "\n  "))
print(f"   Calor transferido al agua = −ΔH = {q_mol:.0f} J/mol de combustible")

print("\n3) LADO DEL AGUA")
print(f"   ṁ = 75 lbm/s = {M_AGUA:.4f} kg/s;  T: {tC[0]:.2f} → {tC[1]:.2f} °C")
print(f"   h_l(25 °C) = {h_liquido_sat(tC[0]):.2f}, h_l(95 °C) = "
      f"{h_liquido_sat(tC[1]):.2f} kJ/kg  → Δh = {dh:.2f} kJ/kg")
print(f"   Q_agua = {Q_agua / 1000:.1f} kW")

print("\n4) FLUJO DE COMBUSTIBLE")
print(f"   ṅ = Q_agua/(−ΔH) = {n_comb:.4f} mol/s")
print(f"   V̇ = ṅRT/P (77 °F, 1 atm) = {V:.4f} m³/s = {V * 35.3147:.2f} ft³/s"
      f"   (maestra {CLAVE_V} m³/s)")
print(f"   [sensibilidad] con Cp_agua = 4.184 kJ/kg·K: V̇ = {V_cp:.4f} m³/s")

print("\nAUTOVALIDACIÓN")
ok = True
for nom, cond in (
        ("balance atómico cierra", max(abs(v) for v in
                                     b["balance_atomos"].values()) < 1e-12),
        ("Q trayectoria = Q formación", abs(b["Q"] - b["Q_formacion"]) < 1e-6),
        ("O2 sobrante > 0", salida["O2"] > 0),
        ("V̇ vs clave (redondeo a 3 cifras)", round(V, 3) == CLAVE_V)):
    ok &= cond
    print(f"  [{'PASS' if cond else 'FAIL'}] {nom}")
print(f"\nVEREDICTO: {'OK' if ok else 'REVISAR'}")

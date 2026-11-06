"""2do parcial práctico (Grupo B, 22/10/2026) — PROBLEMA 2
Combustible gaseoso 60 % CH4 / 40 % H2 (mol) con aire estequiométrico

    Base: 1 mol de combustible a 298.15 K; aire teórico (21 % O2 / 79 % N2)
    a 298.15 K; combustión completa; H2O como VAPOR en los productos.

a) O2 y aire teóricos
b) HHV y LHV de la mezcla (kJ/mol y MJ/kg)
c) Temperatura de flama adiabática

Planteamiento (a mano)
----------------------
R1: CH4 + 2 O2 → CO2 + 2 H2O        ξ1 = 0.6
R2: H2 + 1/2 O2 → H2O               ξ2 = 0.4
a) O2 = 0.6·2 + 0.4·0.5 = 1.40 mol;  aire = 1.40/0.21 = 6.667 mol;
   N2 = 1.40·79/21 = 5.267 mol
b) HHV_i = −ΔH°c,i con H2O(l);  LHV_i = −ΔH°c,i con H2O(g)
   HHV_mezcla = Σ y_i·HHV_i  (ídem LHV);   HHV − LHV = n_H2O·ΔH_lv(298)
   n_H2O = 0.6·2 + 0.4·1 = 1.6 mol  →  HHV − LHV = 1.6·44.012 kJ
   MJ/kg: dividir entre M_mezcla = 0.6·M_CH4 + 0.4·M_H2
c) Productos: CO2 0.6, H2O 1.6, N2 5.267 (sin O2: aire teórico)
   Q = 0 = ΔH°_R + Σ ξ_j ΔH°298,j + ΔH°_P,   ΔH°_R = 0 (todo entra a 298 K)
   ⇒ ΔH°_P = Σ n_i R ∫_{298.15}^{T} (Cp_i/R) dT = −Σ ξ_j ΔH°298,j = LHV
   ⇒ brentq sobre T;  T − 298.15 = LHV/⟨Cp⟩_H (comprobación)
"""

import sys, pathlib, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz

from core.energia import (DHV_H2O_298, T_salida_para_Q, balance_energia,
                          combustion_mezcla, masa_molar, salida_desde_avances)
from core.equilibrium import format_table
from core.reactions import parse_reaction

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
COMBUSTIBLE = {"CH4": 0.6, "H2": 0.4}          # fracciones mol, 1 mol total
T_COMB, T_AIRE = 298.15, 298.15                # K
Y_O2_AIRE = 0.21                               # aire 21/79
EXCESO = 0.0                                   # aire estequiométrico
REACCIONES = ["CH4 + 2 O2 = CO2 + 2 H2O", "H2 + 1/2 O2 = H2O"]
ESPERADO = {"O2": 1.40, "aire": 6.667, "N2": 5.267, "HHV": 648.72,
            "LHV": 578.30, "HHV_kg": 62.18, "LHV_kg": 55.43, "T_ad": 2340.0,
            "Cp_H": 283.2}
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

print("=" * 72)
print("  PROBLEMA 2 — Combustible 0.6 CH4 / 0.4 H2 con aire estequiométrico")
print("=" * 72)

# ---------------------------------------------------------------- a)
xi = [COMBUSTIBLE["CH4"], COMBUSTIBLE["H2"]]   # combustión completa
nus = [parse_reaction(r) for r in REACCIONES]
O2_teo = -sum(x * nu["O2"] for x, nu in zip(xi, nus))
O2 = O2_teo * (1 + EXCESO)
aire = O2 / Y_O2_AIRE
N2 = aire - O2
print("\na) AIRE TEÓRICO  (base 1 mol de combustible)")
for r, x, nu in zip(REACCIONES, xi, nus):
    print(f"   {r:28s} ξ = {x:.1f}  → O2 = {x:.1f}·{-nu['O2']:g} = "
          f"{-x * nu['O2']:.2f} mol")
print(f"   O2 teórico = {O2_teo:.4f} mol")
print(f"   aire = O2/0.21 = {aire:.4f} mol;   N2 = O2·79/21 = {N2:.4f} mol")

# ---------------------------------------------------------------- b)
m = combustion_mezcla(COMBUSTIBLE)
print("\nb) PODER CALORÍFICO (ΔH°f de Tabla C.4; H2O(l) = −285 830 J/mol)")
filas = [[sp, COMBUSTIBLE[sp], masa_molar(sp), r["n_H2O"], r["HHV"] / 1000,
          r["LHV"] / 1000] for sp, r in m["por_especie"].items()]
print("  " + format_table(filas, ["especie", "y_i", "M g/mol", "n_H2O",
                                  "HHV kJ/mol", "LHV kJ/mol"]
                          ).replace("\n", "\n  "))
nH2O_tot = sum(COMBUSTIBLE[sp] * r["n_H2O"] for sp, r in m["por_especie"].items())
print(f"   HHV = Σ y_i·HHV_i = {m['HHV'] / 1000:.2f} kJ/mol de mezcla"
      f"   = {m['HHV_kg']:.2f} MJ/kg")
print(f"   LHV = Σ y_i·LHV_i = {m['LHV'] / 1000:.2f} kJ/mol de mezcla"
      f"   = {m['LHV_kg']:.2f} MJ/kg")
print(f"   M mezcla = {m['M']:.4f} g/mol")
print(f"   HHV − LHV = {(m['HHV'] - m['LHV']) / 1000:.3f} kJ = n_H2O·ΔH_lv = "
      f"{nH2O_tot:.1f}·{DHV_H2O_298 / 1000:.3f} = "
      f"{nH2O_tot * DHV_H2O_298 / 1000:.3f} kJ")

# ---------------------------------------------------------------- c)
entradas = [(dict(COMBUSTIBLE), T_COMB), ({"O2": O2, "N2": N2}, T_AIRE)]
ent_total = {**COMBUSTIBLE, "O2": O2, "N2": N2}
salida = {k: v for k, v in salida_desde_avances(ent_total, REACCIONES,
                                                xi).items() if abs(v) > 1e-12}
T_ad = T_salida_para_Q(entradas, salida, REACCIONES, xi, Q=0.0)
with warnings.catch_warnings(record=True) as avisos:
    warnings.simplefilter("always")
    b = balance_energia(entradas, salida, T_ad, REACCIONES, xi)
fuera = [sp for sp in salida if T_ad > b["datos"][sp]["Tmax"]]

print("\nc) TEMPERATURA DE FLAMA ADIABÁTICA  (Q = 0)")
print("   Salida (mol): " + ", ".join(f"{k} {v:.4f}" for k, v in salida.items()))
filas = [["ΔH°_R (todo entra a 298.15 K)", b["dHR"]]]
filas += [[f"ξ{j + 1}·ΔH°298,{j + 1} = {x:.1f}·({h:.0f})", x * h]
          for j, (x, h) in enumerate(zip(xi, b["dH298_j"]))]
filas += [[f"ΔH°_P (298.15 → {T_ad:.2f} K)", b["dHP"]],
          ["Q = ΔH (trayectoria)", b["Q"]],
          ["Q por entalpías de formación", b["Q_formacion"]]]
print("  " + format_table(filas, ["término", "J"]).replace("\n", "\n  "))
print(f"   ⟨Cp⟩_H productos = {b['Cp_medio_P']:.2f} J/K   →   T − 298.15 = "
      f"{-b['dH298']:.0f}/{b['Cp_medio_P']:.2f} = "
      f"{-b['dH298'] / b['Cp_medio_P']:.2f} K")
print(f"   → T_ad = {T_ad:.2f} K  ({T_ad - 273.15:.2f} °C)")
if fuera:
    Tmax = min(b["datos"][sp]["Tmax"] for sp in fuera)
    print(f"   AVISO: T_ad > Tmax = {Tmax:.0f} K del Cp de "
          f"{', '.join(fuera)} (Tabla C.1) → Cp extrapolado")
for w in avisos:
    print(f"   [warning core.energia] {w.message}")

# ---------------------------------------------------------------- checks
print("\nAUTOVALIDACIÓN")
ok = True


def chk(cond, txt):
    global ok
    ok &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {txt}")


chk(abs(O2 - ESPERADO["O2"]) < 1e-9 and abs(aire - ESPERADO["aire"]) < 5e-4
    and abs(N2 - ESPERADO["N2"]) < 5e-4,
    f"a) O2 = {O2:.3f}, aire = {aire:.3f}, N2 = {N2:.3f} mol "
    f"(esperado {ESPERADO['O2']}, {ESPERADO['aire']}, {ESPERADO['N2']})")
chk(abs(m["HHV"] / 1000 - ESPERADO["HHV"]) < 0.005
    and abs(m["LHV"] / 1000 - ESPERADO["LHV"]) < 0.005,
    f"b) HHV = {m['HHV'] / 1000:.2f}, LHV = {m['LHV'] / 1000:.2f} kJ/mol "
    f"(esperado {ESPERADO['HHV']}, {ESPERADO['LHV']})")
chk(abs(m["HHV_kg"] - ESPERADO["HHV_kg"]) < 0.005
    and abs(m["LHV_kg"] - ESPERADO["LHV_kg"]) < 0.005,
    f"b) HHV = {m['HHV_kg']:.2f}, LHV = {m['LHV_kg']:.2f} MJ/kg "
    f"(esperado {ESPERADO['HHV_kg']}, {ESPERADO['LHV_kg']})")
chk(abs((m["HHV"] - m["LHV"]) - 1.6 * 44012.0) < 1e-6,
    f"b) HHV − LHV = {m['HHV'] - m['LHV']:.3f} J = 1.6·44 012 = "
    f"{1.6 * 44012.0:.3f} J")
chk(max(abs(v) for v in b["balance_atomos"].values()) < 1e-12
    and b["err_materia"] < 1e-12,
    f"c) balance atómico cierra {b['balance_atomos']} y ξ reproduce la salida")
chk(abs(salida["CO2"] - 0.6) < 1e-12 and abs(salida["H2O"] - 1.6) < 1e-12
    and abs(salida["N2"] - N2) < 1e-12 and "O2" not in salida,
    "c) productos = CO2 0.6, H2O 1.6, N2 5.267, sin O2")
chk(abs(b["Q"]) < 1e-6 and abs(b["Q_formacion"]) < 1e-6,
    f"c) Q trayectoria = {b['Q']:.1e} J = Q formación = "
    f"{b['Q_formacion']:.1e} J = 0")
chk(abs(-b["dH298"] - m["LHV"]) < 1e-6,
    f"c) −Σξ·ΔH°298 = {-b['dH298']:.1f} J = LHV de la mezcla")
chk(abs(T_ad - ESPERADO["T_ad"]) < 1.0
    and abs(b["Cp_medio_P"] - ESPERADO["Cp_H"]) < 0.1,
    f"c) T_ad = {T_ad:.2f} K, ⟨Cp⟩_H = {b['Cp_medio_P']:.2f} J/K "
    f"(esperado ≈ {ESPERADO['T_ad']:.0f}, ≈ {ESPERADO['Cp_H']})")
print(f"  [info] Cp extrapolado (T > Tmax): {', '.join(fuera) if fuera else 'no'}")

print("\n>>> RESULTADOS")
print(f"    a) O2 = {O2:.2f} mol, aire = {aire:.3f} mol, N2 = {N2:.3f} mol")
print(f"    b) HHV = {m['HHV'] / 1000:.2f} kJ/mol ({m['HHV_kg']:.2f} MJ/kg);"
      f"  LHV = {m['LHV'] / 1000:.2f} kJ/mol ({m['LHV_kg']:.2f} MJ/kg)")
print(f"    c) T_ad = {T_ad:.1f} K ({T_ad - 273.15:.1f} °C)")
print(f"\nVEREDICTO: {'OK' if ok else 'REVISAR'}")

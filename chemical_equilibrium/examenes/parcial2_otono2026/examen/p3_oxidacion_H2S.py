"""2do parcial práctico (Grupo B, 22/10/2026) — PROBLEMA 3
Oxidación de H2S con aire (dos reacciones), proceso NO adiabático

    R1: H2S + 3/2 O2 → SO2 + H2O          (balanceadas aquí; en el examen
    R2: SO2 + 1/2 O2 → SO3                 vienen sin coeficientes)

Entra SOLO H2S y aire (21 % O2 / 79 % N2) a 100 °C; sale a 830 °C con
    H2S 0.2, O2 1.65, SO2 0.5, H2O 0.8 mol, más SO3 y N2 desconocidos.

a) Alimentación y salida completas
b) Q del proceso
c) T de salida si fuera adiabático (misma composición de salida)

Planteamiento (a mano)
----------------------
a) Balance por avances  n_sal = n_ent + Σ ν_ij ξ_j:
     H2O:  0.8 = 0 + ξ1                 → ξ1 = 0.8  (no entra agua)
     SO2:  0.5 = 0 + ξ1 − ξ2            → ξ2 = ξ1 − 0.5 = 0.3
     SO3:  n_SO3 = ξ2 = 0.3
     H2S:  0.2 = n_H2S,0 − ξ1           → n_H2S,0 = 1.0
     O2:   1.65 = n_O2,0 − 3/2 ξ1 − 1/2 ξ2 → n_O2,0 = 1.65 + 1.2 + 0.15 = 3.0
     N2:   inerte, n_N2 = n_O2,0·79/21 = 11.286
   Base de cálculo: 1 mol de H2S alimentado.
b) Q = ΔH (W_s = 0, ΔEc ≈ ΔEp ≈ 0), trayectoria por 298.15 K:
     Q = ΔH°_R + ξ1·ΔH°298,1 + ξ2·ΔH°298,2 + ΔH°_P
     ΔH°_R = −Σ n_i,0 ∫_{298.15}^{373.15} Cp_i dT    (alimentación 100 → 25 °C)
     ΔH°_P = +Σ n_i  ∫_{298.15}^{1103.15} Cp_i dT    (productos 25 → 830 °C)
     ΔH°298,1 = −296 830 − 241 818 − (−20 630) = −518 018 J
     ΔH°298,2 = −395 720 − (−296 830)          =  −98 890 J
c) Q = 0 con las mismas ξ y salida → brentq sobre T_salida.
"""

import sys, pathlib, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz

from core.energia import (C_a_K, T_salida_para_Q, avances_desde_flujos,
                          balance_atomos, balance_energia)
from core.equilibrium import format_table

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
REACCIONES = ["H2S + 3/2 O2 = SO2 + H2O", "SO2 + 1/2 O2 = SO3"]
SALIDA_CONOCIDA = {"H2S": 0.2, "O2": 1.65, "SO2": 0.5, "H2O": 0.8}   # mol
T_IN, T_OUT = C_a_K(100.0), C_a_K(830.0)          # 373.15 K, 1103.15 K
Y_O2_AIRE = 0.21
ESPERADO = {"xi": (0.8, 0.3), "H2S0": 1.0, "O20": 3.0, "N2": 11.286,
            "dH1": -518018.0, "dH2": -98890.0, "dHR": -34112.0,
            "dH298": -444081.0, "dHP": 390515.0, "Q": -87679.0,
            "Q_mol_sal": -5.95, "T_ad": 1269.8}
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

s = SALIDA_CONOCIDA
# ---- a) balance de materia por avances (despejes del planteamiento)
xi1 = s["H2O"]                                   # no entra agua
xi2 = xi1 - s["SO2"]
n_SO3 = xi2
H2S0 = s["H2S"] + xi1
O20 = s["O2"] + 1.5 * xi1 + 0.5 * xi2
N2 = O20 * (1 - Y_O2_AIRE) / Y_O2_AIRE
alim = {"H2S": H2S0, "O2": O20, "N2": N2}
salida = {"H2S": s["H2S"], "O2": s["O2"], "SO2": s["SO2"], "H2O": s["H2O"],
          "SO3": n_SO3, "N2": N2}
xi = [xi1, xi2]
nT_in, nT_out = sum(alim.values()), sum(salida.values())

print("=" * 72)
print("  PROBLEMA 3 — Oxidación de H2S: R1 H2S + 3/2 O2 → SO2 + H2O,")
print("               R2 SO2 + 1/2 O2 → SO3;  100 °C → 830 °C")
print("=" * 72)

print("\na) BALANCE DE MATERIA  (base: la salida dada ≡ 1 mol de H2S alimentado)")
print(f"   H2O: 0.8 = ξ1                      → ξ1 = {xi1:.4f}")
print(f"   SO2: 0.5 = ξ1 − ξ2                 → ξ2 = {xi2:.4f}")
print(f"   SO3: n = ξ2                        → n_SO3 = {n_SO3:.4f}")
print(f"   H2S: 0.2 = n_H2S,0 − ξ1            → n_H2S,0 = {H2S0:.4f}")
print(f"   O2 : 1.65 = n_O2,0 − 1.5ξ1 − 0.5ξ2 → n_O2,0 = 1.65 + "
      f"{1.5 * xi1:.2f} + {0.5 * xi2:.2f} = {O20:.4f}")
print(f"   N2 : n_O2,0·79/21                  → n_N2 = {N2:.4f}  (aire = "
      f"{O20 + N2:.4f} mol)")
esp = ["H2S", "O2", "N2", "SO2", "H2O", "SO3"]
filas = [[sp, alim.get(sp, 0.0), alim.get(sp, 0.0) / nT_in,
          salida.get(sp, 0.0), salida.get(sp, 0.0) / nT_out] for sp in esp]
filas.append(["TOTAL", nT_in, 1.0, nT_out, 1.0])
print("  " + format_table(filas, ["especie", "entra mol", "y_ent", "sale mol",
                                  "y_sal"]).replace("\n", "\n  "))
bal = balance_atomos(alim, salida)
xi_ls, res_ls = avances_desde_flujos(alim, salida, REACCIONES)
print(f"   Balance atómico (sal − ent): {bal}")
print(f"   avances_desde_flujos → ξ = ({xi_ls[0]:.6f}, {xi_ls[1]:.6f}),"
      f" residual {res_ls:.1e}")

# ---- b) balance de energía
with warnings.catch_warnings(record=True) as avisos:
    warnings.simplefilter("always")
    b = balance_energia([(alim, T_IN)], salida, T_OUT, REACCIONES, xi)
print(f"\nb) BALANCE DE ENERGÍA  Q = ΔH°_R + Σξ·ΔH°298 + ΔH°_P   "
      f"(T_in = {T_IN:.2f} K, T_out = {T_OUT:.2f} K)")
print(f"   ΔH°298,1 = {b['dH298_j'][0]:.0f} J   ΔH°298,2 = "
      f"{b['dH298_j'][1]:.0f} J")
filas = [[f"ΔH°_R (alimentación {T_IN:.2f} → 298.15 K)", b["dHR"]],
         [f"ξ1·ΔH°298,1 = {xi1:.1f}·({b['dH298_j'][0]:.0f})",
          xi1 * b["dH298_j"][0]],
         [f"ξ2·ΔH°298,2 = {xi2:.1f}·({b['dH298_j'][1]:.0f})",
          xi2 * b["dH298_j"][1]],
         ["Σ ξ·ΔH°298", b["dH298"]],
         [f"ΔH°_P (productos 298.15 → {T_OUT:.2f} K)", b["dHP"]],
         ["Q = ΔH (trayectoria)", b["Q"]],
         ["Q por entalpías de formación (comprobación)", b["Q_formacion"]]]
print("  " + format_table(filas, ["término", "J"]).replace("\n", "\n  "))
print(f"   ⟨Cp⟩_H productos = {b['Cp_medio_P']:.3f} J/K")
Q_H2S = b["Q"] / H2S0
Q_sal = b["Q"] / nT_out
print(f"   → Q = {b['Q']:.0f} J = {b['Q'] / 1000:.2f} kJ por {H2S0:g} mol de "
      f"H2S alimentado  ({Q_H2S / 1000:.2f} kJ/mol H2S)")
print(f"     = {Q_sal / 1000:.3f} kJ por mol de gas de salida "
      f"(n_sal = {nT_out:.4f} mol)")
print("   Interpretación: Q < 0 → el sistema LIBERA calor; las reacciones son")
print("   exotérmicas y hay que RETIRAR calor para mantener la salida a 830 °C.")
for w in avisos:
    print(f"   [warning core.energia] {w.message}")

# ---- c) adiabático
T_ad = T_salida_para_Q([(alim, T_IN)], salida, REACCIONES, xi, Q=0.0)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    bc = balance_energia([(alim, T_IN)], salida, T_ad, REACCIONES, xi)
fuera = [sp for sp in salida if T_ad > bc["datos"][sp]["Tmax"]]
print(f"\nc) ADIABÁTICO (Q = 0), misma salida y mismas ξ")
print(f"   ΔH°_P = −(ΔH°_R + Σξ·ΔH°298) = {-(bc['dHR'] + bc['dH298']):.0f} J"
      f"  → brentq sobre T")
print(f"   ⟨Cp⟩_H productos = {bc['Cp_medio_P']:.3f} J/K;  ΔH°_P(T_ad) = "
      f"{bc['dHP']:.0f} J;  Q = {bc['Q']:.1e} J")
print(f"   → T_ad = {T_ad:.2f} K  ({T_ad - 273.15:.2f} °C)  > 1103.15 K: "
      f"sin retirar calor la salida se calienta más")
print(f"   Cp extrapolado (T > Tmax): {', '.join(fuera) if fuera else 'no'}")

# ---------------------------------------------------------------- checks
print("\nAUTOVALIDACIÓN")
ok = True


def chk(cond, txt):
    global ok
    ok &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {txt}")


chk(all(abs(v) < 1e-12 for v in bal.values()),
    f"a) balance de S, H, O, N con balance_atomos: {bal}")
chk(abs(xi_ls[0] - 0.8) < 1e-12 and abs(xi_ls[1] - 0.3) < 1e-12
    and res_ls < 1e-12,
    f"a) avances_desde_flujos da ξ = ({xi_ls[0]:.6f}, {xi_ls[1]:.6f}) = "
    f"(0.8, 0.3)")
chk(abs(H2S0 - ESPERADO["H2S0"]) < 1e-12 and abs(O20 - ESPERADO["O20"]) < 1e-12
    and abs(N2 - ESPERADO["N2"]) < 5e-4,
    f"a) alimentación H2S = {H2S0:.3f}, O2 = {O20:.3f}, N2 = {N2:.3f} "
    f"(esperado 1.0, 3.0, 11.286)")
chk(abs(alim["O2"] / (alim["O2"] + alim["N2"]) - Y_O2_AIRE) < 1e-12,
    "a) aire 21/79: O2/(O2 + N2) = 0.21")
chk(abs(b["dH298_j"][0] - ESPERADO["dH1"]) < 0.5
    and abs(b["dH298_j"][1] - ESPERADO["dH2"]) < 0.5,
    f"b) ΔH°298,1 = {b['dH298_j'][0]:.0f}, ΔH°298,2 = {b['dH298_j'][1]:.0f} J")
for k, v in (("dHR", b["dHR"]), ("dH298", b["dH298"]), ("dHP", b["dHP"]),
             ("Q", b["Q"])):
    chk(abs(v - ESPERADO[k]) / abs(ESPERADO[k]) < 1e-4,
        f"b) {k} = {v:.0f} J vs esperado {ESPERADO[k]:.0f} J")
chk(abs(b["Q"] - b["Q_formacion"]) < 1e-6 and b["err_materia"] < 1e-12,
    f"b) Q trayectoria = Q formación ({b['Q_formacion']:.1f} J) y ξ "
    f"reproduce la salida")
chk(abs(Q_sal / 1000 - ESPERADO["Q_mol_sal"]) < 0.005,
    f"b) Q por mol de salida = {Q_sal / 1000:.3f} kJ/mol (≈ −5.95)")
chk(b["Q"] < 0, "b) Q < 0: proceso exotérmico, hay que retirar calor")
chk(abs(T_ad - ESPERADO["T_ad"]) < 0.1 and abs(bc["Q"]) < 1e-6,
    f"c) T_ad = {T_ad:.2f} K (esperado {ESPERADO['T_ad']}) con Q(T_ad) = "
    f"{bc['Q']:.1e} J")
chk(T_ad > T_OUT, f"c) T_ad = {T_ad:.2f} K > {T_OUT:.2f} K, coherente con "
    f"Q < 0 en b)")

print("\n>>> RESULTADOS")
print(f"    a) ξ1 = {xi1:.1f}, ξ2 = {xi2:.1f}; alimentación: H2S {H2S0:.3f}, "
      f"O2 {O20:.3f}, N2 {N2:.3f} mol; salida: SO3 {n_SO3:.1f}, "
      f"N2 {N2:.3f} mol")
print(f"    b) Q = {b['Q'] / 1000:.2f} kJ por mol de H2S alimentado "
      f"({Q_sal / 1000:.2f} kJ/mol de gas de salida)")
print(f"    c) T_ad = {T_ad:.1f} K ({T_ad - 273.15:.1f} °C)")
print(f"\nVEREDICTO: {'OK' if ok else 'REVISAR'}")

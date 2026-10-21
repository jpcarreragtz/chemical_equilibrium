"""SIMULACRO d) — Q con dos reacciones: producción de syngas (Sesión 19)

    R1: CH4 + H2O → CO + 3 H2        R2: CH4 + 2 H2O → CO2 + 4 H2
    2 mol de vapor por mol de CH4; reactivos precalentados a 600 K;
    productos a 1300 K; CH4 convertido por completo; 17.4 % mol de CO.
    ¿Calor requerido por el reactor?     python3 d_syngas_dos_reacciones.py

PLANTEAMIENTO A MANO   (base 1 mol de CH4, 2 mol de H2O)
--------------------
Sea x la fracción del CH4 que va por R1 → ξ1 = x, ξ2 = 1 − x
    CO = x        CO2 = 1 − x        H2 = 3x + 4(1 − x) = 4 − x
    H2O = 2 − x − 2(1 − x) = x       CH4 = 0
    n_T = x + (1 − x) + (4 − x) + x = 5     (no depende de x)
    y_CO = x/5 = 0.174  →  x = 0.87  →  ξ1 = 0.87 (CO), ξ2 = 0.13 (CO2)
ΔH°298,1 = −110 525 − (−74 520) − (−241 818) = 205 813 J
ΔH°298,2 = −393 509 − (−74 520) − 2(−241 818) = 164 647 J
Q = ΔH°_R + ξ1·ΔH°298,1 + ξ2·ΔH°298,2 + ΔH°_P
    ΔH°_R = −[1·∫Cp_CH4 + 2·∫Cp_H2O] de 298.15 a 600 K
    ΔH°_P = Σ n_sal ∫Cp de 298.15 a 1300 K
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz (simulacro/ está 3 niveles abajo)

from core.energia import (C_a_K, F_a_K, T_salida_para_Q, alimentacion_por_atomos,
                          avances_desde_flujos, balance_energia,
                          combustion_mezcla, salida_desde_avances)
from core.equilibrium import format_table

# ══════════════════════════ INPUTS (edita SOLO esto) ══════════════════════
titulo = "Simulacro d) syngas: 2 reacciones, reactivos 600 K, productos 1300 K, 17.4 % CO"
modo = "Q"                             # "Q" | "T_ad" | "HHV"

reacciones = ["CH4 + H2O = CO + 3 H2",          # R1
              "CH4 + 2 H2O = CO2 + 4 H2"]       # R2

# corrientes de entrada: lista de (moles, T [K]); se enfrían a 298.15 K
entradas = [({"CH4": 1.0, "H2O": 2.0}, 600.0)]
y_CO = 0.174                          # dato: % mol de CO en los productos
x = y_CO * 5.0                        # n_T = 5 para cualquier x (planteamiento)
xi = [x, 1.0 - x]                     # ξ1 = 0.87, ξ2 = 0.13
salida = None                         # None → se calcula con ξ
alimentacion_atomos = None

T_salida = 1300.0                     # modo "Q": T de salida [K]
Q_dado = 0.0                          # modo "T_ad" (no se usa aquí)

mezcla = None                         # modo "HHV" (no se usa aquí)
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

# ─── A partir de aquí: COPIA LITERAL de plantilla/plantilla_balance_energia.py (no editar) ───

print("=" * 72 + f"\n  {titulo}\n" + "=" * 72)

if modo == "HHV":
    m = combustion_mezcla(mezcla)
    filas = [[sp, mezcla[sp], r["HHV"] / 1000, r["LHV"] / 1000]
             for sp, r in m["por_especie"].items()]
    print("  " + format_table(filas, ["especie", "y_i", "HHV kJ/mol",
                                      "LHV kJ/mol"]).replace("\n", "\n  "))
    print(f"\n  HHV = {m['HHV'] / 1000:.3f} kJ/mol  ({m['HHV_kg']:.3f} MJ/kg)")
    print(f"  LHV = {m['LHV'] / 1000:.3f} kJ/mol  ({m['LHV_kg']:.3f} MJ/kg)")
    print(f"  ΔH°c (H2O l) = {m['dHc_l'] / 1000:.3f} kJ/mol;  M = {m['M']:.3f}"
          f" g/mol;  Σy = {m['suma_y']:.4f}")
    raise SystemExit

if salida is not None and alimentacion_atomos:
    alim = alimentacion_por_atomos(salida, alimentacion_atomos)
    T_in = entradas[0][1]
    entradas = [(alim, T_in)]
    print(f"  Alimentación por balance atómico: {alim}")

ent_total = {}
for n, _T in entradas:
    for sp, v in n.items():
        ent_total[sp] = ent_total.get(sp, 0.0) + v

if salida is not None and xi is None:
    xi, res = avances_desde_flujos(ent_total, salida, reacciones)
    print(f"  ξ por mínimos cuadrados: {[round(x, 6) for x in xi]} "
          f"(residual {res:.1e})")
if salida is None:
    salida = {k: v for k, v in salida_desde_avances(ent_total, reacciones,
                                                    xi).items() if abs(v) > 1e-14}

if modo == "T_ad":
    T_salida = T_salida_para_Q(entradas, salida, reacciones, xi, Q=Q_dado)

b = balance_energia(entradas, salida, T_salida, reacciones, xi)
print("\n  Salida (mol): " + ", ".join(f"{k} {v:.4f}" for k, v in salida.items()))
print(f"  Balance atómico: {b['balance_atomos']}")
filas = [["ΔH°_R", b["dHR"]]]
filas += [[f"ξ{j + 1}·ΔH°298 ({h:.0f} J)", x * h]
          for j, (x, h) in enumerate(zip(xi, b["dH298_j"]))]
filas += [[f"ΔH°_P (→ {T_salida:.2f} K)", b["dHP"]], ["Q = ΔH", b["Q"]],
          ["Q por entalpías de formación", b["Q_formacion"]]]
print("  " + format_table(filas, ["término", "J"]).replace("\n", "\n  "))
print(f"\n  T salida = {T_salida:.2f} K ({T_salida - 273.15:.2f} °C);  "
      f"Q = {b['Q']:.1f} J;  ⟨Cp⟩_H = {b['Cp_medio_P']:.3f} J/K")
ok = abs(b["Q"] - b["Q_formacion"]) < 1e-6 and b["err_materia"] < 1e-9
print(f"  [{'PASS' if ok else 'FAIL'}] trayectoria = formación y balance de "
      f"materia consistente")

# ══════════════ AUTOVALIDACIÓN DEL SIMULACRO (cálculo independiente) ══════════
_ok = True
def _chk(cond, txt):
    global _ok
    _ok &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {txt}")

_nT = sum(salida.values())
_chk(abs(xi[0] - 0.87) < 1e-12 and abs(xi[1] - 0.13) < 1e-12, "ξ1 = 0.87, ξ2 = 0.13 (a mano)")
_chk(abs(_nT - 5.0) < 1e-12 and abs(salida["CO"] / _nT - y_CO) < 1e-12,
     f"n_T = {_nT:.4f}, y_CO = {salida['CO'] / _nT:.4f} (dato 0.174)")
_chk(abs(salida.get("CH4", 0.0)) < 1e-12, "CH4 convertido por completo")
_chk(max(abs(v) for v in b["balance_atomos"].values()) < 1e-12,
     f"balance atómico cierra: {b['balance_atomos']}")
_chk(abs(b["dH298_j"][0] - 205813) < 0.5 and abs(b["dH298_j"][1] - 164647) < 0.5,
     f"ΔH°298 = {b['dH298_j'][0]:.0f} y {b['dH298_j'][1]:.0f} J (a mano 205 813 y 164 647)")
_chk(abs(b["Q"] - b["Q_formacion"]) < 1e-6,
     f"Q trayectoria = Q formación ({b['Q']:.2f} vs {b['Q_formacion']:.2f} J)")
_chk(abs(b["Q"] / 1000 - 328.0) < 1.0, f"Q = {b['Q'] / 1000:.2f} kJ/mol CH4 (≈ 328)")
print(f"\nRESULTADO d): Q = {b['Q'] / 1000:.2f} kJ por mol de CH4 (se suministra calor)")
print(f"VEREDICTO SIMULACRO: {'OK' if _ok else 'REVISAR'}")

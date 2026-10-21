"""SIMULACRO b) — Temperatura teórica de flama del etileno (Sesión 18, Ej. 2)

    Etileno a 25 °C, combustión completa, proceso adiabático (Q = 0).
    Incisos (python3 b_flama_etileno.py a|b|c|d|e):
      a) aire teórico a 25 °C        b) 25 % de exceso a 25 °C
      c) 50 % de exceso a 25 °C      d) 100 % de exceso a 25 °C
      e) 50 % de exceso, aire precalentado a 500 °C
    Referencia: excel/s18_ej2_flama_adiabatica_etileno.xlsx (fila 62).

PLANTEAMIENTO A MANO   (base 1 mol de C2H4)
--------------------
R:  C2H4 + 3 O2 → 2 CO2 + 2 H2O(g)       ξ = 1 (combustión completa)
    ΔH°298 = 2(−393 509) + 2(−241 818) − 52 510 = −1 323 164 J
Aire: O2 = 3(1 + exceso);  N2 = O2·79/21
Salida: CO2 = 2, H2O = 2, O2 = 3·exceso, N2 = O2·79/21
Balance de energía (W_s = 0, Q = 0), trayectoria por 298.15 K:
    0 = ΔH°_R + ξ·ΔH°298 + ΔH°_P
    ΔH°_R = −Σ n_aire R·∫_{298.15}^{T_aire} (Cp/R) dT   (0 si el aire entra a 25 °C)
    ΔH°_P = Σ n_sal R·∫_{298.15}^{T} (Cp/R) dT    →  despejar T (brentq)
Más exceso de aire → más moles que calentar → menor T; precalentar la sube.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz (simulacro/ está 3 niveles abajo)

from core.energia import (C_a_K, F_a_K, T_salida_para_Q, alimentacion_por_atomos,
                          avances_desde_flujos, balance_energia,
                          combustion_mezcla, salida_desde_avances)
from core.equilibrium import format_table

# ══════════════════════════ INPUTS (edita SOLO esto) ══════════════════════
INCISO = sys.argv[1] if len(sys.argv) > 1 else "a"
#           exceso, T del aire [°C], T del Excel s18 [K]
CASOS = {"a": (0.00, 25.0, 2533.4878), "b": (0.25, 25.0, 2198.5524),
         "c": (0.50, 25.0, 1950.8989), "d": (1.00, 25.0, 1609.1945),
         "e": (0.50, 500.0, 2282.5381)}
exceso, T_aire_C, T_EXCEL = CASOS[INCISO]
titulo = (f"Simulacro b{INCISO}) flama adiabática C2H4: {exceso:.0%} de exceso, "
          f"aire a {T_aire_C:g} °C")
modo = "T_ad"                          # "Q" | "T_ad" | "HHV"

reacciones = ["C2H4 + 3 O2 = 2 CO2 + 2 H2O"]

# corrientes de entrada: lista de (moles, T [K]); se enfrían a 298.15 K
O2 = 3.0 * (1.0 + exceso)
entradas = [({"C2H4": 1.0}, C_a_K(25)),
            ({"O2": O2, "N2": O2 * 79 / 21}, C_a_K(T_aire_C))]
xi = [1.0]                            # avances (mismo orden que reacciones)
salida = None                         # None → se calcula con ξ
alimentacion_atomos = None

T_salida = None                       # modo "Q": T de salida [K]
Q_dado = 0.0                          # modo "T_ad": Q [J] (0 = adiabático)

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

_chk(max(abs(v) for v in b["balance_atomos"].values()) < 1e-12,
     f"balance atómico cierra: {b['balance_atomos']}")
_chk(abs(b["Q"] - b["Q_formacion"]) < 1e-6 and abs(b["Q"]) < 1e-6,
     f"Q trayectoria = Q formación = 0 ({b['Q']:.1e} J, {b['Q_formacion']:.1e} J)")
_chk(abs(b["dH298"] - (-1323164.0)) < 0.5, f"ξ·ΔH°298 = {b['dH298']:.0f} J (a mano −1 323 164)")
_chk(abs(T_salida - T_EXCEL) < 0.05,
     f"T = {T_salida:.2f} K vs Excel s18 {T_EXCEL:.2f} K (|Δ| = {abs(T_salida - T_EXCEL):.1e} K)")
fuera = [sp for sp in salida if T_salida > b["datos"][sp]["Tmax"]]
print(f"  [info] Cp extrapolado (T > Tmax): {', '.join(fuera) if fuera else 'no'}")
print(f"\nRESULTADO b) inciso {INCISO}: T flama = {T_salida:.1f} K ({T_salida - 273.15:.1f} °C)")
print(f"VEREDICTO SIMULACRO: {'OK' if _ok else 'REVISAR'}")

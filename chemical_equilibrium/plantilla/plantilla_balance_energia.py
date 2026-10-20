"""PLANTILLA — Balance de energía con reacción / combustión / flama adiabática
(2do parcial: Sesiones 16–19).  Copia, edita SOLO el bloque INPUTS y corre:

    cp plantilla/plantilla_balance_energia.py ejercicios/ejN_tema.py
    cd ejercicios && python3 ejN_tema.py

Tres modos (elige uno en `modo`):
  "Q"      → conocidas entradas (con su T), avances ξ y T de salida: calcula
             Q = ΔH°_R + Σ ξ ΔH°298 + ΔH°_P  (y lo comprueba por entalpías
             de formación).  Ej.: Problema 4, syngas Sesión 19.
  "T_ad"   → Q dado (0 = adiabático): calcula la T de salida.
             Ej.: flama adiabática (Sesión 18), NH3→NO (Sesión 19 Ej. 3).
  "HHV"    → poder calorífico de una mezcla combustible (Σ y_i HHV_i).
             Ej.: Problema 3, Sesión 16–17.
Si en vez de ξ te dan la COMPOSICIÓN DE SALIDA (Problema 4), pon
`salida = {...}` y `xi = None`: los ξ salen por mínimos cuadrados y la
alimentación (si `alimentacion_por_atomos = [...]`) por balance atómico.

Reglas
* Especies por FÓRMULA (CH4, C2H6, H2O, CO2...). Datos Cp/ΔH°f del JSON
  verificado (data/sva_tables.json); el H2O de los balances es GAS.
* El exceso de aire se define sobre la combustión COMPLETA (a CO2).
* Temperaturas en K: usa C_a_K(°C) o F_a_K(°F).
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # raíz

from core.energia import (C_a_K, F_a_K, T_salida_para_Q, alimentacion_por_atomos,
                          avances_desde_flujos, balance_energia,
                          combustion_mezcla, salida_desde_avances)
from core.equilibrium import format_table

# ══════════════════════════ INPUTS (edita SOLO esto) ══════════════════════
titulo = "Flama adiabática: CH4 con 20 % de exceso de aire (SVA ej. 4.7)"
modo = "T_ad"                          # "Q" | "T_ad" | "HHV"

reacciones = ["CH4 + 2 O2 = CO2 + 2 H2O"]

# corrientes de entrada: lista de (moles, T [K]); se enfrían a 298.15 K
O2 = 1.2 * 2.0
entradas = [({"CH4": 1.0}, C_a_K(25)),
            ({"O2": O2, "N2": O2 * 79 / 21}, C_a_K(25))]
xi = [1.0]                            # avances (mismo orden que reacciones)
salida = None                         # None → se calcula con ξ
alimentacion_atomos = None            # p. ej. ["CH4", "H2O"] (Problema 4)

T_salida = None                       # modo "Q": T de salida [K]
Q_dado = 0.0                          # modo "T_ad": Q [J] (0 = adiabático)

mezcla = {"CH4": 0.875, "C2H6": 0.07, "C3H8": 0.02, "N2": 0.035}  # modo "HHV"
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

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

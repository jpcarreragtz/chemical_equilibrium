"""SIMULACRO e) — HHV / LHV de mezclas combustibles (Sesiones 16 y 17)

    python3 e_hhv_lhv_mezclas.py 1   → Sesión 17 Ej. 2: 87.5 % CH4, 7.0 % C2H6,
                                        2.0 % C3H8, 3.5 % N2: HHV y LHV en kJ/mol y kJ/kg
    python3 e_hhv_lhv_mezclas.py 2   → Sesión 16 Ej. 1: 75 % CH4, 10 % C2H6,
                                        12 % H2, 3 % N2: O2 y aire estequiométricos
                                        para 100 mol; HHV y LHV en kJ/mol y MJ/kg

PLANTEAMIENTO A MANO
--------------------
Combustión completa (1 mol de cada combustible):
    CH4  + 2 O2   → CO2 + 2 H2O         HHV = 890.649   LHV = 802.625 kJ/mol
    C2H6 + 7/2 O2 → 2 CO2 + 3 H2O       HHV = 1560.688  LHV = 1428.652
    C3H8 + 5 O2   → 3 CO2 + 4 H2O       HHV = 2219.167  LHV = 2043.119
    H2   + 1/2 O2 → H2O                 HHV = 285.830   LHV = 241.818
  (ley de Hess con ΔH°f de la Tabla C.4; HHV con H2O(l) = −285 830 J/mol,
   LHV con H2O(g) = −241 818; HHV = LHV + n_H2O·44.012 kJ/mol)
Mezcla: HHV = Σ y_i·HHV_i (N2 no aporta);  M = Σ y_i·M_i
    kJ/kg = (kJ/mol)/(M/1000);   MJ/kg = kJ/g = (kJ/mol)/M
O2 estequiométrico (100 mol de mezcla) = 100·Σ y_i·ν_O2,i;  aire = O2/0.21
  Sesión 16: 100·(0.75·2 + 0.10·3.5 + 0.12·0.5) = 191 mol O2 → 909.5 mol de aire
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz (simulacro/ está 3 niveles abajo)

from core.energia import (C_a_K, F_a_K, T_salida_para_Q, alimentacion_por_atomos,
                          avances_desde_flujos, balance_energia,
                          combustion_mezcla, salida_desde_avances)
from core.equilibrium import format_table

# ══════════════════════════ INPUTS (edita SOLO esto) ══════════════════════
INCISO = sys.argv[1] if len(sys.argv) > 1 else "1"
CASOS = {
    "1": ("Simulacro e1) Sesión 17 Ej. 2: gas natural 87.5/7/2/3.5",
          {"CH4": 0.875, "C2H6": 0.07, "C3H8": 0.02, "N2": 0.035}),
    "2": ("Simulacro e2) Sesión 16 Ej. 1: 75 % CH4, 10 % C2H6, 12 % H2, 3 % N2",
          {"CH4": 0.75, "C2H6": 0.10, "H2": 0.12, "N2": 0.03}),
}
titulo, mezcla = CASOS[INCISO]
modo = "HHV"                           # "Q" | "T_ad" | "HHV"
# (reacciones/entradas/ξ no se usan en modo "HHV")
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

# ─── A partir de aquí: COPIA LITERAL de plantilla/plantilla_balance_energia.py,
#     solo la rama modo "HHV" (sin el `raise SystemExit` final, para que
#     corra la autovalidación de abajo) ───

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

# ══════════════ AUTOVALIDACIÓN DEL SIMULACRO (cálculo independiente) ══════════
from core.energia import masa_molar, reaccion_combustion
_ok = True
def _chk(cond, txt):
    global _ok
    _ok &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {txt}")

_HHV = {"CH4": 890.649, "C2H6": 1560.688, "C3H8": 2219.167, "H2": 285.830}   # kJ/mol
_nW = {"CH4": 2, "C2H6": 3, "C3H8": 4, "H2": 1}
_hhv = sum(y * _HHV[sp] for sp, y in mezcla.items() if sp in _HHV)
_lhv = _hhv - 44.012 * sum(y * _nW[sp] for sp, y in mezcla.items() if sp in _nW)
_M = m["M"]
print("CHEQUEO A MANO:")
print(f"   HHV = Σ y·HHV_i = {_hhv:.3f} kJ/mol;  LHV = HHV − 44.012·Σ y·n_H2O = {_lhv:.3f} kJ/mol")
_chk(abs(m["HHV"] / 1000 - _hhv) < 0.01, f"HHV plantilla {m['HHV'] / 1000:.3f} = a mano {_hhv:.3f} kJ/mol")
_chk(abs(m["LHV"] / 1000 - _lhv) < 0.01, f"LHV plantilla {m['LHV'] / 1000:.3f} = a mano {_lhv:.3f} kJ/mol")
_chk(abs(m["suma_y"] - 1) < 1e-12, f"Σ y = {m['suma_y']:.4f}")
_chk(m["HHV"] > m["LHV"], "HHV > LHV")
print(f"   M = {_M:.4f} g/mol → HHV = {m['HHV_kg']:.3f} MJ/kg = {m['HHV_kg'] * 1000:.0f} kJ/kg;"
      f"  LHV = {m['LHV_kg']:.3f} MJ/kg = {m['LHV_kg'] * 1000:.0f} kJ/kg")
_O2 = 100 * sum(y * -reaccion_combustion(sp)["O2"]
                for sp, y in mezcla.items() if sp in _HHV)
_O2_mano = 100 * sum(y * {"CH4": 2, "C2H6": 3.5, "C3H8": 5, "H2": 0.5}[sp]
                     for sp, y in mezcla.items() if sp in _HHV)
_chk(abs(_O2 - _O2_mano) < 1e-9, f"O2 estequiométrico (100 mol) = {_O2:.2f} mol (a mano {_O2_mano:.2f})")
print(f"   Aire estequiométrico (100 mol de mezcla) = {_O2:.2f}/0.21 = {_O2 / 0.21:.2f} mol")
print(f"\nRESULTADO e) inciso {INCISO}: HHV = {m['HHV'] / 1000:.2f} kJ/mol ({m['HHV_kg'] * 1000:.0f} kJ/kg, "
      f"{m['HHV_kg']:.2f} MJ/kg); LHV = {m['LHV'] / 1000:.2f} kJ/mol ({m['LHV_kg'] * 1000:.0f} kJ/kg, "
      f"{m['LHV_kg']:.2f} MJ/kg); O2/aire esteq. (100 mol) = {_O2:.1f}/{_O2 / 0.21:.1f} mol")
print(f"VEREDICTO SIMULACRO: {'OK' if _ok else 'REVISAR'}")

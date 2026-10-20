"""
energia.py — Balances de energía con reacción, combustión (HHV/LHV) y
temperatura de flama adiabática (Sesiones 16–19, 2do parcial).

Convenciones (las mismas de thermo.py y de los Excel s17/s18)
-------------------------------------------------------------
* Cp°/R = A + B·T + C·T² + D·T⁻²   (Tabla C.1 SVA), H en J/mol, T en K.
* Calor sensible de una corriente:  H_s(T) = Σ n_i · R · ICPH_i(298.15, T)
* Trayectoria SVA (fig. 4.2 / slides Sesión 18–19) para Q = ΔH:

      reactivos (T_in) ──ΔH°_R──▶ reactivos 298 K ──Σ ξ_j ΔH°298,j──▶
      productos 298 K ──ΔH°_P──▶ productos (T_out)

      ΔH°_R = −Σ_entradas H_s(T_in)     (enfriar los reactivos a 298 K)
      ΔH°_P = +H_s,salida(T_out)         (calentar los productos)
      Q = ΔH = ΔH°_R + Σ ξ_j ΔH°298,j + ΔH°_P        (W_s = 0, ΔEc ≈ ΔEp ≈ 0)

  Comprobación independiente (método de entalpías de formación):
      Q = Σ_sal n_i [ΔH°f,i + H_s,i(T_out)] − Σ_ent n_i [ΔH°f,i + H_s,i(T_in)]
  Ambos deben coincidir a precisión de máquina (si no, hay un error de
  balance de materia).
* Combustión: CxHyOzNwSs + (x + y/4 − z/2 + s) O2 → x CO2 + y/2 H2O + w/2 N2 + s SO2
      HHV = −ΔH°c con H2O(l)       LHV = −ΔH°c con H2O(g)
      HHV = LHV + n_H2O · ΔH_lv(298 K)    (ΔH_lv = 44 012 J/mol, SVA)
      Mezcla: HHV_mezcla = Σ y_i HHV_i   (inertes como N2 aportan 0)

Datos: ΔH°f y Cp de data/sva_tables.json (Tablas C.1/C.4 verificadas);
ΔH°f del agua LÍQUIDA (−285 830 J/mol) y ΔH_lv(298) no están en el JSON
(es de gases) y se toman de la Tabla C.4 de SVA (constantes de abajo).
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq

from .reactions import parse_formula, parse_reaction
from .thermo import R_J, T_REF, cargar_datos, icph

HF_H2O_LIQ = -285830.0      # J/mol, ΔH°f H2O(l) 298.15 K (Tabla C.4 SVA)
HF_H2O_GAS = -241818.0      # J/mol, ΔH°f H2O(g) 298.15 K (Tabla C.4 SVA)
DHV_H2O_298 = HF_H2O_GAS - HF_H2O_LIQ   # 44 012 J/mol (slides: 44.013 kJ/mol)

# Masas atómicas (g/mol) para pasar a MJ/kg
_MASA = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007, "S": 32.06,
         "Cl": 35.45, "Ar": 39.948}


# ----------------------------------------------------------------------------
# Conversión de unidades frecuentes en estos problemas (SVA usa (°F), lbm...)
# ----------------------------------------------------------------------------

def F_a_K(T_F: float) -> float:
    """(°F) → K."""
    return (T_F - 32.0) / 1.8 + 273.15


def C_a_K(T_C: float) -> float:
    return T_C + 273.15


LBM_A_KG = 0.45359237
BTU_A_J = 1055.056
ATM_PA = 101325.0


def masa_molar(formula: str) -> float:
    """g/mol a partir de la fórmula (C, H, O, N, S, Cl, Ar)."""
    return sum(_MASA[el] * n for el, n in parse_formula(formula).items())


# ----------------------------------------------------------------------------
# Calor sensible
# ----------------------------------------------------------------------------

def H_sensible(moles: Dict[str, float], T: float, datos: Dict[str, dict],
               T_ref: float = T_REF, avisar_Tmax: bool = True) -> float:
    """Σ n_i·R·∫_{T_ref}^{T} (Cp_i/R) dT   [J] (positivo si T > T_ref)."""
    total = 0.0
    for sp, n in moles.items():
        if n == 0.0:
            continue
        d = datos[sp]
        if avisar_Tmax and T > d.get("Tmax", float("inf")) + 1e-9:
            warnings.warn(f"H_sensible: T = {T:.1f} K > Tmax = {d['Tmax']:.0f}"
                          f" K del Cp de {sp} (extrapolación)", stacklevel=2)
        total += n * R_J * icph(T_ref, T, d["A"], d["B"], d["C"], d["D"])
    return total


def Cp_medio_H(moles: Dict[str, float], T: float, datos: Dict[str, dict],
               T_ref: float = T_REF) -> float:
    """⟨Cp⟩_H = Σ n_i ⟨Cp_i⟩_H = H_s/(T − T_ref)  [J/K] (slides Sesión 18)."""
    if abs(T - T_ref) < 1e-9:
        return sum(n * R_J * (datos[sp]["A"] + datos[sp]["B"] * T
                              + datos[sp]["C"] * T * T + datos[sp]["D"] / T ** 2)
                   for sp, n in moles.items())
    return H_sensible(moles, T, datos, T_ref, avisar_Tmax=False) / (T - T_ref)


# ----------------------------------------------------------------------------
# Calor de reacción
# ----------------------------------------------------------------------------

def dH298(nu: Dict[str, float], datos: Dict[str, dict],
          agua_liquida: bool = False) -> float:
    """ΔH°298 = Σ ν_i ΔH°f,i  [J por 'mol de reacción'].
    agua_liquida=True usa ΔH°f del H2O(l) para la especie 'H2O'."""
    total = 0.0
    for sp, v in nu.items():
        hf = (HF_H2O_LIQ if (sp == "H2O" and agua_liquida)
              else float(datos[sp]["Hf298"]))
        total += v * hf
    return total


def reaccion_combustion(formula: str, a_CO: bool = False) -> Dict[str, float]:
    """ν de la combustión de 1 mol de `formula` (completa a CO2, o a CO si
    a_CO=True). N → N2, S → SO2, H → H2O. Devuelve {especie: ν}."""
    at = parse_formula(formula)
    x, y = at.get("C", 0), at.get("H", 0)
    z, w, s = at.get("O", 0), at.get("N", 0), at.get("S", 0)
    prodC = "CO" if a_CO else "CO2"
    oC = 1 if a_CO else 2
    nO2 = (oC * x + y / 2.0 + 2 * s - z) / 2.0
    nu = {formula: -1.0, "O2": -nO2}
    if x:
        nu[prodC] = nu.get(prodC, 0.0) + x
    if y:
        nu["H2O"] = nu.get("H2O", 0.0) + y / 2.0
    if w:
        nu["N2"] = nu.get("N2", 0.0) + w / 2.0
    if s:
        nu["SO2"] = nu.get("SO2", 0.0) + s
    if formula in ("O2",):
        raise ValueError("O2 no es combustible")
    # limpia ceros (p. ej. H2 + 1/2 O2 -> H2O)
    return {k: v for k, v in nu.items() if abs(v) > 1e-12}


def calor_combustion(formula: str, datos: Optional[Dict[str, dict]] = None,
                     a_CO: bool = False) -> Dict[str, float]:
    """ΔH°c de 1 mol de combustible a 298.15 K por ley de Hess con ΔH°f.
    Devuelve dict: nu, n_H2O, dHc_l, dHc_g (J/mol), HHV, LHV (J/mol),
    M (g/mol), HHV_kg, LHV_kg (MJ/kg)."""
    nu = reaccion_combustion(formula, a_CO=a_CO)
    especies = [sp for sp in nu]
    if datos is None:
        datos = cargar_datos(especies)
    nH2O = nu.get("H2O", 0.0)
    dHl = dH298(nu, datos, agua_liquida=True)
    dHg = dH298(nu, datos, agua_liquida=False)
    M = masa_molar(formula)
    return {"nu": nu, "n_H2O": nH2O, "dHc_l": dHl, "dHc_g": dHg,
            "HHV": -dHl, "LHV": -dHg, "M": M,
            "HHV_kg": -dHl / M / 1000.0, "LHV_kg": -dHg / M / 1000.0,
            "check_HHV_LHV": (-dHl) - (-dHg) - nH2O * DHV_H2O_298}


def combustion_mezcla(y: Dict[str, float],
                      inertes: Sequence[str] = ("N2", "CO2", "H2O", "Ar"),
                      ) -> Dict[str, object]:
    """HHV/LHV de una mezcla: Σ y_i HHV_i (inertes aportan 0).
    Devuelve dict con 'por_especie' {sp: calor_combustion(sp)} y totales
    HHV, LHV (J/mol de mezcla), M (g/mol), HHV_kg, LHV_kg (MJ/kg),
    dHc_l = −HHV (J/mol de mezcla)."""
    suma = sum(y.values())
    if abs(suma - 1.0) > 1e-6:
        warnings.warn(f"combustion_mezcla: Σ y_i = {suma:.6f} ≠ 1",
                      stacklevel=2)
    por = {sp: calor_combustion(sp) for sp in y if sp not in inertes}
    HHV = sum(y[sp] * por[sp]["HHV"] for sp in por)
    LHV = sum(y[sp] * por[sp]["LHV"] for sp in por)
    M = sum(yi * masa_molar(sp) for sp, yi in y.items())
    return {"por_especie": por, "HHV": HHV, "LHV": LHV, "dHc_l": -HHV,
            "dHc_g": -LHV, "M": M, "HHV_kg": HHV / M / 1000.0,
            "LHV_kg": LHV / M / 1000.0, "suma_y": suma}


# ----------------------------------------------------------------------------
# Balances de materia auxiliares
# ----------------------------------------------------------------------------

def _atomos(moles: Dict[str, float]) -> Dict[str, float]:
    tot: Dict[str, float] = {}
    for sp, n in moles.items():
        for el, k in parse_formula(sp).items():
            tot[el] = tot.get(el, 0.0) + n * k
    return tot


def balance_atomos(entrada: Dict[str, float], salida: Dict[str, float]
                   ) -> Dict[str, float]:
    """{elemento: átomos_salida − átomos_entrada} (debe ser ~0)."""
    a_in, a_out = _atomos(entrada), _atomos(salida)
    return {el: a_out.get(el, 0.0) - a_in.get(el, 0.0)
            for el in sorted(set(a_in) | set(a_out))}


def alimentacion_por_atomos(salida: Dict[str, float],
                            especies_alimentacion: Sequence[str]
                            ) -> Dict[str, float]:
    """Resuelve las moles alimentadas de `especies_alimentacion` que cierran
    el balance atómico con la `salida` (mínimos cuadrados; avisa si el
    residual no es ~0). Ej. P4: salida CO2/CO/H2O/H2 → CH4 y H2O."""
    a_out = _atomos(salida)
    elementos = sorted(set(a_out) | {el for sp in especies_alimentacion
                                     for el in parse_formula(sp)})
    A = np.array([[parse_formula(sp).get(el, 0) for sp in especies_alimentacion]
                  for el in elementos], dtype=float)
    b = np.array([a_out.get(el, 0.0) for el in elementos])
    n, *_ = np.linalg.lstsq(A, b, rcond=None)
    res = A @ n - b
    if np.max(np.abs(res)) > 1e-9:
        warnings.warn(f"alimentacion_por_atomos: el balance no cierra "
                      f"(residual {res})", stacklevel=2)
    return {sp: float(v) for sp, v in zip(especies_alimentacion, n)}


def avances_desde_flujos(entrada: Dict[str, float], salida: Dict[str, float],
                         reacciones: Sequence[str]) -> Tuple[List[float], float]:
    """ξ_j tales que n_sal = n_ent + Σ_j ν_ij ξ_j (mínimos cuadrados).
    Devuelve (ξ, residual_max). Residual ≠ 0 → las reacciones dadas no
    explican el cambio de composición."""
    nus = [parse_reaction(r) for r in reacciones]
    especies = sorted(set(entrada) | set(salida) | {s for nu in nus for s in nu})
    A = np.array([[nu.get(sp, 0.0) for nu in nus] for sp in especies])
    b = np.array([salida.get(sp, 0.0) - entrada.get(sp, 0.0) for sp in especies])
    xi, *_ = np.linalg.lstsq(A, b, rcond=None)
    return [float(v) for v in xi], float(np.max(np.abs(A @ xi - b)))


def salida_desde_avances(entrada: Dict[str, float], reacciones: Sequence[str],
                         xi: Sequence[float]) -> Dict[str, float]:
    """n_sal,i = n_ent,i + Σ_j ν_ij ξ_j."""
    out = dict(entrada)
    for r, x in zip(reacciones, xi):
        for sp, v in parse_reaction(r).items():
            out[sp] = out.get(sp, 0.0) + v * x
    return out


# ----------------------------------------------------------------------------
# Balance de energía con reacción (Q = ΔH)
# ----------------------------------------------------------------------------

def balance_energia(entradas: Sequence[Tuple[Dict[str, float], float]],
                    salida: Dict[str, float], T_salida: float,
                    reacciones: Sequence[str], xi: Sequence[float],
                    datos: Optional[Dict[str, dict]] = None,
                    T_ref: float = T_REF) -> Dict[str, object]:
    """Q = ΔH°_R + Σ ξ_j ΔH°298,j + ΔH°_P   (trayectoria SVA por 298.15 K).

    entradas : lista de corrientes (moles, T) — p. ej. [(combustible, 298.15),
               (aire, 773.15)]; cada una se enfría a 298.15 K por separado.
    salida   : moles de la corriente de productos a T_salida.
    reacciones, xi : reacciones como texto y sus avances (mismo orden).

    Devuelve dict con dHR, dH298 (por reacción y total Σξ·ΔH), dHP, Q,
    Q_formacion (comprobación), err_materia (máx |n_sal − n_ent − νξ|),
    balance_atomos, Cp_medio_P."""
    especies = set(salida)
    for n, _T in entradas:
        especies |= set(n)
    nus = [parse_reaction(r) for r in reacciones]
    for nu in nus:
        especies |= set(nu)
    if datos is None:
        datos = cargar_datos(sorted(especies))

    ent_total: Dict[str, float] = {}
    for n, _T in entradas:
        for sp, v in n.items():
            ent_total[sp] = ent_total.get(sp, 0.0) + v

    dHR = -sum(H_sensible(n, T, datos, T_ref) for n, T in entradas) + 0.0
    dH_j = [dH298(nu, datos) for nu in nus]
    dH298_tot = sum(x * h for x, h in zip(xi, dH_j))
    dHP = H_sensible(salida, T_salida, datos, T_ref)
    Q = dHR + dH298_tot + dHP

    # comprobación por entalpías de formación (independiente de ξ)
    H_out = sum(v * datos[sp]["Hf298"] for sp, v in salida.items()) + dHP
    H_in = sum(sum(v * datos[sp]["Hf298"] for sp, v in n.items())
               + H_sensible(n, T, datos, T_ref) for n, T in entradas)
    Q_form = H_out - H_in

    pred = salida_desde_avances(ent_total, reacciones, xi)
    err = max(abs(pred.get(sp, 0.0) - salida.get(sp, 0.0))
              for sp in set(pred) | set(salida))
    return {"dHR": dHR, "dH298_j": dH_j, "dH298": dH298_tot, "dHP": dHP,
            "Q": Q, "Q_formacion": Q_form, "err_materia": err,
            "balance_atomos": balance_atomos(ent_total, salida),
            "Cp_medio_P": Cp_medio_H(salida, T_salida, datos, T_ref),
            "datos": datos}


def T_salida_para_Q(entradas: Sequence[Tuple[Dict[str, float], float]],
                    salida: Dict[str, float], reacciones: Sequence[str],
                    xi: Sequence[float], Q: float = 0.0,
                    T_lo: float = 298.15, T_hi: float = 4000.0,
                    datos: Optional[Dict[str, dict]] = None) -> float:
    """T de salida tal que Q(T) = Q dado (Q = 0 → temperatura de flama
    adiabática). brentq sobre Q(T) − Q; Q(T) es monótona creciente."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = lambda T: balance_energia(entradas, salida, T, reacciones, xi,
                                      datos)["Q"] - Q
        return brentq(f, T_lo, T_hi, xtol=1e-10)


# ----------------------------------------------------------------------------
# Agua líquida (lado del servicio: calentadores, calderas)
# ----------------------------------------------------------------------------

# h del líquido saturado [kJ/kg] vs T [°C] — Tabla A-4 de Çengel (tablas de
# vapor, misma referencia que las de SVA: líquido en el punto triple)
_HF_LIQ = [(0.01, 0.001), (5, 21.020), (10, 42.022), (15, 62.982),
           (20, 83.915), (25, 104.83), (30, 125.74), (35, 146.64),
           (40, 167.53), (45, 188.44), (50, 209.34), (55, 230.26),
           (60, 251.18), (65, 272.12), (70, 293.07), (75, 314.03),
           (80, 335.02), (85, 356.02), (90, 377.04), (95, 398.09),
           (100, 419.17)]


def h_liquido_sat(T_C: float) -> float:
    """h del agua líquida saturada [kJ/kg] por interpolación lineal
    (y = y1 + (x − x1)(y2 − y1)/(x2 − x1), slides Sesión 19)."""
    for (x1, y1), (x2, y2) in zip(_HF_LIQ, _HF_LIQ[1:]):
        if x1 <= T_C <= x2:
            return y1 + (T_C - x1) * (y2 - y1) / (x2 - x1)
    raise ValueError(f"T = {T_C} °C fuera de la tabla (0.01–100 °C)")

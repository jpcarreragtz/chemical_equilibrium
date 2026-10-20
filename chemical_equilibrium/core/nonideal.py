"""
nonideal.py — Correlaciones de NO idealidad del curso (Sesiones 13 y 15).

* phi_pitzer:  coeficiente de fugacidad de especie PURA por la correlación
  virial generalizada de Pitzer (SVA ec. 11.68 / slides Sesión 13):

      ln phi = (Pr/Tr)·(B0 + omega·B1)
      B0 = 0.083 - 0.422/Tr^1.6        B1 = 0.139 - 0.172/Tr^4.2

  Validez (criterio Fig. 3.14 de SVA): la correlación virial de dos
  términos es adecuada cuando Vr >= 2 (equivalente aproximado: Pr < 0.4
  para Tr < 1). Aquí Vr = Z·Tr/Pr con V* = R·Tc/Pc y
  Z = 1 + (B0 + omega·B1)·Pr/Tr. Si no se cumple se emite un warning y
  el resultado marca valido = False.

* wilson_gamma: coeficientes de actividad binarios por Wilson (slide 7,
  Sesión 15):

      L12 = (V2/V1)·exp(-a12/(R·T))     L21 = (V1/V2)·exp(-a21/(R·T))
      Omega = L12/(x1 + x2·L12) - L21/(x2 + x1·L21)
      ln g1 = -ln(x1 + x2·L12) + x2·Omega
      ln g2 = -ln(x2 + x1·L21) - x1·Omega

  UNIDADES: a12, a21 en cal/mol (convención DECHEMA) -> R = 1.987
  cal/(mol·K); V1, V2 en cm³/mol (solo importa su cociente).
  Límites correctos: g_i -> 1 cuando x_i -> 1.
"""

from __future__ import annotations

import warnings
from math import exp, log
from typing import Dict, Optional, Tuple

R_CAL = 1.987          # cal/(mol·K) — para los a_ij de Wilson (DECHEMA)


def phi_pitzer(T: float, P: float, Tc: float, Pc: float,
               omega: float) -> Dict[str, float]:
    """phi de especie pura por la virial generalizada.

    Devuelve dict con: phi, Tr, Pr, B0, B1, Z, Vr, valido.
    T [K], P y Pc en la MISMA unidad (bar), Tc [K], omega adimensional."""
    Tr, Pr = T / Tc, P / Pc
    B0 = 0.083 - 0.422 / Tr ** 1.6
    B1 = 0.139 - 0.172 / Tr ** 4.2
    lnphi = (Pr / Tr) * (B0 + omega * B1)
    Z = 1.0 + (B0 + omega * B1) * Pr / Tr
    Vr = Z * Tr / Pr                      # V/(R·Tc/Pc), criterio Fig. 3.14
    valido = Vr >= 2.0
    if not valido:
        warnings.warn(
            f"phi_pitzer: (Tr = {Tr:.3f}, Pr = {Pr:.3f}, Vr = {Vr:.2f} < 2) "
            f"fuera de la zona de validez de la correlación virial "
            f"(Fig. 3.14 SVA); el phi = {exp(lnphi):.4f} es una "
            f"extrapolación.", stacklevel=2)
    return {"phi": exp(lnphi), "Tr": Tr, "Pr": Pr, "B0": B0, "B1": B1,
            "Z": Z, "Vr": Vr, "valido": valido}


def wilson_lambdas(T: float, V1: float, V2: float, a12: float,
                   a21: float) -> Tuple[float, float]:
    """(L12, L21) de Wilson; a_ij en cal/mol, T en K, V_i en cm³/mol."""
    L12 = (V2 / V1) * exp(-a12 / (R_CAL * T))
    L21 = (V1 / V2) * exp(-a21 / (R_CAL * T))
    return L12, L21


def wilson_gamma(x1: float, T: float, V1: float, V2: float, a12: float,
                 a21: float) -> Tuple[float, float]:
    """(gamma1, gamma2) por Wilson para un binario en x1 (x2 = 1 - x1)."""
    if not 0.0 <= x1 <= 1.0:
        raise ValueError(f"x1 = {x1} fuera de [0, 1]")
    x2 = 1.0 - x1
    L12, L21 = wilson_lambdas(T, V1, V2, a12, a21)
    d1 = x1 + x2 * L12
    d2 = x2 + x1 * L21
    Om = L12 / d1 - L21 / d2
    ln_g1 = -log(d1) + x2 * Om
    ln_g2 = -log(d2) - x1 * Om
    return exp(ln_g1), exp(ln_g2)


def wilson_gamma_inf(T: float, V1: float, V2: float, a12: float,
                     a21: float) -> Tuple[float, float]:
    """Coeficientes a dilución infinita: ln g1_inf = 1 - ln L12 - L21;
    ln g2_inf = 1 - ln L21 - L12 (forma cerrada, para verificación)."""
    L12, L21 = wilson_lambdas(T, V1, V2, a12, a21)
    return exp(1.0 - log(L12) - L21), exp(1.0 - log(L21) - L12)


# ════════════════════════════════════════════════════════════════════════════
# 2do parcial: Margules (EQ. Físico) y regla de Lewis-Randall para K_phi
# ════════════════════════════════════════════════════════════════════════════

def margules_gamma(x1: float, A12: float, A21: Optional[float] = None
                   ) -> Tuple[float, float]:
    """(gamma1, gamma2) por Margules.

    2 parámetros:  G^E/RT = x1·x2·(A21·x1 + A12·x2)
        ln g1 = x2²·[A12 + 2(A21 − A12)·x1]
        ln g2 = x1²·[A21 + 2(A12 − A21)·x2]
    1 parámetro (A21 = None → A21 = A12 = A):  G^E/RT = A·x1·x2
        ln g1 = A·x2²      ln g2 = A·x1²
    Límites: g_i → 1 cuando x_i → 1;  ln g1_inf = A12, ln g2_inf = A21."""
    if not 0.0 <= x1 <= 1.0:
        raise ValueError(f"x1 = {x1} fuera de [0, 1]")
    if A21 is None:
        A21 = A12
    x2 = 1.0 - x1
    ln_g1 = x2 * x2 * (A12 + 2.0 * (A21 - A12) * x1)
    ln_g2 = x1 * x1 * (A21 + 2.0 * (A12 - A21) * x2)
    return exp(ln_g1), exp(ln_g2)


def margules_GE_RT(x1: float, A12: float, A21: Optional[float] = None) -> float:
    """G^E/RT = x1·x2·(A21·x1 + A12·x2)."""
    if A21 is None:
        A21 = A12
    x2 = 1.0 - x1
    return x1 * x2 * (A21 * x1 + A12 * x2)


def cargar_criticas(especies, ruta=None) -> Dict[str, Tuple[float, float, float]]:
    """{especie: (Tc [K], Pc [bar], omega)} desde data/criticas_b1.json."""
    import json
    from . import DATA
    ruta = ruta or str(DATA / "criticas_b1.json")
    with open(ruta, encoding="utf-8") as f:
        tab = json.load(f)
    faltan = [sp for sp in especies if sp not in tab]
    if faltan:
        raise KeyError(f"Sin Tc/Pc/omega en criticas_b1.json: {faltan} "
                       f"(agrégalas de la Tabla B.1 de SVA)")
    return {sp: (tab[sp]["Tc"], tab[sp]["Pc"], tab[sp]["omega"])
            for sp in especies}


def K_phi_lewis(nu: Dict[str, float], T: float, P: float,
                criticas: Dict[str, Tuple[float, float, float]],
                avisar: bool = True) -> Tuple[float, Dict[str, dict]]:
    """K_phi = Π phi_i^nu_i con phi_i de especie PURA a (T, P) — regla de
    Lewis-Randall = "solución ideal de gases" (SVA sec. 13.7 / Sesión 13):

        K = K_phi · K_y · (P/P°)^δ      →     K_y = K·(P/P°)^(−δ) / K_phi
        (forma rigurosa, SVA ec. 13.28)

    CÓMO USARLO: esta función solo devuelve Π phi_i^nu_i. De qué lado de la
    ecuación entra (convención "maestra" o "rigurosa") lo decide ÚNICAMENTE
    core/convenciones.py (factor_Kphi, K_efectiva); no lo decidas aquí ni
    en los scripts.

    Devuelve (K_phi, {especie: dict de phi_pitzer})."""
    import warnings as _w
    detalle = {}
    with _w.catch_warnings():
        if not avisar:
            _w.simplefilter("ignore")
        for sp in nu:
            detalle[sp] = phi_pitzer(T, P, *criticas[sp])
    Kphi = 1.0
    for sp, v in nu.items():
        Kphi *= detalle[sp]["phi"] ** v
    return Kphi, detalle


def phi_virial_vdw(T: float, P: float, Tc: float, Pc: float) -> Dict[str, float]:
    """phi de especie pura por la virial TRUNCADA con B de van der Waals
    (Sesión 13, ej. 9.25):  Z = 1 + B·P/(R·T),  B = b − a/(R·T),
        a = 27 R² Tc²/(64 Pc),  b = R·Tc/(8 Pc),   ln phi = B·P/(R·T)
    T, Tc [K]; P, Pc [bar]. Devuelve dict con phi, B [m³/mol], a, b, Z."""
    Rb = 8.314e-5                        # bar·m³/(mol·K)
    a = 27.0 * Rb * Rb * Tc * Tc / (64.0 * Pc)
    b = Rb * Tc / (8.0 * Pc)
    B = b - a / (Rb * T)
    lnphi = B * P / (Rb * T)
    return {"phi": exp(lnphi), "B": B, "a": a, "b": b, "Z": 1.0 + lnphi}

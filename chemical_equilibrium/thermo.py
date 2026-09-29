"""
thermo.py — Cp(T), ΔH°(T), ΔG°(T) y K(T) con las convenciones de
Smith, Van Ness & Abbott (SVA).

Convenciones
------------
* Cp°/R = A + B·T + C·T² + D·T⁻²          (Tabla C.1 de SVA, T en K)
* ΔH°(T) = ΔH°298 + R·ICPH(298.15, T)     con ICPH = ∫ (ΔCp°/R) dT
* ΔG°(T)/RT = (ΔG°298 − ΔH°298)/(R·T0) + ΔH°298/(R·T)
              + (1/T)·ICPH − ICPS         (ec. 13.18 de SVA)
  con ICPS = ∫ (ΔCp°/R) dT/T  y  T0 = 298.15 K.
* K(T) = exp(−ΔG°(T)/RT)
* Van't Hoff con ΔH° constante:  ln(K/K0) = −(ΔH°/R)(1/T − 1/T0)

Datos por especie (dict `datos`):
    {"NH3": {"Hf298": -46110.0, "Gf298": -16450.0,
             "A": 3.578, "B": 3.020e-3, "C": 0.0, "D": -0.186e5}, ...}
Hf298 y Gf298 en J/mol (Tabla C.4 de SVA); A, B, C, D adimensionales
(Cp/R, Tabla C.1). Si faltan A..D se toman como 0 (equivale a ΔCp = 0,
es decir Van't Hoff con ΔH constante).
"""

from __future__ import annotations

import json
import os
from math import exp, log
from typing import Dict, Optional, Sequence

R_J = 8.314          # J/(mol·K)
T_REF = 298.15       # K

# Tablas C.1/C.4 capturadas en JSON (ver data/raw/LEEME.txt sobre su estatus
# de verificación)
RUTA_TABLAS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "data", "sva_tables.json")


# ----------------------------------------------------------------------------
# Cp y sus integrales (formas ICPH / ICPS de SVA)
# ----------------------------------------------------------------------------

def cp_R(T: float, A: float, B: float = 0.0, C: float = 0.0,
         D: float = 0.0) -> float:
    """Cp°/R = A + B·T + C·T² + D·T⁻² (adimensional)."""
    return A + B * T + C * T * T + D / (T * T)


def icph(T0: float, T: float, A: float, B: float = 0.0, C: float = 0.0,
         D: float = 0.0) -> float:
    """ICPH = ∫_{T0}^{T} (Cp/R) dT  [K]."""
    return (A * (T - T0) + B / 2.0 * (T * T - T0 * T0)
            + C / 3.0 * (T ** 3 - T0 ** 3) - D * (1.0 / T - 1.0 / T0))


def icps(T0: float, T: float, A: float, B: float = 0.0, C: float = 0.0,
         D: float = 0.0) -> float:
    """ICPS = ∫_{T0}^{T} (Cp/R) dT/T  [adimensional]."""
    return (A * log(T / T0) + B * (T - T0)
            + C / 2.0 * (T * T - T0 * T0)
            - D / 2.0 * (1.0 / (T * T) - 1.0 / (T0 * T0)))


# ----------------------------------------------------------------------------
# Propiedades de reacción
# ----------------------------------------------------------------------------

_CAMPOS_CP = ("A", "B", "C", "D")


def delta_property(nu: Dict[str, float], datos: Dict[str, dict],
                   campo: str, obligatorio: bool = False) -> float:
    """Δ(propiedad) = Σ_i ν_i · propiedad_i sobre las especies de la reacción."""
    total = 0.0
    for sp, coef in nu.items():
        if sp not in datos:
            raise KeyError(f"Faltan datos termodinámicos de la especie '{sp}'")
        if campo not in datos[sp]:
            if obligatorio:
                raise KeyError(f"Falta '{campo}' para la especie '{sp}'")
            continue
        total += coef * float(datos[sp][campo])
    return total


def reaction_thermo(nu: Dict[str, float], datos: Dict[str, dict],
                    T: float) -> Dict[str, float]:
    """
    Para UNA reacción (nu = {especie: ν_i}) y una T [K] devuelve:
        dH298, dG298   [J/mol]  (a 298.15 K)
        dCp_R          ΔCp°/R evaluado en T
        dH             ΔH°(T)  [J/mol]
        dG             ΔG°(T)  [J/mol]
        K              exp(−ΔG°(T)/RT)
    """
    dH0 = delta_property(nu, datos, "Hf298", obligatorio=True)
    dG0 = delta_property(nu, datos, "Gf298", obligatorio=True)
    dA, dB, dC, dD = (delta_property(nu, datos, c) for c in _CAMPOS_CP)

    idh = icph(T_REF, T, dA, dB, dC, dD)
    ids = icps(T_REF, T, dA, dB, dC, dD)

    dH = dH0 + R_J * idh
    dG_RT = ((dG0 - dH0) / (R_J * T_REF) + dH0 / (R_J * T)
             + idh / T - ids)                      # ec. 13.18 SVA
    return {
        "dH298": dH0, "dG298": dG0,
        "dCp_R": cp_R(T, dA, dB, dC, dD),
        "dH": dH, "dG": dG_RT * R_J * T,
        "K": exp(-dG_RT),
    }


def K_de_datos(nu: Dict[str, float], datos: Dict[str, dict],
               T: float) -> float:
    """K(T) de una reacción a partir de datos de formación + Cp (ec. 13.18)."""
    return reaction_thermo(nu, datos, T)["K"]


def K_vant_hoff(T: float, K0: float, dH: float, T0: float = T_REF) -> float:
    """Van't Hoff integrado con ΔH° constante:
    K(T) = K0 · exp(−(ΔH°/R)(1/T − 1/T0))."""
    return K0 * exp(-(dH / R_J) * (1.0 / T - 1.0 / T0))


def tabla_termo(nu_reacciones: Sequence[Dict[str, float]],
                datos: Dict[str, dict],
                temperaturas: Sequence[float]) -> Dict[float, list]:
    """{T: [reaction_thermo(R1), reaction_thermo(R2), ...]} para un barrido."""
    return {T: [reaction_thermo(nu, datos, T) for nu in nu_reacciones]
            for T in temperaturas}


# ----------------------------------------------------------------------------
# Carga de especies desde data/sva_tables.json
# ----------------------------------------------------------------------------

def cargar_datos(especies: Sequence[str],
                 ruta: Optional[str] = None) -> Dict[str, dict]:
    """
    Carga las especies pedidas desde data/sva_tables.json y devuelve el dict
    {especie: {Hf298, Gf298, A, B, C, D, Tmax, nombre, verificado}} que
    consumen reaction_thermo / run_exercise.

    Lanza KeyError con la lista de faltantes (y de disponibles) si alguna
    especie no está en el JSON.
    """
    ruta = ruta or RUTA_TABLAS
    if not os.path.exists(ruta):
        raise FileNotFoundError(f"No existe {ruta} — captura ahí las Tablas "
                                f"C.1/C.4 (ver data/raw/LEEME.txt)")
    with open(ruta, encoding="utf-8") as f:
        tablas = json.load(f)
    disponibles = sorted(k for k in tablas if not k.startswith("_"))
    faltan = [sp for sp in especies if sp not in tablas]
    if faltan:
        raise KeyError(f"Especies sin datos en {os.path.basename(ruta)}: "
                       f"{', '.join(faltan)}. Disponibles: "
                       f"{', '.join(disponibles)}")
    datos: Dict[str, dict] = {}
    for sp in especies:
        e = tablas[sp]
        cp = e.get("Cp", {})
        datos[sp] = {
            "Hf298": float(e["Hf298"]), "Gf298": float(e["Gf298"]),
            "A": float(cp.get("A", 0.0)), "B": float(cp.get("B", 0.0)),
            "C": float(cp.get("C", 0.0)), "D": float(cp.get("D", 0.0)),
            "Tmax": float(cp.get("Tmax", float("inf"))),
            "nombre": e.get("nombre", sp),
            "verificado": e.get("verificado", "?"),
        }
    return datos

"""
gibbs.py — Equilibrio por minimización de la energía de Gibbs (Sesión 25).

No necesita reacciones ni valor inicial: solo las especies, sus ΔG°f y la
alimentación. Resuelve, para gas en una fase (SVA ec. 13.41 en adelante):

    min  nG/RT = Σ n_i [ ΔG°f,i/RT + ln( y_i φ_i P/P0 ) ]
    sujeto a  A n = b,   n_i > 0

con A = átomos por especie y b = A n0. Sus condiciones de Lagrange son

    ΔG°f,i/RT + ln( y_i φ_i P/P0 ) + Σ_e A_ei λ_e/RT = 0    (una por especie)

Algoritmo: Newton sobre esas condiciones en ln n_i (RAND de White, Johnson y
Dantzig; forma de Gordon y McBride), con paso amortiguado. Trabajar en ln n_i
evita ln(0), y converger TODAS las ln n_i es lo que deja bien las trazas: un
minimizador general se detiene antes porque casi no pesan en G.

Alcance: una fase gaseosa ideal con φ_i constante por especie, más
  * restricciones lineales extra (`restricciones`, p. ej. ASF en
    Fischer–Tropsch): se agregan como filas a A n = b, cada una con su
    multiplicador;
  * fases condensadas puras con actividad 1 (`condensadas`, p. ej. "C(s)",
    "Si(s)", "S(ref)", "H2O(l)"), con la prueba de fase de la Sesión 28:
    se resuelve el gas, Δ_k = g_k/RT + Σ_e A_ek λ_e/RT; si Δ_k < 0 la fase
    entra con g_k/RT + Σ A λ = 0 y n_k ≥ 0, y si su n sale negativo se quita.
Con φ_i dependiente de y, o con soluciones líquidas, nG/RT cambia y este
módulo no aplica.

Datos (`dGf_janaf`, `dGf_sva`, `cargar_dgf`):
    * Primero data/janaf_dgf.json, capturado de https://janaf.nist.gov/.
    * ec. 13.18 de SVA solo para especies declaradas `no_en_janaf`.

Uso:
    res = minimizar_gibbs(especies, dGf, n0, T=1000, P=1)
    imprimir_resultado(res)
    validar_gibbs(res)
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from . import DATA
from .independencia import (elegir_componentes, escalar_columnas,
                            matriz_atomica, nu_completa, nu_por_especie,
                            rango, reaccion_desde_columna,
                            reducir_por_componentes, transpuesta)
from .reactions import parse_formula
from .selfcheck import Report, _box
from .thermo import R_J, cargar_datos, reaction_thermo

RUTA_JANAF = DATA / "janaf_dgf.json"

# Estado de referencia por elemento: (especie, átomos del elemento por molécula)
_REFERENCIA = {"C": ("C(grafito)", 1), "H": ("H2", 2), "O": ("O2", 2),
               "N": ("N2", 2), "Cl": ("Cl2", 2)}
# Tabla C.2 de SVA, carbono (grafito); no está en sva_tables.json (solo gases)
_GRAFITO = {"Hf298": 0.0, "Gf298": 0.0, "A": 1.771, "B": 0.771e-3,
            "C": 0.0, "D": -0.867e5}
# Sufijo de fase de una especie condensada: "Si(s)", "H2O(l)", "S(ref)"
_FASE = re.compile(r"\((?:s|l|cr|liq|ref)\)$")


def _composicion(especie: str) -> Dict[str, int]:
    """Átomos de una especie; ignora el sufijo de fase ("Si(s)" -> {Si: 1})."""
    return parse_formula(_FASE.sub("", especie))


# ----------------------------------------------------------------------------
# Datos de ΔG°f
# ----------------------------------------------------------------------------

def dGf_sva(especie: str, T: float) -> float:
    """ΔG°f(T) [J/mol] con la ec. 13.18 de SVA (reaction_thermo) para la
    formación de `especie` desde sus elementos."""
    nu: Dict[str, float] = {especie: 1.0}
    for el, c in parse_formula(especie).items():
        if el not in _REFERENCIA:
            raise KeyError(f"Sin estado de referencia para el elemento '{el}'")
        ref, atomos = _REFERENCIA[el]
        nu[ref] = nu.get(ref, 0.0) - c / atomos
    datos = cargar_datos([sp for sp in nu if sp != "C(grafito)"])
    datos["C(grafito)"] = dict(_GRAFITO)
    return reaction_thermo(nu, datos, T)["dG"]


def dGf_janaf(especie: str, T: float) -> Tuple[float, str]:
    """(ΔG°f [J/mol], tabla) capturado de JANAF. KeyError si no está a esa T."""
    with open(RUTA_JANAF, encoding="utf-8") as f:
        tabla = json.load(f)
    clave = f"{T:g}"
    if especie not in tabla or clave not in tabla[especie]:
        raise KeyError(f"Falta ΔG°f de {especie} a {T:g} K en "
                       f"data/janaf_dgf.json: búscalo en "
                       f"https://janaf.nist.gov/ y agrégalo")
    return float(tabla[especie][clave]), tabla[especie]["tabla"]


def cargar_dgf(especies: Sequence[str], T: float
               ) -> Tuple[Dict[str, float], Dict[str, str]]:
    """({especie: ΔG°f [J/mol]}, {especie: fuente}). Usa JANAF; solo recurre
    a SVA para las especies declaradas `no_en_janaf` en el JSON."""
    with open(RUTA_JANAF, encoding="utf-8") as f:
        no_en_janaf = json.load(f)["_meta"]["no_en_janaf"]
    dGf, fuentes = {}, {}
    for sp in especies:
        if sp in no_en_janaf:
            dGf[sp] = dGf_sva(sp, T)
            fuentes[sp] = "SVA ec. 13.18 (no está en JANAF)"
        else:
            dGf[sp], tabla = dGf_janaf(sp, T)
            fuentes[sp] = f"JANAF {tabla}"
    return dGf, fuentes


# ----------------------------------------------------------------------------
# Solver
# ----------------------------------------------------------------------------

@dataclass
class ResultadoGibbs:
    especies: List[str]
    elementos: List[str]
    n: Dict[str, float]            # moles en el equilibrio
    y: Dict[str, float]            # fracciones mol
    n_total: float
    lam_RT: Dict[str, float]       # λ_e/RT de la convención de la Sesión 25
    G_RT: float                    # (nG)_ref/RT en el equilibrio (ec. 11)
    G_RT_alim: float               # (nG)_ref/RT de la alimentación sin reaccionar
    iteraciones: int
    T: float
    P: float
    P0: float
    dGf: Dict[str, float]
    n0: Dict[str, float]
    phi: Dict[str, float]
    excluidas: List[str] = field(default_factory=list)
    # restricciones lineales extra: [({especie: coef}, rhs)] y sus λ/RT
    restricciones: List[Tuple[Dict[str, float], float]] = field(default_factory=list)
    lam_restr: List[float] = field(default_factory=list)
    # fases condensadas: candidatas, presentes {especie: n}, Δ_k final y Δ_k
    # del gas solo (prueba de fase), y la bitácora de altas y bajas
    condensadas: List[str] = field(default_factory=list)
    n_cond: Dict[str, float] = field(default_factory=dict)
    delta_fase: Dict[str, float] = field(default_factory=dict)
    delta_fase_gas: Dict[str, float] = field(default_factory=dict)
    bitacora_fases: List[str] = field(default_factory=list)


def _G_RT(n: np.ndarray, g: np.ndarray) -> float:
    """nG/RT; los términos con n_i = 0 valen 0 (límite de n ln n)."""
    N = n.sum()
    m = n > 0
    return float(np.sum(n[m] * (g[m] + np.log(n[m] / N))))


def _rand(g: np.ndarray, A: np.ndarray, b: np.ndarray, lnn: np.ndarray,
          tol: float, max_iter: int, gc: Optional[np.ndarray] = None,
          Ac: Optional[np.ndarray] = None, nc: Optional[np.ndarray] = None
          ) -> Tuple[np.ndarray, np.ndarray, int, np.ndarray]:
    """Newton de Gordon-McBride. A = filas (elementos y restricciones) ×
    especies gaseosas. Con fases condensadas puras (g_c, columnas Ac, n_c)
    se agrega una ecuación lineal por fase, Σ_e A_ec π_e = g_c, y su Δn_c
    entra en los balances. Devuelve (ln n, λ/RT, iteraciones, n_c).
    `lstsq` tolera A sin rango completo."""
    m = A.shape[0]
    k = 0 if gc is None else len(gc)
    if k:
        nc = np.asarray(nc, dtype=float)
    else:
        gc, Ac, nc = np.zeros(0), np.zeros((m, 0)), np.zeros(0)
    lnN = float(np.log(np.exp(lnn).sum()))
    tol_b = tol * max(float(np.abs(b).max()), 1.0)
    for it in range(1, max_iter + 1):
        n = np.exp(lnn)
        N = np.exp(lnN)
        mu = g + lnn - lnN
        An = A @ n
        M = np.zeros((m + 1 + k, m + 1 + k))
        M[:m, :m] = (A * n) @ A.T
        M[:m, m] = An
        M[m, :m] = An
        M[m, m] = n.sum() - N
        r = np.zeros(m + 1 + k)
        r[:m] = b - An + (A * n) @ mu
        r[m] = N - n.sum() + n @ mu
        if k:
            M[:m, m + 1:] = Ac
            M[m + 1:, :m] = Ac.T
            r[:m] -= Ac @ nc
            r[m + 1:] = gc
        sol = np.linalg.lstsq(M, r, rcond=None)[0]
        pi, dlnN, dnc = sol[:m], sol[m], sol[m + 1:]
        dlnn = dlnN + A.T @ pi - mu
        if (max(np.abs(dlnn).max(), abs(dlnN)) < tol
                and np.abs(An + Ac @ nc - b).max() < tol_b
                and (not k or np.abs(dnc).max() < tol_b)):
            return lnn, -pi, it, nc
        den = max(5.0 * abs(dlnN), np.abs(dlnn).max())
        paso = 1.0 if den <= 2.0 else 2.0 / den
        lnn = lnn + paso * dlnn
        lnN = lnN + paso * dlnN
        nc = nc + paso * dnc
    raise RuntimeError(f"minimizar_gibbs no convergió en {max_iter} "
                       f"iteraciones")


def minimizar_gibbs(especies: Sequence[str],
                    dGf: Mapping[str, float],
                    n0: Mapping[str, float],
                    T: float,
                    P: float = 1.0,
                    P0: float = 1.0,
                    phi: Optional[Mapping[str, float]] = None,
                    arranque: Optional[Mapping[str, float]] = None,
                    tol: float = 1e-10,
                    max_iter: int = 20000,
                    restricciones: Optional[Sequence[Tuple[Mapping[str, float], float]]] = None,
                    condensadas: Optional[Sequence[str]] = None,
                    fases_iniciales: Optional[Sequence[str]] = None
                    ) -> ResultadoGibbs:
    """
    Composición de equilibrio de un gas en una fase a T y P, opcionalmente
    con restricciones lineales extra y fases condensadas puras.

    especies : todas las especies gaseosas candidatas.
    dGf      : {especie: ΔG°f(T) [J/mol]}, también de las condensadas.
    n0       : {especie gaseosa: mol alimentados}; las ausentes valen 0.
    T, P, P0 : K; P y P0 en las mismas unidades (solo entra P/P0).
    phi      : {especie: φ_i constante}; las ausentes valen 1.
    arranque : {especie: mol} opcional y positivo. Sin él se usa un reparto
               uniforme, así que no hace falta adivinar nada.
    restricciones : [({especie gaseosa: coef}, rhs), ...]; cada una agrega
               Σ coef·n = rhs a A n = b (p. ej. ASF: n_Ck − α^(k−1) n_C1 = 0).
    condensadas : especies puras con actividad 1 ("C(s)", "Si(s)", "S(ref)",
               "H2O(l)"); entran solo si la prueba de fase lo pide.
    fases_iniciales : condensadas que se incluyen desde el principio (se
               quitan si su n sale negativo).

    Una especie con un elemento que no está en la alimentación no puede
    existir: se excluye (n = 0) y queda en `excluidas`.
    """
    especies = list(especies)
    if len(set(especies)) != len(especies):
        raise ValueError("Hay especies repetidas")
    if T <= 0 or P <= 0 or P0 <= 0:
        raise ValueError("T, P y P0 deben ser positivas")
    faltan = [sp for sp in especies if sp not in dGf]
    if faltan:
        raise KeyError("Falta ΔG°f de: " + ", ".join(faltan))
    ajenas = [sp for sp in n0 if sp not in especies]
    if ajenas:
        raise ValueError("Alimentación fuera de `especies`: " + ", ".join(ajenas))
    if any(v < 0 for v in n0.values()) or sum(n0.values()) <= 0:
        raise ValueError("La alimentación debe ser no negativa y no nula")
    phi = dict(phi or {})
    if any(v <= 0 for v in phi.values()):
        raise ValueError("Los φ_i deben ser positivos")
    condensadas = list(condensadas or [])
    fases_iniciales = list(fases_iniciales or [])
    malas = [k for k in condensadas if k in especies or k not in dGf]
    if malas or len(set(condensadas)) != len(condensadas):
        raise ValueError("Condensadas repetidas, sin ΔG°f o también en "
                         "`especies`: " + ", ".join(malas))
    if any(k not in condensadas for k in fases_iniciales):
        raise ValueError("`fases_iniciales` debe estar dentro de `condensadas`")
    restricciones = [({sp: float(c) for sp, c in coef.items()}, float(rhs))
                     for coef, rhs in (restricciones or [])]
    ajenas = sorted({sp for coef, _ in restricciones for sp in coef
                     if sp not in especies})
    if ajenas:
        raise ValueError("Restricciones con especies fuera de `especies`: "
                         + ", ".join(ajenas))

    A_ex, elementos = matriz_atomica(especies)
    A_tot = np.array(transpuesta(A_ex), dtype=float)      # elementos × especies
    n0_tot = np.array([float(n0.get(sp, 0.0)) for sp in especies])
    b_tot = A_tot @ n0_tot

    # Una especie con un elemento ausente de la alimentación no puede formarse
    posible = ~np.any((A_tot > 0) & (b_tot[:, None] == 0), axis=0)
    activos = b_tot > 0
    A, b = A_tot[activos][:, posible], b_tot[activos]
    elem = [el for el, a in zip(elementos, activos) if a]
    sp_act = [sp for sp, p in zip(especies, posible) if p]
    excluidas = [sp for sp, p in zip(especies, posible) if not p]

    g_tot = (np.array([dGf[sp] for sp in especies]) / (R_J * T)
             + np.log(P / P0)
             + np.log([phi.get(sp, 1.0) for sp in especies]))
    g = g_tot[posible]

    if arranque is None:
        lnn = np.full(len(sp_act), np.log(b.sum() / (10.0 * len(sp_act))))
    else:
        x = np.array([float(arranque[sp]) for sp in sp_act])
        if np.any(x <= 0):
            raise ValueError("El arranque debe ser positivo")
        lnn = np.log(x)

    # Restricciones extra: filas adicionales de A n = b
    nr = len(restricciones)
    if nr:
        A = np.vstack([A, [[coef.get(sp, 0.0) for sp in sp_act]
                           for coef, _ in restricciones]])
        b = np.concatenate([b, [rhs for _, rhs in restricciones]])

    # Fases condensadas: solo elementos que ya están en alguna especie gaseosa
    comp_c = {k: _composicion(k) for k in condensadas}
    sin_gas = [k for k in condensadas if any(el not in elementos for el in comp_c[k])]
    if sin_gas:
        raise ValueError("Condensadas con elementos que ninguna especie "
                         "gaseosa tiene: " + ", ".join(sin_gas))
    cond = [k for k in condensadas
            if all(b_tot[elementos.index(el)] > 0 for el in comp_c[k])]
    excluidas += [k for k in condensadas if k not in cond]
    g_c = {k: dGf[k] / (R_J * T) for k in cond}
    Ac_col = {k: np.array([comp_c[k].get(el, 0) for el in elem] + [0.0] * nr)
              for k in cond}

    def resolver(presentes, lnn0, nc0):
        if not presentes:
            return _rand(g, A, b, lnn0, tol, max_iter)
        return _rand(g, A, b, lnn0, tol, max_iter,
                     np.array([g_c[k] for k in presentes]),
                     np.column_stack([Ac_col[k] for k in presentes]), nc0)

    def deltas(lam):
        return {k: g_c[k] + sum(c * lam[elem.index(el)]
                                for el, c in comp_c[k].items()) for k in cond}

    lnn, lam, it, nc = resolver([], lnn, None)
    delta_gas = deltas(lam)
    presentes, bitacora = [k for k in fases_iniciales if k in cond], []
    if presentes:
        bitacora.append("se incluyen desde el principio: " + ", ".join(presentes))
        lnn, lam, it2, nc = resolver(presentes, lnn, np.zeros(len(presentes)))
        it += it2
    tol_n = tol * max(float(np.abs(b).max()), 1.0)
    for _ in range(4 * len(cond) + 4):
        if presentes and nc.min() < -tol_n:                 # fase con n < 0: fuera
            j = int(np.argmin(nc))
            bitacora.append(f"{presentes[j]}: n = {nc[j]:.4g} < 0, se quita")
            presentes.pop(j)
            nc = np.delete(nc, j)
        else:
            D = deltas(lam)
            ausentes = {k: D[k] for k in cond if k not in presentes}
            if not ausentes or min(ausentes.values()) >= -1e-9:
                break
            k = min(ausentes, key=ausentes.get)             # Δ_k < 0: dentro
            bitacora.append(f"{k}: Δ = {ausentes[k]:.6g} < 0, se incluye")
            presentes.append(k)
            nc = np.append(nc, 0.0)
        lnn, lam, it2, nc = resolver(presentes, lnn, nc)
        it += it2
    else:
        raise RuntimeError("La prueba de fase no se estabilizó: "
                           + "; ".join(bitacora))

    n_act = np.exp(lnn)
    N = float(n_act.sum())
    n = {sp: 0.0 for sp in especies}
    n.update(zip(sp_act, map(float, n_act)))
    m_el = len(elem)
    G_cond = sum(float(v) * g_c[k] for k, v in zip(presentes, nc))
    return ResultadoGibbs(
        especies=especies, elementos=elem, n=n,
        y={sp: v / N for sp, v in n.items()}, n_total=N,
        lam_RT=dict(zip(elem, map(float, lam[:m_el]))),
        G_RT=_G_RT(n_act, g) + G_cond, G_RT_alim=_G_RT(n0_tot, g_tot),
        iteraciones=it, T=T, P=P, P0=P0, dGf=dict(dGf),
        n0={sp: float(n0.get(sp, 0.0)) for sp in especies},
        phi={sp: float(phi.get(sp, 1.0)) for sp in especies},
        excluidas=excluidas,
        restricciones=restricciones, lam_restr=[float(v) for v in lam[m_el:]],
        condensadas=condensadas, n_cond=dict(zip(presentes, map(float, nc))),
        delta_fase={k: float(v) for k, v in deltas(lam).items()},
        delta_fase_gas={k: float(v) for k, v in delta_gas.items()},
        bitacora_fases=bitacora)


def imprimir_resultado(res: ResultadoGibbs) -> None:
    print(f"\n  T = {res.T:g} K, P = {res.P:g}; {res.iteraciones} iteraciones")
    print(f"  {'especie':<9}{'n_i (mol)':>14}{'y_i':>14}")
    for sp in sorted(res.especies, key=lambda s: -res.y[s]):
        print(f"  {sp:<9}{res.n[sp]:>14.5e}{res.y[sp]:>14.5e}")
    print(f"  {'total':<9}{res.n_total:>14.5e}")
    print("  λ/RT: " + ", ".join(f"{k} = {v:.4f}" for k, v in res.lam_RT.items()))
    if res.restricciones:
        print("  λ/RT de las restricciones extra: "
              + ", ".join(f"r{i + 1} = {v:.4f}" for i, v in enumerate(res.lam_restr)))
    if res.condensadas:
        print(f"  {'fase condensada':<16}{'Δ_k solo gas':>14}{'Δ_k final':>12}{'n (mol)':>14}")
        for k in res.condensadas:
            if k in res.excluidas:
                print(f"  {k:<16}{'excluida (elemento ausente)':>40}")
                continue
            print(f"  {k:<16}{res.delta_fase_gas[k]:>14.5f}{res.delta_fase[k]:>12.5f}"
                  f"{res.n_cond.get(k, 0.0):>14.5e}")
        for linea in res.bitacora_fases:
            print(f"    prueba de fase: {linea}")
    print(f"  (nG)ref/RT: {res.G_RT_alim:.4f} (alimentación) -> "
          f"{res.G_RT:.4f} (equilibrio)")


# ----------------------------------------------------------------------------
# Reacciones independientes y validación
# ----------------------------------------------------------------------------

def reacciones_independientes(res: ResultadoGibbs
                              ) -> Tuple[List[str], np.ndarray, List[str]]:
    """Reacciones independientes con el criterio de CLAUDE.md (componentes:
    la alimentación primero). Devuelve (especies, ν por reacción, textos).
    Si no hay reacciones, la ν es una matriz vacía. Las fases condensadas
    presentes entran como especies (con actividad 1 en las K)."""
    especies = ([sp for sp in res.especies if sp not in res.excluidas]
                + list(res.n_cond))
    A, _ = matriz_atomica([_FASE.sub("", sp) for sp in especies])
    alim = [sp for sp in especies if res.n0.get(sp, 0.0) > 0]
    comp, _ = elegir_componentes(especies, A, alim)
    red = reducir_por_componentes(transpuesta(A), especies, comp)
    if not red["no_componentes"]:
        return especies, np.zeros((0, len(especies))), []
    nu_frac, filas = nu_completa(red)
    nu, _ = escalar_columnas(nu_frac)
    textos = [reaccion_desde_columna(col, filas, especies)
              for col in transpuesta(nu)]
    nu_esp = np.array(nu_por_especie(nu, filas, especies), dtype=float)
    return especies, nu_esp, textos


def validar_gibbs(res: ResultadoGibbs, n_arranques: int = 20,
                  semilla: int = 0, tol: float = 1e-7) -> Report:
    """
    Autovalidación, sin clave de respuestas:
      1. balance de átomos in/out
      2. n_i > 0 y Σ y_i = 1
      3. condiciones de Lagrange, especie por especie
      4. K de cada reacción independiente: ln K(ΔG°) vs Σ ν ln(y φ P/P0)
      5. robustez: mismas ln y (trazas incluidas) desde arranques al azar
      6. nG/RT baja respecto a la alimentación
      7. carbono: a_C < 1, o el resultado solo-gas no es estable
    Solo si aplican:
      8. restricciones extra: Σ coef·n = rhs (con ellas el punto 4 es n/a,
         porque el equilibrio está restringido y las K no se cumplen)
      9. fases condensadas: Δ_k ≥ 0 en las ausentes y n_k ≥ 0 en las presentes
    """
    _box(f"Validación Gibbs ({res.T:g} K, P = {res.P:g})")
    rep = Report()
    especies = [sp for sp in res.especies if sp not in res.excluidas]
    A_ex, elementos = matriz_atomica(res.especies)
    A_tot = np.array(transpuesta(A_ex), dtype=float)
    n = np.array([res.n[sp] for sp in res.especies])
    n0 = np.array([res.n0[sp] for sp in res.especies])
    b = A_tot @ n0
    sale_cond = np.zeros(len(elementos))
    for k, nk in res.n_cond.items():
        for el, c in _composicion(k).items():
            sale_cond[elementos.index(el)] += c * nk

    # 1
    rel = float(np.max(np.abs(A_tot @ n + sale_cond - b) / np.maximum(b, 1.0)))
    rep.add("pass" if rel < tol else "fail", 1,
            f"balance de átomos: máx |entra − sale| relativo = {rel:.1e}")
    # 2
    y = np.array([res.y[sp] for sp in especies])
    ok = bool(np.all(y > 0)) and abs(y.sum() - 1.0) < tol
    rep.add("pass" if ok else "fail", 2,
            f"n_i > 0 y Σy_i = {y.sum():.12f}")
    # 3 (con los multiplicadores de las restricciones y las fases presentes)
    g = {sp: res.dGf[sp] / (R_J * res.T) + np.log(res.phi[sp])
         + np.log(res.P / res.P0) for sp in especies}
    comp = {sp: parse_formula(sp) for sp in especies}
    resid = [g[sp] + np.log(res.y[sp])
             + sum(c * res.lam_RT[el] for el, c in comp[sp].items())
             + sum(coef.get(sp, 0.0) * lr
                   for (coef, _), lr in zip(res.restricciones, res.lam_restr))
             for sp in especies]
    resid += [res.dGf[k] / (R_J * res.T)
              + sum(c * res.lam_RT[el] for el, c in _composicion(k).items())
              for k in res.n_cond]
    r3 = float(np.max(np.abs(resid)))
    rep.add("pass" if r3 < tol else "fail", 3,
            f"condiciones de Lagrange: máx |residuo| = {r3:.1e}"
            + (f" (con {len(res.restricciones)} multiplicadores de restricciones)"
               if res.restricciones else "")
            + (f" (con {len(res.n_cond)} fase(s) condensada(s))"
               if res.n_cond else ""))
    # 4
    sp4, nu, textos = reacciones_independientes(res)
    if res.restricciones:
        rep.add("na", 4, "equilibrio restringido por "
                f"{len(res.restricciones)} restricciones extra: las K no se cumplen")
    elif len(textos):
        ln_a = np.array([0.0 if sp in res.n_cond else
                         np.log(res.y[sp] * res.phi[sp] * res.P / res.P0)
                         for sp in sp4])
        lnK = -nu @ np.array([res.dGf[sp] for sp in sp4]) / (R_J * res.T)
        r4 = float(np.max(np.abs(nu @ ln_a - lnK)))
        rep.add("pass" if r4 < tol else "fail", 4,
                [f"{len(textos)} reacciones independientes: máx "
                 f"|Σν ln(y φ P/P0) − ln K| = {r4:.1e}"]
                + [f"R{i + 1}: {t}  (ln K = {lk:.4g})"
                   for i, (t, lk) in enumerate(zip(textos, lnK))])
    else:
        rep.add("na", 4, "sin reacciones independientes")
    # 5
    rng = np.random.default_rng(semilla)
    peor = 0.0
    ref = np.log(y)
    no_conv = 0
    for _ in range(n_arranques):
        x0 = dict(zip(especies, 10.0 ** rng.uniform(-8, 1, len(especies))))
        try:
            otro = minimizar_gibbs(res.especies, res.dGf, res.n0, res.T, res.P,
                                   res.P0, res.phi, arranque=x0,
                                   restricciones=res.restricciones,
                                   condensadas=res.condensadas)
        except RuntimeError:          # el Newton se estancó desde ese arranque
            no_conv += 1
            continue
        peor = max(peor, float(np.max(np.abs(
            np.log([otro.y[sp] for sp in especies]) - ref))))
        if res.condensadas:                     # misma fase y mismos moles
            if set(otro.n_cond) != set(res.n_cond):
                peor = float("inf")
            else:
                peor = max([peor] + [abs(otro.n_cond[k] - v) / max(abs(v), 1e-300)
                                     for k, v in res.n_cond.items()])
    rep.add("fail" if peor >= 1e-6 else ("warn" if no_conv else "pass"), 5,
            f"robustez: {n_arranques} arranques al azar, máx |Δ ln y| = "
            f"{peor:.1e}"
            + (f"; {no_conv} no convergieron (el Newton se estancó; no es "
               "otra solución, G de un gas ideal es convexa)" if no_conv else ""))
    # 6
    rep.add("pass" if res.G_RT <= res.G_RT_alim + 1e-9 else "fail", 6,
            f"nG/RT baja de {res.G_RT_alim:.4f} a {res.G_RT:.4f}")
    # 7
    A_act, elem_act = matriz_atomica(especies)
    lleno = rango(A_act) == len(elem_act)      # si no, λ no es única
    carbono = [k for k in res.n_cond if _composicion(k) == {"C": 1}]
    if carbono:
        rep.add("pass", 7, f"{carbono[0]} presente: a_C = "
                f"{float(np.exp(-res.lam_RT['C'])):.6g}")
    elif "C" in res.lam_RT and lleno:
        a_c = float(np.exp(-res.lam_RT["C"]))
        rep.add("pass" if a_c < 1.0 else "warn", 7,
                f"a_C = {a_c:.3g}"
                + ("" if a_c < 1.0 else ": puede formarse carbono sólido y "
                   "este resultado solo-gas no sería estable"))
    else:
        rep.add("na", 7, "a_C no aplica")
    # 8
    if res.restricciones:
        r8 = max(abs(sum(c * res.n[sp] for sp, c in coef.items()) - rhs)
                 / max(abs(rhs), 1.0) for coef, rhs in res.restricciones)
        rep.add("pass" if r8 < tol else "fail", 8,
                f"{len(res.restricciones)} restricciones extra: máx "
                f"|Σ coef·n − rhs| = {r8:.1e}")
    # 9
    if res.condensadas:
        malas = [k for k in res.condensadas if k not in res.excluidas
                 and ((k in res.n_cond and res.n_cond[k] < -tol)
                      or (k not in res.n_cond and res.delta_fase[k] < -tol))]
        lineas = [f"prueba de fase: {len(res.n_cond)} de "
                  f"{len(res.condensadas)} fase(s) condensada(s) presentes"]
        lineas += [f"{k}: excluida (elemento ausente)" if k in res.excluidas else
                   f"{k}: Δ solo gas = {res.delta_fase_gas[k]:+.4f}, "
                   f"Δ final = {res.delta_fase[k]:+.2e}, "
                   + (f"n = {res.n_cond[k]:.6g}" if k in res.n_cond else "ausente")
                   for k in res.condensadas]
        rep.add("fail" if malas else "pass", 9, lineas)
    rep.verdict()
    return rep


def prueba_de_fase(res: ResultadoGibbs, especie: str,
                   dGf: "float | Mapping[str, float]") -> float:
    """Δ_k = ΔG°f,k/RT + Σ_e A_ek λ_e/RT de una fase condensada pura con los λ
    de `res`, sin recalcular (convención de SVA Ej. 13.14, λ_C/RT = 0.7635).
    Δ_k < 0: la fase baja G y debe incluirse; Δ_k ≥ 0: no se forma.
    `dGf` es ΔG°f [J/mol] o un mapeo que contiene a `especie`."""
    g = float(dGf[especie]) if isinstance(dGf, Mapping) else float(dGf)
    comp = _composicion(especie)
    faltan = [el for el in comp if el not in res.lam_RT]
    if faltan:
        raise KeyError(f"{especie}: sin λ para {', '.join(faltan)} "
                       f"(elemento ausente de la alimentación)")
    return g / (R_J * res.T) + sum(c * res.lam_RT[el] for el, c in comp.items())

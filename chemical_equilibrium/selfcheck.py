"""
selfcheck.py — Capa de autovalidación para ejercicios NUEVOS sin clave de
respuestas (Método A / avance de reacción).

Uso típico (ver plantilla_ejercicio.py):

    from selfcheck import Problem, solve_and_validate, validate_sweep

    prob = Problem(species=..., n0=..., nu=..., K=[...], P=..., P0=1.0,
                   atoms=..., inerts=[...], T=..., pressure_unit="bar")
    res = solve_and_validate(prob)     # resuelve y llama validate() al final

o, si ya tienes el resultado de solve_extents:

    validate(prob, res)

y tras un barrido de temperaturas:

    validate_sweep(casos, xi_por_T)    # casos: {T: [K_j]}, xi_por_T: {T: (xi_j,)}

Niveles del reporte (cada punto numerado 1-12):
    [PASS]  la comprobación se cumplió
    [WARN]  sospechoso: revísalo a mano (no invalida por sí solo)
    [FAIL]  violación dura: NO confíes en el resultado
    [info]  información para comparar con tu derivación a mano
    [n/a ]  no aplica a este problema

Puntos verificados
------------------
Sanidad de la entrada (errores de captura):
  1. balance de elementos por reacción (solve_and_validate SE NIEGA a
     resolver si falla, nombrando reacción y elemento)
  2. columnas nulas de nu con n0 > 0 no declaradas inertes
  3. número de K == número de reacciones
  4. delta_j y expresión de n_T(xi) — compárala con tu derivación
  5. expresión simbólica de cada K — compárala con la tuya
Calidad de la solución:
  6. residuos |K_calc - K|/K por reacción
  7. balance de átomos in/out, n_i >= 0, sum y_i = 1
  8. robustez: re-resuelve desde >= 5 arranques aleatorios factibles y
     exige el mismo xi (dispersión < 1e-6)
  9. contraste independiente Método A vs Método B (solo 1 reacción)
 10. avances en su límite físico (reactivo agotado -> conversión ~ 1)
Plausibilidad física (solo WARN):
 11. dirección de Le Chatelier en un barrido de T (validate_sweep)
 12. factor de presión (P/P0)^delta_j realmente usado
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import least_squares

from equilibrium import (
    AtomComposition, ExtentResult, StoichiometricTable, Ky_to_Kc,
    atom_totals, reaction_element_imbalance, solve_extents,
)

_TAGS = {"pass": "[PASS]", "warn": "[WARN]", "fail": "[FAIL]",
         "info": "[info]", "na": "[n/a ]"}
_IND = " " * 15          # continuación alineada bajo el texto del punto


# ----------------------------------------------------------------------------
# Especificación del problema y reporte
# ----------------------------------------------------------------------------

@dataclass
class Problem:
    """Todo lo que define un problema de Método A (lo mismo que recibe
    solve_extents, más metadatos que permiten validar)."""
    species: List[str]
    n0: Dict[str, float]
    nu: Sequence[Sequence[float]]
    K: Sequence[float]
    P: float
    P0: float = 1.0
    atoms: Optional[AtomComposition] = None
    inerts: Sequence[str] = ()        # especies inertes DECLARADAS (punto 2)
    T: Optional[float] = None         # temperatura del caso (puntos 9 y 12)
    pressure_unit: str = "bar"        # unidad de P y P0 (punto 9)
    label: str = ""


@dataclass
class Report:
    """Contador de niveles; imprime cada punto al agregarlo."""
    n_pass: int = 0
    n_warn: int = 0
    n_fail: int = 0

    @property
    def ok(self) -> bool:
        """True si no hubo ningún FAIL (los WARN se revisan a mano)."""
        return self.n_fail == 0

    def add(self, level: str, num, text) -> None:
        lines = [text] if isinstance(text, str) else list(text)
        if level == "pass":
            self.n_pass += 1
        elif level == "warn":
            self.n_warn += 1
        elif level == "fail":
            self.n_fail += 1
        print(f"    {_TAGS[level]} {num:>2}. {lines[0]}")
        for extra in lines[1:]:
            print(_IND + extra)

    def verdict(self) -> str:
        v = "FAIL" if self.n_fail else ("WARN" if self.n_warn else "PASS")
        print(f"  Veredicto: {v}  ({self.n_pass} PASS, {self.n_warn} WARN, "
              f"{self.n_fail} FAIL)")
        return v


def _box(title: str) -> None:
    print(f"\n  {'-' * 62}\n  {title}\n  {'-' * 62}")


# ----------------------------------------------------------------------------
# Punto 1 — balance de elementos (también usado para negarse a resolver)
# ----------------------------------------------------------------------------

def element_problems(problem: Problem) -> List[str]:
    """Desbalances por reacción: ['R2: elemento O neto +1', ...]; [] si todo
    balancea (o si no hay dict de átomos con qué verificar)."""
    out: List[str] = []
    if problem.atoms is None:
        return out
    for j, row in enumerate(np.atleast_2d(np.asarray(problem.nu, float))):
        net = reaction_element_imbalance(row, problem.species, problem.atoms)
        for el in sorted(net):
            if abs(net[el]) > 1e-9:
                out.append(f"R{j + 1}: elemento {el} neto {net[el]:+g}")
    return out


# ----------------------------------------------------------------------------
# Punto 8 — re-solver independiente (mismas ecuaciones, arranques aleatorios)
# ----------------------------------------------------------------------------

_EPS2, _PEN = 1e-28, 1e6     # piso suave de moles: (1e-14)^2


def _make_residuals(nu, n0_vec, lnK, delta, ln_Pr):
    """Residuos en forma log como en equilibrium.py, pero con un piso de
    moles SUAVE, n_pos = (n + sqrt(n^2 + 4*eps^2))/2, en lugar del recorte
    duro max(n, 1e-12). Mismas raíces (el suavizado desplaza n en ~eps^2/n,
    despreciable para n físicos), pero el gradiente nunca se anula cerca de
    la frontera n = 0, así que un arranque aleatorio no se queda atascado
    en el 'escalón' plano del recorte duro."""
    def _parts(xi):
        n = n0_vec + nu.T @ xi
        s = np.sqrt(n * n + 4.0 * _EPS2)
        n_pos = np.maximum(0.5 * (n + s), 1e-300)
        w = 0.5 * (1.0 + n / s)                 # d n_pos / d n
        return n, n_pos, w

    def residuals(xi):
        n, n_pos, _ = _parts(xi)
        r = (nu @ np.log(n_pos) - delta * np.log(n_pos.sum())
             + delta * ln_Pr - lnK)
        return np.concatenate([r, _PEN * np.minimum(n, 0.0)])

    def jacobian(xi):
        n, n_pos, w = _parts(xi)
        J_K = (nu * (w / n_pos)) @ nu.T - np.outer(delta, nu @ w) / n_pos.sum()
        J_pen = _PEN * np.where(n[:, None] < 0.0, nu.T, 0.0)
        return np.vstack([J_K, J_pen])

    return residuals, jacobian


def _random_feasible_xi(nu, n0_vec, rng) -> np.ndarray:
    """Avanza las reacciones una por una (orden aleatorio) dentro de los
    límites de moles disponibles: siempre produce un xi factible."""
    n_rxn = nu.shape[0]
    n_feed = max(float(n0_vec.sum()), 1.0)
    xi = np.zeros(n_rxn)
    n = n0_vec.astype(float).copy()
    for j in rng.permutation(n_rxn):
        cons = np.where(nu[j] < 0)[0]
        prod = np.where(nu[j] > 0)[0]
        hi = min([n[i] / -nu[j, i] for i in cons], default=n_feed)
        lo = -min([n[i] / nu[j, i] for i in prod], default=n_feed)
        xi[j] = lo + rng.random() * (hi - lo)
        n = n + nu[j] * xi[j]
    return xi


def _check_robustness(rep: Report, nu, n0_vec, K, delta, P, P0, xi_ref,
                      n_obj: int = 5, max_try: int = 25) -> None:
    n_rxn = nu.shape[0]
    K_arr = np.asarray(K, dtype=float)
    if K_arr.shape != (n_rxn,):
        rep.add("na", 8, "robustez: omitida (número de K != número de "
                         "reacciones, ver punto 3)")
        return
    residuals, jacobian = _make_residuals(nu, n0_vec, np.log(K_arr), delta,
                                          np.log(P / P0))
    x_scale = np.maximum(np.abs(np.asarray(xi_ref, dtype=float)), 1e-3)

    def _resolve_from(x0):
        """Hasta 4 rondas de least_squares (cada reinicio refresca la región
        de confianza); devuelve (xi, convergió)."""
        x = np.asarray(x0, dtype=float)
        for _ in range(4):
            try:
                sol = least_squares(residuals, x, jac=jacobian, method="trf",
                                    x_scale=x_scale, xtol=1e-15, ftol=1e-15,
                                    gtol=1e-15, max_nfev=3000)
            except Exception:
                return x, False
            x = sol.x
            n = n0_vec + nu.T @ x
            if (np.max(np.abs(residuals(x)[:n_rxn])) < 5e-7
                    and float(np.min(n)) > -1e-9):
                return x, True
        return x, False

    rng = np.random.default_rng(20260901)
    sols: List[np.ndarray] = []
    attempts = 0
    while len(sols) < n_obj and attempts < max_try:
        attempts += 1
        x, ok = _resolve_from(_random_feasible_xi(nu, n0_vec, rng))
        if ok:
            sols.append(x)

    todos = sols + [np.asarray(xi_ref, dtype=float)]
    spread = max(float(np.max(np.abs(a - b))) for a in todos for b in todos)

    if len(sols) >= n_obj and spread < 1e-6:
        rep.add("pass", 8,
                f"robustez: {len(sols)} arranques aleatorios factibles "
                f"convergen al mismo xi (dispersión máx = {spread:.1e} < 1e-6, "
                f"{attempts} intentos)")
    elif spread >= 1e-6:
        clusters: List[Tuple[np.ndarray, int]] = []
        for x in todos:
            for k, (cx, cnt) in enumerate(clusters):
                if np.max(np.abs(x - cx)) < 1e-6:
                    clusters[k] = (cx, cnt + 1)
                    break
            else:
                clusters.append((x, 1))
        lines = [f"robustez: soluciones DISTINTAS según el arranque "
                 f"(dispersión = {spread:.1e}) — posibles equilibrios "
                 f"múltiples o mal condicionamiento:"]
        for k, (cx, cnt) in enumerate(clusters):
            lines.append(f"solución {k + 1} (x{cnt}): xi = ["
                         + ", ".join(f"{v:.6f}" for v in cx) + "]")
        rep.add("warn", 8, lines)
    else:
        rep.add("warn", 8,
                f"robustez: solo {len(sols)}/{n_obj} arranques aleatorios "
                f"convergieron en {attempts} intentos (los que convergieron "
                f"coinciden, dispersión {spread:.1e}) — posible mal "
                f"condicionamiento")


# ----------------------------------------------------------------------------
# Punto 9 — contraste Método B (solo 1 reacción)
# ----------------------------------------------------------------------------

def _check_method_b(rep: Report, problem: Problem, result: ExtentResult,
                    nu, n0_vec) -> None:
    if nu.shape[0] != 1:
        rep.add("na", 9, f"contraste Método B: solo aplica a problemas de "
                         f"1 reacción (aquí hay {nu.shape[0]})")
        return
    if problem.T is None:
        rep.add("na", 9, "contraste Método B: requiere problem.T para "
                         "convertir K -> Kc")
        return
    row = nu[0]
    sp = problem.species
    reactants = [i for i, v in enumerate(row) if v < 0 and n0_vec[i] > 0]
    if not reactants:
        rep.add("na", 9, "contraste Método B: ningún reactivo con n0 > 0")
        return
    i_lim = min(reactants, key=lambda i: n0_vec[i] / -row[i])
    A, a, nA0 = sp[i_lim], -row[i_lim], n0_vec[i_lim]
    nu_b = {sp[i]: row[i] / a for i in range(len(sp)) if row[i] != 0}
    theta = {sp[i]: n0_vec[i] / nA0 for i in range(len(sp)) if n0_vec[i] > 0}
    try:
        tbl = StoichiometricTable(nu=nu_b, theta=theta, limiting=A,
                                  T0=problem.T, P0=problem.P,
                                  pressure_unit=problem.pressure_unit,
                                  variable_volume=True, n_A0=nA0,
                                  atoms=problem.atoms)
        # Kc de la reacción NORMALIZADA (dividida entre a = |nu_A|):
        #   K' = K^(1/a), delta' = delta/a, Kc' = K' * C_T0(P0)^delta'.
        Kc = Ky_to_Kc(float(problem.K[0]) ** (1.0 / a), float(row.sum()) / a,
                      problem.T, P0=problem.P0, unit=problem.pressure_unit)
        Xe = tbl.solve_Xe(Kc)
    except Exception as exc:
        rep.add("warn", 9, f"contraste Método B no se pudo ejecutar: {exc}")
        return
    n_b = tbl.moles(Xe)
    n_tot = sum(n_b.values())
    y_b = {s: v / n_tot for s, v in n_b.items()}
    y_a = result.mole_fractions
    dmax = max(abs(y_a[s] - y_b.get(s, 0.0)) for s in sp)
    X_a = (nA0 - result.moles[A]) / nA0
    head = (f"contraste Método B (limitante {A}): X_e(B) = {Xe:.6f} vs "
            f"X(A) = {X_a:.6f};  máx |y_A - y_B| = {dmax:.1e}")
    if dmax < 1e-6:
        rep.add("pass", 9, head)
    else:
        rep.add("warn", 9, [head + "  — los métodos NO coinciden:",
                            "y por especie (A / B): " + ", ".join(
                                f"{s}: {y_a[s]:.6f}/{y_b.get(s, 0.0):.6f}"
                                for s in sp)])


# ----------------------------------------------------------------------------
# validate — el reporte completo (puntos 1-10 y 12; el 11 vive en el barrido)
# ----------------------------------------------------------------------------

def validate(problem: Problem, result: ExtentResult,
             n_starts: int = 5, k_tol: float = 1e-6) -> Report:
    """Corre TODAS las comprobaciones sobre un resultado de solve_extents y
    imprime el reporte numerado. Devuelve el Report (report.ok = sin FAIL).

    n_starts: arranques aleatorios del punto 8 (robustez).
    k_tol   : tolerancia relativa |K_calc-K|/K del punto 6."""
    nu = np.atleast_2d(np.asarray(problem.nu, dtype=float))
    n_rxn, n_sp = nu.shape
    rep = Report()

    title = "AUTOVALIDACIÓN (selfcheck.validate)"
    extras = [s for s in (problem.label,
                          f"T = {problem.T:g} K" if problem.T is not None
                          else "") if s]
    if extras:
        title += " — " + ", ".join(extras)
    _box(title)

    if n_sp != len(problem.species):
        rep.add("fail", 1, f"nu tiene {n_sp} columnas pero hay "
                           f"{len(problem.species)} especies — revisa la matriz")
        rep.verdict()
        return rep

    n0_vec = np.array([float(problem.n0.get(s, 0.0)) for s in problem.species])
    delta = nu.sum(axis=1)

    # ---------------- Sanidad de la entrada ----------------
    print("  Sanidad de la entrada:")

    # 1. balance de elementos por reacción
    if problem.atoms is None:
        rep.add("warn", 1, "balance de elementos: sin dict `atoms`, "
                           "no se puede verificar")
    else:
        try:
            probs = element_problems(problem)
        except KeyError as exc:
            probs = [f"especie sin composición atómica: {exc}"]
        if probs:
            rep.add("fail", 1, ["reacciones NO balanceadas — no confíes en "
                                "nada de lo que sigue:"] + probs)
        else:
            els = sorted({el for s in problem.species
                          for el in problem.atoms.get(s, {})})
            rep.add("pass", 1, f"balance de elementos ({', '.join(els)}): "
                               f"R1..R{n_rxn} balanceadas")

    # 2. columnas nulas de nu vs inertes declarados
    zero_cols = [problem.species[i] for i in range(n_sp)
                 if np.all(np.abs(nu[:, i]) < 1e-14)]
    inerts = set(problem.inerts or ())
    sospechosas = [s for s in zero_cols
                   if problem.n0.get(s, 0.0) > 0 and s not in inerts]
    mal_flag = [s for s in inerts
                if s in problem.species and s not in zero_cols]
    if sospechosas or mal_flag:
        lines = []
        if sospechosas:
            lines.append("especies con n0 > 0 y columna NULA en nu sin "
                         "declarar inertes: " + ", ".join(sospechosas)
                         + "  (¿falta en alguna reacción?)")
        if mal_flag:
            lines.append("declaradas inertes pero con coeficientes != 0: "
                         + ", ".join(mal_flag))
        rep.add("warn", 2, lines)
    elif zero_cols:
        rep.add("pass", 2, "columna nula en nu: " + ", ".join(zero_cols)
                           + " — declarada(s) inerte(s), consistente")
    else:
        rep.add("pass", 2, "todas las especies participan en alguna reacción")

    # 3. número de K vs número de reacciones
    if len(problem.K) == n_rxn:
        rep.add("pass", 3, f"número de K ({len(problem.K)}) == número de "
                           f"reacciones ({n_rxn})")
    else:
        rep.add("warn", 3, f"número de K ({len(problem.K)}) != número de "
                           f"reacciones ({n_rxn})")

    # 4. delta_j y n_T(xi)
    dtxt = ",  ".join(f"R{j + 1}: {d:+g}" for j, d in enumerate(delta))
    terms = "".join(f" {'+' if d > 0 else '-'} {abs(d):g}·xi_{j + 1}"
                    for j, d in enumerate(delta) if abs(d) > 1e-12)
    expr = f"n_T = {float(n0_vec.sum()):g}{terms}"
    if not terms:
        expr += "   (constante: todas las delta_j = 0)"
    rep.add("info", 4, [f"delta_j:  {dtxt}", expr])

    # 5. expresiones de equilibrio simbólicas
    lines5 = []
    for j in range(n_rxn):
        num = " ".join(f"y_{problem.species[i]}"
                       + (f"^{nu[j, i]:g}" if nu[j, i] != 1 else "")
                       for i in range(n_sp) if nu[j, i] > 0)
        den = " ".join(f"y_{problem.species[i]}"
                       + (f"^{-nu[j, i]:g}" if nu[j, i] != -1 else "")
                       for i in range(n_sp) if nu[j, i] < 0)
        s = f"K{j + 1} = {num if num else '1'}"
        if den:
            s += f" / ({den})"
        if abs(delta[j]) > 1e-12:
            s += f" * (P/P0)^{delta[j]:g}"
        lines5.append(s)
    rep.add("info", 5, lines5)

    # ---------------- Calidad de la solución ----------------
    print("  Calidad de la solución:")

    # 6. residuos |K_calc - K|/K
    rel = np.abs(result.K_calc - result.K_given) / result.K_given
    txt6 = ("residuos |K_calc-K|/K: ["
            + ", ".join(f"{r:.1e}" for r in rel)
            + f"]  (máx = {rel.max():.1e} < {k_tol:g})")
    if result.converged and rel.max() < k_tol:
        rep.add("pass", 6, txt6)
    else:
        lines = [txt6.replace(f" < {k_tol:g}", f" — NO cumple {k_tol:g}")]
        if not result.converged:
            lines.append(f"solver: {result.message}")
        rep.add("fail", 6, lines)

    # 7. átomos in/out, n_i >= 0, sum y = 1
    ok7, lines7 = True, []
    if problem.atoms is not None:
        in_tot = atom_totals({s: problem.n0.get(s, 0.0)
                              for s in problem.species}, problem.atoms)
        out_tot = atom_totals(result.moles, problem.atoms)
        worst = max((abs(in_tot.get(e, 0.0) - out_tot.get(e, 0.0))
                     for e in set(in_tot) | set(out_tot)), default=0.0)
        ok7 &= worst < 1e-8
        lines7.append("átomos in vs out: máx |Δ| = "
                      f"{worst:.1e} mol  ("
                      + ", ".join(f"{e}: {in_tot.get(e, 0.0):g} -> "
                                  f"{out_tot.get(e, 0.0):.6g}"
                                  for e in sorted(set(in_tot) | set(out_tot)))
                      + ")")
    else:
        lines7.append("átomos in vs out: omitido (sin dict `atoms`)")
    nmin = float(result.n.min())
    ok7 &= nmin > -1e-9
    lines7.append(f"n_i >= 0: mín n_i = {nmin:.2e} mol")
    sy = float(result.y.sum())
    ok7 &= abs(sy - 1.0) < 1e-9
    lines7.append(f"sum y_i = {sy:.12f}")
    rep.add("pass" if ok7 else "fail", 7, lines7)

    # 8. robustez multi-arranque
    _check_robustness(rep, nu, n0_vec, problem.K, delta,
                      problem.P, problem.P0, result.xi,
                      n_obj=n_starts, max_try=max(25, 3 * n_starts))

    # 9. contraste con Método B (1 sola reacción)
    _check_method_b(rep, problem, result, nu, n0_vec)

    # 10. avances en su límite físico
    lines10 = []
    for j in range(n_rxn):
        cons = [i for i in range(n_sp) if nu[j, i] < 0]
        if not cons:
            continue
        slack, i_min = min((result.n[i] / -nu[j, i], i) for i in cons)
        if slack < 1e-6:
            s = problem.species[i_min]
            if n0_vec[i_min] > 0:
                X_s = (n0_vec[i_min] - result.n[i_min]) / n0_vec[i_min]
                x_txt = f"conversión de {s} = {X_s:.10f} ≈ 1"
            else:
                x_txt = f"{s} no estaba en la alimentación"
            lines10.append(f"R{j + 1} en su límite físico: {s} agotado "
                           f"(n = {result.n[i_min]:.2e} mol) -> {x_txt}")
    if lines10:
        rep.add("info", 10, lines10)
    else:
        rep.add("pass", 10, "ningún avance en su límite físico "
                            "(todas las holguras > 1e-6 mol)")

    # ---------------- Plausibilidad física ----------------
    print("  Plausibilidad física:")

    # 11. (vive en validate_sweep: necesita el barrido completo)
    rep.add("na", 11, "Le Chatelier: requiere el barrido de T completo — "
                      "llama validate_sweep(casos, xi_por_T) al final")

    # 12. factor de presión realmente usado
    Pr = problem.P / problem.P0
    if abs(Pr - 1.0) < 1e-12:
        rep.add("info", 12, f"P = P0 = {problem.P:g} {problem.pressure_unit} "
                            f"-> factor (P/P0)^delta_j = 1 en todas las "
                            f"reacciones")
    else:
        lines12 = [f"P/P0 = {problem.P:g}/{problem.P0:g} "
                   f"{problem.pressure_unit} = {Pr:g}   "
                   f"(¡P y P0 deben estar en la MISMA unidad!)"]
        con_delta = [f"R{j + 1}: (P/P0)^{d:g} = {Pr ** d:.6g}"
                     for j, d in enumerate(delta) if abs(d) > 1e-12]
        lines12 += con_delta or ["todas las delta_j = 0: la presión no "
                                 "desplaza ninguna reacción"]
        rep.add("info", 12, lines12)

    rep.verdict()
    return rep


# ----------------------------------------------------------------------------
# Punto 11 — Le Chatelier sobre un barrido de temperaturas
# ----------------------------------------------------------------------------

def _trend(vals: Sequence[float], tol: float = 1e-12) -> str:
    if all(b > a + tol for a, b in zip(vals, vals[1:])):
        return "creciente"
    if all(b < a - tol for a, b in zip(vals, vals[1:])):
        return "decreciente"
    return "no monótona"


def validate_sweep(K_por_T: Dict[float, Sequence[float]],
                   xi_por_T: Dict[float, Sequence[float]],
                   reaction_names: Optional[Sequence[str]] = None) -> Report:
    """Punto 11: para cada reacción compara la tendencia de K(T) con la de
    xi(T) a lo largo del barrido. K y xi en dirección OPUESTA monótona =>
    WARN (revisa K(T) o atribúyelo al acoplamiento entre reacciones)."""
    rep = Report()
    _box("AUTOVALIDACIÓN (selfcheck.validate_sweep) — Le Chatelier")
    temps = sorted(set(K_por_T) & set(xi_por_T))
    if len(temps) < 2:
        rep.add("na", 11, "se necesitan >= 2 temperaturas para ver tendencias")
        rep.verdict()
        return rep
    print("  Temperaturas: " + ", ".join(f"{t:g} K" for t in temps))
    n_rxn = min(len(K_por_T[t]) for t in temps)
    for j in range(n_rxn):
        Ks = [float(K_por_T[t][j]) for t in temps]
        xis = [float(xi_por_T[t][j]) for t in temps]
        tK, tx = _trend(Ks), _trend(xis)
        name = (reaction_names[j] if reaction_names else f"R{j + 1}")
        head = f"{name + ':':<5s}"
        lines = [head + "K : " + " -> ".join(f"{v:.4e}" for v in Ks)
                 + f"   ({tK})",
                 " " * len(head) + "xi: " + " -> ".join(f"{v:.5f}" for v in xis)
                 + f"   ({tx})"]
        if {tK, tx} == {"creciente", "decreciente"}:
            lines.append("K y xi van en direcciones OPUESTAS: revisa K(T), o "
                         "atribúyelo (a mano) al acoplamiento entre reacciones")
            rep.add("warn", 11, lines)
        elif "no monótona" in (tK, tx):
            lines.append("tendencia no monótona: sin conclusión directa "
                         "(típico de reacciones acopladas compitiendo)")
            rep.add("info", 11, lines)
        else:
            lines.append("coherente con Le Chatelier (K y xi en la misma "
                         "dirección)")
            rep.add("pass", 11, lines)
    rep.verdict()
    return rep


# ----------------------------------------------------------------------------
# Resolver + validar en un solo paso
# ----------------------------------------------------------------------------

def solve_and_validate(problem: Problem, **solver_kwargs) -> ExtentResult:
    """Resuelve por Método A y ejecuta validate() automáticamente al final.

    Antes de resolver corre el punto 1 (balance de elementos) y SE NIEGA a
    resolver (ValueError) si alguna reacción está desbalanceada, nombrando
    la reacción y el elemento."""
    if problem.atoms is not None:
        try:
            probs = element_problems(problem)
        except KeyError as exc:
            raise ValueError(f"Composición atómica faltante: {exc}") from None
        if probs:
            raise ValueError("Reacción(es) sin balancear — me niego a "
                             "resolver:\n  " + "\n  ".join(probs))
    res = solve_extents(problem.species, problem.n0, problem.nu, problem.K,
                        problem.P, P0=problem.P0, atoms=problem.atoms,
                        **solver_kwargs)
    validate(problem, res)
    return res

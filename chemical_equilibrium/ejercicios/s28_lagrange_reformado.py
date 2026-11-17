"""Sesión 28 — Reformado de metano resuelto como sistema de Lagrange.

Koretsky ec. 9.45–9.47: m + l ecuaciones (CNS en GAMS)
    Δg°f,i/RT + ln(y_i P/P°) + Σ_j β_ij λ_j/RT = 0     (m = 5 especies)
    Σ_i n_i β_ij = b_j                               (l = 3 elementos)
incógnitas ln n_i y λ_j/RT. Se resuelven tal cual con scipy.optimize.root,
sin el Newton de core/gibbs.py, y se contrastan con minimizar_gibbs.

Ejemplo de la diapositiva: reformado 4:1 a 1 bar, 600–1100 K. OJO: la
diapositiva escribe los balances C = 2, H = 14, O = 3, que son los de
2 CH4 + 3 H2O (SVA Ej. 13.14), no los de H2O:CH4 = 4:1 (C = 1, H = 12, O = 4).
Se resuelven los dos casos. La maestra confirmó 4:1 (17-nov-2026): es el
caso que se entrega; el 2:3 queda como referencia.

Nota de la maestra: con los λ de la fase gas se evalúa si otra fase baja G.
Aquí, grafito: Δ_C = g_C(s)/RT + λ_C/RT = λ_C/RT (ΔfG = 0, JANAF C-002). Si
Δ_C < 0 se agrega C(s) con la ecuación g_C(s)/RT + λ_C/RT = 0 y su n_C en el
balance de C, y se recalcula; si n_C sale < 0, se quita.
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto
import core.entorno  # noqa: F401

import numpy as np
from scipy.optimize import root

from core.gibbs import cargar_dgf, minimizar_gibbs
from core.independencia import matriz_atomica, transpuesta
from core.thermo import R_J

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
ESPECIES = ["CH4", "H2O", "CO", "CO2", "H2"]
CASOS = {"4:1 (enunciado)": {"CH4": 1.0, "H2O": 4.0},
         "2:3 (balances de la diapositiva)": {"CH4": 2.0, "H2O": 3.0}}
P, P0 = 1.0, 1.0                                   # bar
TEMPERATURAS = [600, 700, 800, 900, 1000, 1100]    # K
DEMO = ("demostración H2O:CH4 = 1:1", {"CH4": 1.0, "H2O": 1.0}, 800)

A_ex, ELEMENTOS = matriz_atomica(ESPECIES)
BETA = np.array(transpuesta(A_ex), dtype=float)    # β: elementos × especies
COMP_C = np.array([1.0 if el == "C" else 0.0 for el in ELEMENTOS])   # C(s)


def titulo(t: str) -> None:
    print(f"\n{'=' * 84}\n{t}\n{'=' * 84}")


def resolver_lagrange(dGf, n0, T, con_grafito=False, x0=None):
    """Sistema de Koretsky 9.45–9.47 (y la ecuación de C(s) si va). Devuelve
    dict con n, y, λ/RT, n_C(s) y el máximo residuo."""
    g = np.array([dGf[sp] for sp in ESPECIES]) / (R_J * T)
    b = BETA @ np.array([n0.get(sp, 0.0) for sp in ESPECIES])
    m, l = len(ESPECIES), len(ELEMENTOS)

    def F(x):
        lnn, lam = x[:m], x[m:m + l]
        n = np.exp(lnn)
        nC = x[m + l] if con_grafito else 0.0
        f = [g + lnn - np.log(n.sum()) + np.log(P / P0) + BETA.T @ lam,      # 9.45
             (BETA @ n + nC * COMP_C - b) / b]                               # 9.47
        if con_grafito:
            f.append([dGf["C(s)"] / (R_J * T) + COMP_C @ lam])               # C(s)
        return np.concatenate(f)

    if x0 is None:
        x0 = np.concatenate([np.full(m, np.log(b.sum() / (10.0 * m))), np.zeros(l)])
    if con_grafito and len(x0) == m + l:
        x0 = np.append(x0, 0.0)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        s = root(F, x0, method="hybr", options={"xtol": 1e-14})
        if not s.success or np.abs(s.fun).max() > 1e-10:
            s = root(F, x0, method="lm", options={"xtol": 1e-15, "ftol": 1e-15})
    n = np.exp(s.x[:m])
    return {"n": dict(zip(ESPECIES, n)), "y": dict(zip(ESPECIES, n / n.sum())),
            "lam": dict(zip(ELEMENTOS, s.x[m:m + l])), "x": s.x,
            "n_C": float(s.x[m + l]) if con_grafito else 0.0,
            "res": float(np.abs(s.fun).max()), "ok": bool(s.success)}


def resolver_con_prueba(dGf, n0, T):
    """Gas solo -> Δ_C -> (si Δ_C < 0) con C(s) -> (si n_C < 0) sin C(s)."""
    gas = resolver_lagrange(dGf, n0, T)
    delta = dGf["C(s)"] / (R_J * T) + gas["lam"]["C"]
    final, decision = gas, "no se forma (Δ_C ≥ 0)"
    if delta < 0:
        con = resolver_lagrange(dGf, n0, T, con_grafito=True, x0=gas["x"])
        if con["n_C"] >= 0:
            final, decision = con, f"se incluye C(s): n_C = {con['n_C']:.5f} mol"
        else:
            decision = f"n_C = {con['n_C']:.3g} < 0: se quita C(s)"
    return gas, delta, final, decision


# %% 1. El sistema
titulo("1. Sistema de Lagrange (Koretsky 9.45–9.47): m = 5 especies, l = 3 elementos")
print("  " + "\n  ".join(
    f"{sp:<4}: Δg°f/RT + ln(y_{sp}·P/P°) + "
    + " + ".join(f"{int(BETA[j, i])}·λ_{el}" for j, el in enumerate(ELEMENTOS) if BETA[j, i])
    + " = 0" for i, sp in enumerate(ESPECIES)))
print("  " + "\n  ".join(
    f"{el:<4}: " + " + ".join(f"{int(BETA[j, i])}·n_{sp}" for i, sp in enumerate(ESPECIES) if BETA[j, i])
    + f" = b_{el}" for j, el in enumerate(ELEMENTOS)))
for nombre, n0 in CASOS.items():
    b = BETA @ np.array([n0.get(sp, 0.0) for sp in ESPECIES])
    print(f"  {nombre}: b = " + ", ".join(f"{el} {v:g}" for el, v in zip(ELEMENTOS, b)))
print("  INCONSISTENCIA: la diapositiva dice '4:1' pero escribe C = 2, H = 14, O = 3,")
print("  que son los balances de 2 CH4 + 3 H2O. Con 4:1 son C = 1, H = 12, O = 4.")

# %% 2. Barrido para los dos casos
DGF = {T: cargar_dgf(ESPECIES + ["C(s)"], float(T))[0] for T in TEMPERATURAS}
RESULTADOS = {}
for nombre, n0 in CASOS.items():
    for T in TEMPERATURAS:
        gas, delta, final, decision = resolver_con_prueba(DGF[T], n0, float(T))
        rand = minimizar_gibbs(ESPECIES, DGF[T], n0, float(T), P, P0,
                               condensadas=["C(s)"])
        dif = max(abs(np.log(final["y"][sp]) - np.log(rand.y[sp])) for sp in ESPECIES)
        RESULTADOS[(nombre, T)] = {"final": final, "delta": delta, "decision": decision,
                                   "rand": rand, "dif": dif}
    titulo(f"2. {nombre}, {P:g} bar: composición, λ_j/RT y prueba de fase de C(s)")
    print(f"{'T':>5}" + "".join(f"{'y_' + sp:>11}" for sp in ESPECIES)
          + "".join(f"{'λ_' + el + '/RT':>10}" for el in ELEMENTOS)
          + f"{'Δ_C':>9}{'|res|':>9}{'vs RAND':>9}  C(s)")
    for T in TEMPERATURAS:
        r = RESULTADOS[(nombre, T)]
        f = r["final"]
        print(f"{T:>5}" + "".join(f"{f['y'][sp]:>11.4e}" for sp in ESPECIES)
              + "".join(f"{f['lam'][el]:>10.4f}" for el in ELEMENTOS)
              + f"{r['delta']:>9.4f}{f['res']:>9.1e}{r['dif']:>9.1e}  {r['decision']}")

# %% 3. Demostración del paso de la maestra (fuera del enunciado)
nombre, n0, T = DEMO
gas, delta, final, decision = resolver_con_prueba(DGF[T], n0, float(T))
rand = minimizar_gibbs(ESPECIES, DGF[T], n0, float(T), P, P0, condensadas=["C(s)"])
titulo(f"3. {nombre} a {T} K (fuera del enunciado): la fase sí entra")
print(f"  gas solo: λ_C/RT = {gas['lam']['C']:+.4f} -> Δ_C = {delta:+.4f} < 0: el grafito baja G")
print(f"  con C(s): {decision}; λ_C/RT = {final['lam']['C']:+.2e} (= −g_C(s)/RT = 0)")
print(f"  contraste minimizar_gibbs(condensadas=['C(s)']): n_C = {rand.n_cond.get('C(s)', 0):.5f}, "
      f"máx |Δ ln y| = {max(abs(np.log(final['y'][s]) - np.log(rand.y[s])) for s in ESPECIES):.1e}")
print("  " + ", ".join(f"y_{sp} = {final['y'][sp]:.4e}" for sp in ESPECIES))

"""Sesión 26, Ejercicio 1 — Minimización de Gibbs: barrido en T.

Reformado de metano con H2O:CH4 = 4:1 a 1 bar, de 600 a 1100 K. Especies
CH4, H2O, CO, CO2, H2. ΔG°f de JANAF (data/janaf_dgf.json) en cada T.

Por cada T: minimización de G (core/gibbs.py), validación de 7 puntos,
contraste con el Método A (solve_extents con las K de los mismos ΔG°f) y
la prueba de fase del carbono, Δ_C = λ_C/RT (= −ln a_C): si fuera negativa,
el grafito bajaría G y el resultado solo-gas no sería estable.
"""

import contextlib
import io
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto
import core.entorno  # noqa: F401

import numpy as np

from core import SALIDAS
from core.equilibrium import solve_extents
from core.gibbs import (cargar_dgf, minimizar_gibbs, prueba_de_fase,
                        reacciones_independientes, validar_gibbs)
from core.thermo import R_J

GUARDAR = __name__ == "__main__"     # importado desde tests/: no escribe archivos

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
ESPECIES = ["CH4", "H2O", "CO", "CO2", "H2"]
ALIMENTACION = {"CH4": 1.0, "H2O": 4.0}   # mol, H2O:CH4 = 4:1
P = 1.0                                   # bar
TEMPERATURAS = [600, 700, 800, 900, 1000, 1100]   # K


def titulo(t: str) -> None:
    print(f"\n{'=' * 72}\n{t}\n{'=' * 72}")


# %% 1. ΔG°f de JANAF en cada T
DGF, FUENTES = {}, {}
for T in TEMPERATURAS:
    DGF[T], FUENTES[T] = cargar_dgf(ESPECIES + ["C(s)"], float(T))

titulo("1. ΔG°f (J/mol) de JANAF")
print(f"{'T (K)':>6}" + "".join(f"{sp:>11}" for sp in ESPECIES))
for T in TEMPERATURAS:
    print(f"{T:>6}" + "".join(f"{DGF[T][sp]:>11.0f}" for sp in ESPECIES))
print("Tablas: " + ", ".join(f"{sp} {FUENTES[TEMPERATURAS[0]][sp]}" for sp in ESPECIES)
      + "; C(s) JANAF C-002 (grafito, referencia)")

# %% 2. Minimización de Gibbs y validación en cada T
RESULTADOS, VALIDACION = {}, {}
for T in TEMPERATURAS:
    r = minimizar_gibbs(ESPECIES, DGF[T], ALIMENTACION, float(T), P)
    with contextlib.redirect_stdout(io.StringIO()):
        VALIDACION[T] = validar_gibbs(r)
    RESULTADOS[T] = r

titulo(f"2. Composición de equilibrio, {P:g} bar, H2O:CH4 = 4:1")
print(f"{'T (K)':>6}" + "".join(f"{'y_' + sp:>12}" for sp in ESPECIES)
      + f"{'n_T':>9}{'validación':>14}")
for T, r in RESULTADOS.items():
    v = VALIDACION[T]
    print(f"{T:>6}" + "".join(f"{r.y[sp]:>12.4e}" for sp in ESPECIES)
          + f"{r.n_total:>9.4f}{f'{v.n_pass} PASS {v.n_fail} FAIL':>14}")

# %% 3. Multiplicadores y prueba de fase del carbono
DELTA_C = {T: prueba_de_fase(r, "C(s)", DGF[T]) for T, r in RESULTADOS.items()}

titulo("3. λ_j/RT y prueba de fase de C(s): Δ_C = λ_C/RT = −ln a_C")
print(f"{'T (K)':>6}{'λ_C/RT':>10}{'λ_H/RT':>10}{'λ_O/RT':>10}{'Δ_C':>10}{'a_C':>9}  C(s)")
for T, r in RESULTADOS.items():
    print(f"{T:>6}{r.lam_RT['C']:>10.4f}{r.lam_RT['H']:>10.4f}{r.lam_RT['O']:>10.4f}"
          f"{DELTA_C[T]:>10.4f}{np.exp(-DELTA_C[T]):>9.4f}  "
          + ("no se forma (Δ_C > 0)" if DELTA_C[T] > 0 else "SE FORMARÍA (Δ_C < 0)"))

# %% 4. Contraste con el Método A (mismos ΔG°f -> mismas K)
METODO_A = {}
for T, r in RESULTADOS.items():
    sp, nu, textos = reacciones_independientes(r)
    K = np.exp(-nu @ np.array([DGF[T][s] for s in sp]) / (R_J * T))
    m = solve_extents(sp, ALIMENTACION, nu, K, P)
    METODO_A[T] = (m.converged, max(abs(np.log(r.y[s]) - np.log(m.mole_fractions[s]))
                                    for s in sp))

titulo("4. Contraste con el Método A (solve_extents)")
print("Reacciones: " + "; ".join(textos))
for T, (conv, d) in METODO_A.items():
    print(f"  {T} K: Método A {'convergió' if conv else 'NO convergió'}, "
          f"máx |Δ ln y| Gibbs vs Método A = {d:.1e}")

# %% 5. Gráfica y_i vs T
PNG_PATH = SALIDAS / "grafica_s26_ej1_reformado_gibbs.png"
if GUARDAR:
    from matplotlib.figure import Figure

    TINTA, TINTA_2, REJILLA, FONDO = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    COLORES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
    MARCAS = ["o", "s", "^", "D", "v"]
    fig = Figure(figsize=(7.5, 4.8), dpi=150, facecolor=FONDO)
    ax = fig.add_subplot(facecolor=FONDO)
    for sp, c, mk in zip(ESPECIES, COLORES, MARCAS):
        ys = [RESULTADOS[T].y[sp] for T in TEMPERATURAS]
        ax.plot(TEMPERATURAS, ys, lw=2, color=c, marker=mk, ms=6,
                markeredgecolor=FONDO, markeredgewidth=1, label=sp)
        ax.annotate(sp, (TEMPERATURAS[-1], ys[-1]), xytext=(6, 0),
                    textcoords="offset points", va="center", color=TINTA, fontsize=9)
    ax.set_yscale("log")
    ax.set_xlim(570, 1160)
    ax.set_xlabel("T (K)", color=TINTA)
    ax.set_ylabel("y_i (escala log)", color=TINTA)
    ax.set_title(f"Reformado de CH4, H2O:CH4 = 4:1, {P:g} bar (minimización de G)",
                 color=TINTA, fontsize=11, loc="left")
    ax.grid(True, which="major", color=REJILLA, lw=0.8)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color(TINTA_2)
    ax.tick_params(colors=TINTA_2)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13),
              labelcolor=TINTA, ncol=5, fontsize=9)
    fig.tight_layout()
    SALIDAS.mkdir(exist_ok=True)
    fig.savefig(PNG_PATH, facecolor=fig.get_facecolor())
    titulo("5. Salidas")
    print(f"Gráfica: {PNG_PATH.relative_to(SALIDAS.parent)}")

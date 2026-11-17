"""Sesión 28, Ejercicio — Si epitaxial a 1300 K y 1 bar por Lagrange con prueba de fase.

Mismo problema que ejercicios/s24_ej3_si_epitaxial.py (de ahí se importan los
datos y los resultados, no se copian), ahora como en la Sesión 28:
  1. Se resuelve solo el gas (SiCl2, SiCl4, SiCl3H, SiCl2H2, SiClH3, SiH4, H2,
     HCl) con las condiciones de Lagrange.
  2. Prueba de fase con los λ del gas: Δ_Si = g_Si(s)/RT + λ_Si/RT = λ_Si/RT
     (Si(s) de referencia, ΔfG = 0).
  3. Si Δ_Si < 0, Si(s) entra con g_Si(s)/RT + λ_Si/RT = 0 y n_Si en el
     balance de Si, y se recalcula (core.gibbs, `condensadas=["Si(s)"]`).
R mínimo para depositar desde SiCl4: Δ_Si(R) = 0 con el gas solo.

core usa R = 8.314 J/(mol·K); s24_ej3 usa 8.314462618 (CODATA), así que los
n_Si(s) difieren ~3e-5.
"""

import contextlib
import importlib.util
import io
import pathlib
import sys

_RAIZ = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_RAIZ))  # raíz del proyecto
import core.entorno  # noqa: F401

from scipy.optimize import brentq

from core.gibbs import minimizar_gibbs, prueba_de_fase, validar_gibbs


def importar_s24():
    """Carga ejercicios/s24_ej3_si_epitaxial.py sin imprimir ni escribir salidas."""
    spec = importlib.util.spec_from_file_location(
        "s24_ej3", _RAIZ / "ejercicios" / "s24_ej3_si_epitaxial.py")
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


S24 = importar_s24()

# ══════════════════════════ DATOS (de s24_ej3) ════════════════════════════
T, P = S24.T, S24.P                                  # 1300 K, 1 bar
GAS = list(S24.GAS)
DGF = {sp: 1000.0 * v for sp, v in S24.DGF.items() if sp != "Si"}   # J/mol
DGF["Si(s)"] = 1000.0 * S24.DGF["Si"]                # 0: Si de referencia (JANAF Si-001)
CLOROSILANOS = list(S24.CLOROSILANOS)
RELACIONES = list(S24.RELACIONES)


def titulo(t: str) -> None:
    print(f"\n{'=' * 84}\n{t}\n{'=' * 84}")


def equilibrio(clorosilano, R):
    return minimizar_gibbs(GAS, DGF, {clorosilano: 1.0, "H2": float(R)}, T, P,
                           condensadas=["Si(s)"])


def delta_gas(clorosilano, R):
    """Δ_Si de la prueba de fase con los λ del gas solo."""
    gas = minimizar_gibbs(GAS, DGF, {clorosilano: 1.0, "H2": float(R)}, T, P)
    return prueba_de_fase(gas, "Si(s)", DGF)


# %% 1. Barrido con prueba de fase en cada punto
RESULTADOS, VALIDACION = {}, {}
for cl in CLOROSILANOS:
    for R in RELACIONES:
        r = equilibrio(cl, R)
        with contextlib.redirect_stdout(io.StringIO()):
            VALIDACION[(cl, R)] = validar_gibbs(r, n_arranques=5)
        RESULTADOS[(cl, R)] = r

for cl in CLOROSILANOS:
    titulo(f"1. {cl} + R H2 a {T:g} K y {P:g} bar: prueba de fase y n_Si(s)")
    print(f"{'R':>5}{'Δ_Si gas solo':>15}{'decisión':>22}{'n_Si(s)':>11}{'s24_ej3':>11}"
          f"{'|Δn|':>10}{'λ_Si/RT':>10}{'λ_Cl/RT':>10}{'λ_H/RT':>9}  valid.")
    for R in RELACIONES:
        r = RESULTADOS[(cl, R)]
        n_si = r.n_cond.get("Si(s)", 0.0)
        s24 = S24.RESULTADOS[(cl, R)]["n_Si"]
        v = VALIDACION[(cl, R)]
        dec = "Si(s) entra" if "Si(s)" in r.n_cond else "solo gas (Δ > 0)"
        print(f"{R:>5g}{r.delta_fase_gas['Si(s)']:>15.4f}{dec:>22}{n_si:>11.5f}{s24:>11.5f}"
              f"{abs(n_si - s24):>10.1e}{r.lam_RT['Si']:>10.4f}{r.lam_RT['Cl']:>10.4f}"
              f"{r.lam_RT['H']:>9.4f}  {v.n_pass}/{v.n_warn}/{v.n_fail}")

# %% 2. R mínimo para depositar desde SiCl4
R_MIN = brentq(lambda R: delta_gas("SiCl4", R), 5.0, 10.0, xtol=1e-12)
titulo("2. R mínimo para depositar Si(s) desde SiCl4 (Δ_Si del gas solo = 0)")
print(f"R_min = {R_MIN:.4f}  (s24_ej3, con n_Si(s) = 0: {S24.R_MIN:.4f}; "
      f"diferencia {R_MIN - S24.R_MIN:+.1e}, por R = 8.314 vs 8.314462618)")
print("Desde SiCl3H, Δ_Si < 0 en todo el intervalo: siempre hay depósito.")

# %% 3. Contraste con s24_ej3
DIF = max(abs(RESULTADOS[k].n_cond.get("Si(s)", 0.0) - S24.RESULTADOS[k]["n_Si"])
          for k in RESULTADOS)
REGIMEN = all(("Si(s)" in RESULTADOS[k].n_cond) == (S24.RESULTADOS[k]["regimen"] == "Si(s)")
              for k in RESULTADOS)
titulo("3. Contraste con s24_ej3 (base SiCl4, SiCl3H, H2 y solver propio)")
print(f"máx |n_Si(s) Lagrange − n_Si(s) s24| = {DIF:.1e} en {len(RESULTADOS)} puntos")
print(f"mismo régimen (con o sin sólido) en todos los puntos: {'SÍ' if REGIMEN else 'NO'}")
print("La validación (PASS/WARN/FAIL) incluye el punto 9: Δ ≥ 0 sin sólido y n ≥ 0 con sólido.")

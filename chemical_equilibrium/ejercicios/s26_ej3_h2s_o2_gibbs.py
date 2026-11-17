"""Sesión 26, Ejercicio 3 — H2S y O2 en un reactor isotérmico (minimización de G).

Especies (todas gas): H2S, O2, H2O, SO2, S2, H2. ΔG°f de JANAF a 500 y 700 K
(data/janaf_dgf.json: H2S H-080, SO2 O-034, S2 S-012, H2O H-064).
  a) 2 mol H2S + 1 mol O2, 1 bar, 500 y 700 K
  b) 2:1, 500 K, 5 y 10 bar (gas ideal)
  c) 500 K, 1 bar, H2S:O2 = 1:1 y 1:2
  d) ¿Qué condición favorece el azufre elemental y cuál el SO2?

Aparte (no cambia la respuesta pedida): prueba de fase con azufre condensado.
En JANAF el S(ref) (S-001) tiene ΔfG = 0, así que Δ_S = λ_S/RT.
"""

import contextlib
import io
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto
import core.entorno  # noqa: F401

from core.gibbs import cargar_dgf, minimizar_gibbs, prueba_de_fase, validar_gibbs

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
ESPECIES = ["H2S", "O2", "H2O", "SO2", "S2", "H2"]
CASOS = {                      # inciso: (T K, P bar, mol H2S, mol O2)
    "a 500 K 1 bar 2:1":  (500.0, 1.0, 2.0, 1.0),
    "a 700 K 1 bar 2:1":  (700.0, 1.0, 2.0, 1.0),
    "b 500 K 5 bar 2:1":  (500.0, 5.0, 2.0, 1.0),
    "b 500 K 10 bar 2:1": (500.0, 10.0, 2.0, 1.0),
    "c 500 K 1 bar 1:1":  (500.0, 1.0, 1.0, 1.0),
    "c 500 K 1 bar 1:2":  (500.0, 1.0, 1.0, 2.0),
}


def titulo(t: str) -> None:
    print(f"\n{'=' * 78}\n{t}\n{'=' * 78}")


# %% 1. ΔG°f de JANAF
DGF, FUENTES = {}, {}
for T in sorted({c[0] for c in CASOS.values()}):
    DGF[T], FUENTES[T] = cargar_dgf(ESPECIES + ["S(ref)"], T)

titulo("1. ΔG°f (J/mol) de JANAF")
print(f"{'T (K)':>6}" + "".join(f"{sp:>10}" for sp in ESPECIES))
for T in DGF:
    print(f"{T:>6.0f}" + "".join(f"{DGF[T][sp]:>10.0f}" for sp in ESPECIES))
print("Tablas: " + ", ".join(f"{sp} {FUENTES[500.0][sp]}" for sp in ESPECIES))

# %% 2. Los 6 casos
RESULTADOS, VALIDACION = {}, {}
for caso, (T, P, nH2S, nO2) in CASOS.items():
    r = minimizar_gibbs(ESPECIES, DGF[T], {"H2S": nH2S, "O2": nO2}, T, P)
    with contextlib.redirect_stdout(io.StringIO()):
        VALIDACION[caso] = validar_gibbs(r)
    RESULTADOS[caso] = r

titulo("2. Composición de equilibrio (y_i) y moles totales")
print(f"{'caso':<20}" + "".join(f"{sp:>11}" for sp in ESPECIES) + f"{'n_T':>8}"
      + "  validación (PASS/WARN/FAIL)")
for caso, r in RESULTADOS.items():
    v = VALIDACION[caso]
    print(f"{caso:<20}" + "".join(f"{r.y[sp]:>11.4e}" for sp in ESPECIES)
          + f"{r.n_total:>8.4f}  {v.n_pass}/{v.n_warn}/{v.n_fail}")
if any(v.n_warn for v in VALIDACION.values()):
    print("WARN del punto 5 en 1:2: algunos arranques al azar no convergen (el Newton se")
    print("estanca cuando el O2, mayoritario, arranca como traza); la solución no cambia.")

titulo("3. Azufre: fracción de S alimentado que sale como S2 y como SO2")
FRAC = {}
for caso, r in RESULTADOS.items():
    S0 = r.n0["H2S"]
    FRAC[caso] = (2 * r.n["S2"] / S0, r.n["SO2"] / S0, r.n["H2S"] / S0)
    print(f"{caso:<20} S2: {FRAC[caso][0]:7.2%}   SO2: {FRAC[caso][1]:7.2%}   "
          f"H2S sin reaccionar: {FRAC[caso][2]:7.2%}")

# %% 4. d) Respuesta, leída de los resultados
a5, a7 = FRAC["a 500 K 1 bar 2:1"], FRAC["a 700 K 1 bar 2:1"]
b5, b10 = FRAC["b 500 K 5 bar 2:1"], FRAC["b 500 K 10 bar 2:1"]
c11, c12 = FRAC["c 500 K 1 bar 1:1"], FRAC["c 500 K 1 bar 1:2"]
titulo("4. d) ¿Qué favorece el azufre elemental y qué el SO2?")
print(f"""Azufre elemental (S2):
  * Temperatura alta: de 500 a 700 K el S que sale como S2 sube de {a5[0]:.1%} a {a7[0]:.1%}.
  * Presión baja: a 500 K baja de {a5[0]:.1%} (1 bar) a {b5[0]:.1%} (5 bar) y {b10[0]:.1%} (10 bar),
    porque el paso que forma S2, 2 H2S + SO2 -> 3/2 S2 + 2 H2O, aumenta los moles
    de gas (3 -> 3.5) y la presión lo desfavorece (Le Chatelier).
  * Relación H2S:O2 = 2:1 (la estequiométrica de 2 H2S + O2 -> S2 + 2 H2O): con 1:1 cae a {c11[0]:.1%}
    y con 1:2 desaparece ({c12[0]:.1e}).
SO2:
  * Exceso de O2: con 1:2 todo el azufre sale como SO2 ({c12[1]:.1%}); con 1:1, {c11[1]:.1%};
    con 2:1, {a5[1]:.1%}.
  * Temperatura baja: con 2:1 y 1 bar, {a5[1]:.1%} a 500 K contra {a7[1]:.1%} a 700 K.
  * Presión alta (efecto pequeño): {a5[1]:.1%} a 1 bar, {b5[1]:.1%} a 5 bar, {b10[1]:.1%} a 10 bar.""")

# %% 5. Aparte: prueba de fase con azufre condensado (no cambia la respuesta)
DELTA_S = {caso: prueba_de_fase(r, "S(ref)", DGF[r.T]) for caso, r in RESULTADOS.items()}
titulo("5. Aparte: prueba de fase con azufre condensado, Δ_S = λ_S/RT (S(ref): ΔfG = 0)")
print("Comentario aparte: los incisos a–d se resolvieron solo con gases, como pide el")
print("enunciado. Con los λ de ese gas, Δ_S < 0 dice que el azufre condensado bajaría G.")
for caso, r in RESULTADOS.items():
    extra = ""
    if DELTA_S[caso] < 0:
        con = minimizar_gibbs(ESPECIES, DGF[r.T], r.n0, r.T, r.P, condensadas=["S(ref)"])
        extra = (f" -> con S(ref): n_S = {con.n_cond.get('S(ref)', 0):.4f} mol, "
                 f"y_S2 = {con.y['S2']:.3e}")
    print(f"  {caso:<20} Δ_S = {DELTA_S[caso]:+8.3f}"
          + ("  (condensado favorable)" if DELTA_S[caso] < 0 else "  (no se forma)") + extra)

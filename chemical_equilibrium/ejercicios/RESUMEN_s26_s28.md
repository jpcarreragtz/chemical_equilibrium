# Resumen de las Sesiones 26, 27 y 28 (minimización de Gibbs, Fischer–Tropsch y Lagrange)

Generado el 2026-11-17 a partir de los ejercicios de `ejercicios/` (Python es la versión canónica; GAMS es la segunda vía). ΔG°f de JANAF consultados el 2026-11-17 y cotejados contra log Kf (|Δ| ≤ 0.018 kJ/mol). R = 8.314 J/(mol·K) en `core`.

Para correr: `.venv/bin/python ejercicios/<archivo>.py` desde `chemical_equilibrium/`; GAMS: `/Library/Frameworks/GAMS.framework/Versions/Current/Resources/gams <archivo>.gms` (opcional `--salida=archivo.csv`).

## S26 Ej1 — Reformado de CH4, H2O:CH4 = 4:1, 1 bar

| T (K) | y CH4 | y H2O | y CO | y CO2 | y H2 | n_T | a_C | Δ_C = λ_C/RT | vs Método A |
|---|---|---|---|---|---|---|---|---|---|
| 600 | 1.6951e-01 | 7.2168e-01 | 9.2173e-05 | 2.1689e-02 | 8.7031e-02 | 5.2277 | 0.2277 | +1.4799 | 1.6e-12 |
| 700 | 1.2379e-01 | 6.0602e-01 | 1.9773e-03 | 5.2455e-02 | 2.1575e-01 | 5.6108 | 0.3029 | +1.1943 | 6.9e-14 |
| 800 | 6.1639e-02 | 4.6032e-01 | 1.6108e-02 | 8.2721e-02 | 3.7921e-01 | 6.2318 | 0.3119 | +1.1651 | 7.9e-14 |
| 900 | 1.2034e-02 | 3.6586e-01 | 4.9203e-02 | 8.5059e-02 | 4.8784e-01 | 6.8355 | 0.1599 | +1.8331 | 9.8e-13 |
| 1000 | 9.2291e-04 | 3.5806e-01 | 6.9969e-02 | 7.2229e-02 | 4.9882e-01 | 6.9871 | 0.0387 | +3.2525 | 3.2e-14 |
| 1100 | 8.3675e-05 | 3.6785e-01 | 8.1922e-02 | 6.0876e-02 | 4.8927e-01 | 6.9988 | 0.0098 | +4.6269 | 3.5e-11 |

Coincide con los valores de referencia en todas las cifras dadas. Δ_C > 0 en todo el intervalo: no se forma grafito. Gráfica: `salidas/grafica_s26_ej1_reformado_gibbs.png`.

## S26 Ej2 — Steam cracker (ya resuelto en `s25_ej2_*`)

Confirmado: `s25_ej2_steam_cracker_gibbs.py` termina en *Veredicto: PASS (7 PASS)* y `s25_ej2_gibbs.gms` (IPOPT) da lo mismo. y_i: H2 0.604, H2O 0.1704, CO 0.1551, CO2 0.0629, CH4 0.007634, C2H6 1.695e-08, C2H4 1.041e-08, C2H2 3.972e-11, O2 6.027e-22.

## S26 Ej3 — H2S + O2 (reactor isotérmico, todo gas)

| caso | y H2S | y O2 | y H2O | y SO2 | y S2 | y H2 | n_T | Δ_S = λ_S/RT |
|---|---|---|---|---|---|---|---|---|
| a 500 K 1 bar 2:1 | 4.3801e-01 | 1.0936e-34 | 3.0166e-01 | 2.1901e-01 | 4.1326e-02 | 3.7598e-07 | 2.7039 | -4.303 |
| a 700 K 1 bar 2:1 | 3.1063e-01 | 2.8400e-24 | 4.0773e-01 | 1.5535e-01 | 1.2622e-01 | 6.3301e-05 | 2.7839 | -0.904 |
| b 500 K 5 bar 2:1 | 4.5744e-01 | 6.1643e-35 | 2.8546e-01 | 2.2872e-01 | 2.8372e-02 | 2.1193e-07 | 2.6921 | -4.920 |
| b 500 K 10 bar 2:1 | 4.6426e-01 | 4.8274e-35 | 2.7978e-01 | 2.3213e-01 | 2.3826e-02 | 1.6597e-07 | 2.6880 | -5.179 |
| c 500 K 1 bar 1:1 | 1.8268e-01 | 3.5599e-34 | 4.1485e-01 | 3.9010e-01 | 1.2373e-02 | 2.8658e-07 | 1.6736 | -3.700 |
| c 500 K 1 bar 1:2 | 1.3563e-50 | 2.0000e-01 | 4.0000e-01 | 4.0000e-01 | 4.1217e-68 | 1.1658e-23 | 2.5000 | +71.684 |

**d) Respuesta.** El **azufre elemental** (S2) se favorece con **T alta** (el S que sale
como S2 sube de 11.2% a 500 K a 35.1% a 700 K), **P baja** (11.2% a 1 bar,
7.6% a 5 bar, 6.4% a 10 bar; el paso que forma S2, 2 H2S + SO2 → 3/2 S2 + 2 H2O, pasa de 3 a 3.5 moles de gas, y la presión lo desfavorece) y
**H2S:O2 = 2:1** (con 1:1 cae a 4.1% y con 1:2 desaparece). El **SO2** se favorece con
**exceso de O2** (1:2: 100.0% del S como SO2; 1:1: 65.3%; 2:1: 29.6%), con
**T baja** (29.6% a 500 K contra 21.6% a 700 K) y, poco, con P alta.

*Aparte (comentario, no cambia la respuesta):* con los λ del gas, Δ_S < 0 en los casos 2:1
y 1:1, así que el azufre condensado S(ref) bajaría G: a 500 K y 2:1 se condensarían 1.96 de
los 2 mol de S. En 1:2, Δ_S = +71.7 (no se forma).

## S27 — Fischer–Tropsch a 500 K y 25 bar (datos de la maestra)

**SUPUESTO:** los Análisis 1 y 2 no dicen qué b usar; se toma H2/CO = 2 (C 1, H 4, O 1), en la variable `B_BASE` de `s27_fischer_tropsch.py`.

| alimentación | ASF | φ | n CH4 | n CO2 | n H2O | n C5–C15 | n H2 | % C1 | % C2–C4 | % C5–C11 | % C12–C15 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| base H2/CO=2 (supuesto) | no | no | 0.74542 | 0.25001 | 0.49996 | 3.043e-11 | 0.002364 | 99.43 | 0.57 | 0.00 | 0.00 |
| base H2/CO=2 (supuesto) | no | sí | 0.74508 | 0.24992 | 0.50015 | 4.937e-11 | 0.002208 | 99.38 | 0.62 | 0.00 | 0.00 |
| base H2/CO=2 (supuesto) | sí | no | 0.043117 | 0.073776 | 0.85244 | 0.08072 | 0.01335 | 5.15 | 26.37 | 53.18 | 15.29 |
| base H2/CO=2 (supuesto) | sí | sí | 0.043135 | 0.073379 | 0.85323 | 0.08075 | 0.01207 | 5.15 | 26.37 | 53.18 | 15.29 |
| biomasa H2/CO=1 | no | no | 0.47643 | 0.49658 | 0.0067587 | 1.283e-07 | 0.0001685 | 94.98 | 5.02 | 0.00 | 0.00 |
| biomasa H2/CO=1 | no | sí | 0.47504 | 0.49637 | 0.0071745 | 1.836e-07 | 0.0001619 | 94.68 | 5.32 | 0.00 | 0.00 |
| biomasa H2/CO=1 | sí | no | 0.028833 | 0.38059 | 0.23878 | 0.05398 | 0.002747 | 5.15 | 26.37 | 53.18 | 15.29 |
| biomasa H2/CO=1 | sí | sí | 0.028837 | 0.38051 | 0.23895 | 0.05398 | 0.00248 | 5.15 | 26.37 | 53.18 | 15.29 |
| gas natural H2/CO=2.5 | no | no | 0.8707 | 0.12549 | 0.74902 | 9.291e-12 | 0.003862 | 99.59 | 0.41 | 0.00 | 0.00 |
| gas natural H2/CO=2.5 | no | sí | 0.87046 | 0.12538 | 0.74924 | 1.508e-11 | 0.003606 | 99.55 | 0.45 | 0.00 | 0.00 |
| gas natural H2/CO=2.5 | sí | no | 0.046551 | 1.0171e-05 | 0.99998 | 0.08715 | 0.2755 | 5.15 | 26.37 | 53.18 | 15.29 |
| gas natural H2/CO=2.5 | sí | sí | 0.046551 | 7.302e-06 | 0.99999 | 0.08715 | 0.2755 | 5.15 | 26.37 | 53.18 | 15.29 |

Con ASF el % en peso es el analítico (C1 5.15, C2–C4 26.37, C5–C11 53.18, C12–C15 15.29) y no depende de la alimentación ni de φ.

**Pregunta A.** Sin ASF el equilibrio es casi puro metano (y_C1 ≈ 0.50, 99.4 % en peso
de las parafinas) y C5–C15 queda en ~1e-11 mol: cada CH2 adicional sube ΔG y la termodinámica
manda todo el carbono a C1, CO2 y H2O. Con ASF, la distribución la fija la cinética de
crecimiento de cadena (α), no el equilibrio. Un reactor de equilibrio puro no sirve para
simular combustibles líquidos en FT porque predice metano; el producto real es la
distribución ASF que el catalizador produce lejos del equilibrio entre parafinas. Hay que
imponer la ASF (u otro modelo cinético) como restricción.

## S28 — Reformado por Lagrange (Koretsky 9.45–9.47, sistema explícito)

**4:1 (enunciado)** (λ_j/RT y prueba de fase de C(s)):

| T (K) | λ_C/RT | λ_H/RT | λ_O/RT | Δ_C | C(s) | vs RAND |
|---|---|---|---|---|---|---|
| 600 | 1.4799 | 1.2207 | 40.7856 | +1.4799 | no se forma (Δ_C ≥ 0) | 1.6e-12 |
| 700 | 1.1943 | 0.7668 | 34.8468 | +1.1943 | no se forma (Δ_C ≥ 0) | 6.6e-14 |
| 800 | 1.1651 | 0.4848 | 30.4015 | +1.1651 | no se forma (Δ_C ≥ 0) | 8.0e-14 |
| 900 | 1.8331 | 0.3589 | 26.7602 | +1.8331 | no se forma (Δ_C ≥ 0) | 9.8e-13 |
| 1000 | 3.2525 | 0.3478 | 23.4961 | +3.2525 | no se forma (Δ_C ≥ 0) | 2.4e-14 |
| 1100 | 4.6269 | 0.3574 | 20.7363 | +4.6269 | no se forma (Δ_C ≥ 0) | 3.5e-11 |

**2:3 (balances de la diapositiva)** (λ_j/RT y prueba de fase de C(s)):

| T (K) | λ_C/RT | λ_H/RT | λ_O/RT | Δ_C | C(s) | vs RAND |
|---|---|---|---|---|---|---|
| 600 | 0.7798 | 1.2077 | 41.1232 | +0.7798 | no se forma (Δ_C ≥ 0) | 1.4e-13 |
| 700 | 0.3827 | 0.7499 | 35.2424 | +0.3827 | no se forma (Δ_C ≥ 0) | 1.4e-14 |
| 800 | 0.0860 | 0.4531 | 30.9543 | +0.0860 | no se forma (Δ_C ≥ 0) | 1.9e-13 |
| 900 | 0.1014 | 0.2764 | 27.7502 | +0.1014 | no se forma (Δ_C ≥ 0) | 1.8e-14 |
| 1000 | 0.7586 | 0.2002 | 25.0785 | +0.7586 | no se forma (Δ_C ≥ 0) | 2.3e-14 |
| 1100 | 1.9941 | 0.1871 | 22.4971 | +1.9941 | no se forma (Δ_C ≥ 0) | 2.0e-13 |

Demostración fuera del enunciado: con H2O:CH4 = 1:1 a 800 K el gas da Δ_C = −0.314 y C(s) entra con n_C = 0.14511 mol (igual con `minimizar_gibbs`).

## S28 — Si epitaxial a 1300 K y 1 bar por Lagrange con prueba de fase

| R | SiCl4: Δ_Si gas | SiCl4: n_Si(s) | SiCl3H: Δ_Si gas | SiCl3H: n_Si(s) | s24 SiCl4 | s24 SiCl3H |
|---|---|---|---|---|---|---|
| 1 | +1.4424 | 0.00000 | -3.7267 | 0.21310 | 0.00000 | 0.21309 |
| 2 | +1.0143 | 0.00000 | -3.4383 | 0.21928 | 0.00000 | 0.21927 |
| 3 | +0.7399 | 0.00000 | -3.2672 | 0.22707 | 0.00000 | 0.22707 |
| 5 | +0.3738 | 0.00000 | -3.0794 | 0.24412 | 0.00000 | 0.24412 |
| 7 | +0.1231 | 0.00000 | -2.9869 | 0.26153 | 0.00000 | 0.26153 |
| 10 | -0.1492 | 0.01537 | -2.9253 | 0.28708 | 0.01537 | 0.28708 |
| 15 | -0.4658 | 0.05774 | -2.9084 | 0.32721 | 0.05774 | 0.32721 |
| 20 | -0.6948 | 0.09779 | -2.9330 | 0.36415 | 0.09780 | 0.36416 |
| 30 | -1.0241 | 0.17080 | -3.0197 | 0.42957 | 0.17081 | 0.42959 |
| 50 | -1.4516 | 0.29282 | -3.2103 | 0.53372 | 0.29284 | 0.53374 |
| 75 | -1.8015 | 0.41128 | -3.4164 | 0.62828 | 0.41130 | 0.62830 |
| 100 | -2.0549 | 0.50225 | -3.5847 | 0.69609 | 0.50228 | 0.69611 |
| 125 | -2.2537 | 0.57314 | -3.7238 | 0.74586 | 0.57317 | 0.74588 |
| 150 | -2.4170 | 0.62911 | -3.8407 | 0.78323 | 0.62914 | 0.78325 |

R mínimo (SiCl4) = 8.2328 (s24_ej3: 8.2327). Máx |n_Si Lagrange − s24| = 2.8e-05 (R = 8.314 contra 8.314462618).

## GAMS (segunda vía)

| modelo | casos | NLP | CNS | máx rel. y_i > 1e-6 | trazas |
|---|---|---|---|---|---|
| `s26_ej1_reformado.gms` | 6 T | 2 (Locally Optimal) | 16 (Solved) | 3.5e-11 | 0 |
| `s26_ej3_h2s.gms` | 6 casos | 2 | 16 | 1.7e-9 | 1.6e-9 |
| `s27_fischer_tropsch.gms` | 12 escenarios | 2 | 16 | 5.4e-9 | 3.1e-8 |
| `s28_si_lagrange.gms` | 28 puntos | 2 (3 con reintento) | 16 gas / 16 con Si(s) | 8.2e-10 | 1.1e-9 |

En `s28_si_lagrange.gms`, el NLP desde el reparto uniforme no llega en 3 puntos (SiCl4 R = 125
y SiCl3H R = 75: "local infeasibility"; SiCl3H R = 125: "search direction too small"). Se
reintenta desde el gas del R anterior y queda Locally Optimal; el CNS da Solved en todos.

## Supuestos

* S27, Análisis 1 y 2: b = C 1, H 4, O 1 (H2/CO = 2), en `B_BASE`.
* S27: alimentación como CO + H2 a partir de b (exige C = O, que se cumple en las tres).
* S27 %: PM = 12k + (2k + 2), como en la diapositiva (con 12.011 y 1.008, C12–C15 = 15.287).
* S28 reformado: la prueba de fase usa grafito de referencia (JANAF C-002, ΔfG = 0).

## Inconsistencias detectadas

1. **S28, diapositiva del ejemplo:** dice "4:1" pero escribe los balances C = 2, H = 14,
   O = 3, que son los de 2 CH4 + 3 H2O (SVA Ej. 13.14). Con 4:1 son C = 1, H = 12, O = 4. Se
   resolvieron ambos.
   **Resuelto (17-nov-2026):** la maestra confirmó 4:1 (C = 1, H = 12, O = 4); ese es el caso que se entrega. El 2:3 queda solo como referencia.
2. **S27, ΔG a 500 K y φ:** Datos de la maestra (Sesión 27), confirmados por ella el
   17-nov-2026; se usan tal cual.

## Lo que no se pudo verificar

* La plantilla de GAMS de la maestra (con y.lo = 1e-8) no estaba disponible: la nota de
  CLAUDE.md sobre las trazas en 1e-8 se deduce de la cota, no se corrió esa plantilla.
* Las ecuaciones 9.45–9.47 de Koretsky se tomaron del enunciado, no del libro.


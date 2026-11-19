# Ejercicios — índice

Un folder por ejercicio, `sNN_ejN_tema/`, con su `.py` y, cuando lo hay, su `.gms`.
Los resultados clave salen de correr cada script; el detalle de S26–S28 está en
[RESUMEN_s26_s28.md](RESUMEN_s26_s28.md).

| folder | sesión | tema | método | resultado clave |
|---|---|---|---|---|
| [s06_ej2_reformado_propano/](s06_ej2_reformado_propano/) | 6 | Reformado de propano con vapor, 3 reacciones (avance ξ), 700–1000 K y 1 bar | Python | X_C3H8 = 1.000 en las tres T; n_H2 = 0.855, 3.366 y 5.232 mol a 700, 900 y 1000 K |
| [s09_ej3_butano/](s09_ej3_butano/) | 9 | Craqueo de n-butano a propileno y metano, 500 K (Koretsky) | Python | K(500 K) = 1.129; ξ = 7.283 mol/s a 1 bar y 2.079 mol/s a 25 bar |
| [s15_etileno_etanol/](s15_etileno_etanol/) | 15 | Hidratación de etileno con fases L + V, 200 °C y 34.5 bar (SVA Ej. 13.13) | Python | K = 0.03104; x_EtOH = 0.0798; y_C2H4 = 0.350, y_H2O = 0.463, y_EtOH = 0.187 |
| [s24_ej1_reformado_matrices/](s24_ej1_reformado_matrices/) | 23–24 | Reacciones independientes por matrices: reformado de metano | Python | r = 2: CH4 + H2O → CO + 3 H2 y CH4 + 2 H2O → CO2 + 4 H2; F = 2 |
| [s24_ej2_oxidacion_h2s_s8/](s24_ej2_oxidacion_h2s_s8/) | 23–24 | Reacciones independientes con sólido puro: H2S + O2 con S8(s) | Python | r = 3 (SO2, SO3 y H2O desde H2S, O2 y S8); π = 2, F = 2 |
| [s24_ej3_si_epitaxial/](s24_ej3_si_epitaxial/) | 23–24 | Si(s) depositado desde SiCl4 o SiCl3H con H2, 1300 K y 1 bar | Python | R mínimo (SiCl4) = 8.2327; n_Si(s) a R = 150: 0.629 (SiCl4) y 0.783 (SiCl3H) |
| [s24_ej4_perclorato/](s24_ej4_perclorato/) | 23–24 | Reacciones independientes: descomposición del perclorato de amonio | Python | 9 especies, rango 4, r = 5 reacciones con las componentes del enunciado |
| [s25_ej1_gibbs_reformado/](s25_ej1_gibbs_reformado/) | 25 | Minimización de Gibbs: 2 CH4 + 3 H2O a 1000 K y 1 bar (SVA Ej. 13.14) | Python | y_H2 = 0.670, y_CO = 0.174, y_H2O = 0.0988, y_CO2 = 0.0369, y_CH4 = 0.0202 (JANAF) |
| [s25_ej2_steam_cracker/](s25_ej2_steam_cracker/) | 25 (y S26 Ej2) | Steam cracker: 1 C2H6 + 4 H2O, 9 especies, 1000 K y 1 bar | Python y GAMS | y_H2 = 0.604, y_H2O = 0.170, y_CO = 0.155, y_CO2 = 0.0629, y_CH4 = 0.00763 |
| [s26_ej1_reformado_barrido/](s26_ej1_reformado_barrido/) | 26 | Reformado H2O:CH4 = 4:1 a 1 bar, barrido 600–1100 K, prueba de fase de C(s) | Python y GAMS | y_H2 sube de 0.087 (600 K) a 0.499 (1000 K); y_CH4 = 9.229e-04 a 1000 K; no se forma C(s) |
| [s26_ej3_h2s_o2/](s26_ej3_h2s_o2/) | 26 | H2S + O2 en reactor isotérmico, 6 casos de T, P y relación | Python y GAMS | S que sale como S2: 11.2 % (500 K, 1 bar, 2:1), 35.1 % a 700 K y 0 % con 1:2 (todo SO2) |
| [s27_fischer_tropsch/](s27_fischer_tropsch/) | 27 | Fischer–Tropsch a 500 K y 25 bar, con y sin ASF (α = 0.80) y φ | Python y GAMS | Base H2/CO = 2: sin ASF n_CH4 = 0.745 mol y C5–C15 = 3.0e-11 mol; con ASF, 0.0431 y 0.0807 mol |
| [s28_ej1_lagrange_reformado/](s28_ej1_lagrange_reformado/) | 28 | Reformado como sistema de Lagrange explícito (Koretsky 9.45–9.47), 4:1 y 2:3 | Python | Mismos y_i que S26 Ej1 (y_CH4 = 9.229e-04 a 1000 K con 4:1); con 1:1 a 800 K entra C(s): n_C = 0.145 mol |
| [s28_ej2_si_epitaxial/](s28_ej2_si_epitaxial/) | 28 | Si epitaxial por Lagrange con prueba de fase, 1300 K y 1 bar | Python y GAMS | R mínimo (SiCl4) = 8.2328; n_Si(s) a R = 150: 0.629 (SiCl4) y 0.783 (SiCl3H) |

`s15_etileno_etanol/entrega/` guarda el entregable de la sesión 15 (resumen y CSV de
sensibilidad); el script lo vuelve a escribir igual cada vez que corre.

## Cómo correr un ejercicio

**Python.** Con el botón ▶ ("Run Python File") de VS Code sobre el `.py`, o desde la
terminal, parado en `chemical_equilibrium/`:

```bash
python ejercicios/<folder>/<archivo>.py
```

También corre desde su propio folder (`cd ejercicios/<folder> && python <archivo>.py`):
cada script busca hacia arriba la carpeta que contiene `core/`, así que no depende de
dónde esté parado ni de qué tan profundo esté el archivo.

**GAMS.** Con el `.gms` abierto en VS Code, `Cmd+Shift+P` → "GAMS: Run this file". O desde
la terminal, dentro del folder del ejercicio:

```bash
cd ejercicios/<folder>
gams <archivo>.gms          # con lo=0 la terminal muestra solo la tabla
```

La tabla de resultados sale al final de la terminal y se guarda en
`<archivo>_resultados.txt`, junto al `.gms`. Ese archivo, el `.lst` y `ipopt.opt` se
regeneran en cada corrida y no se versionan.

## Ejercicio nuevo

Crea `ejercicios/sNN_ejN_tema/`, copia ahí `plantilla/plantilla_ejercicio.py` y, si
lleva GAMS, `plantilla/plantilla_gibbs.gms`, y agrega una fila a la tabla de arriba.

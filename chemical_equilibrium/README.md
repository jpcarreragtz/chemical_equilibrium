# chemical_equilibrium — mapa del proyecto

Solver de equilibrio químico del curso (Método A: avance de reacción ξ,
Smith–Van Ness; Método B: tabla estequiométrica de Fogler) con
autovalidación, más los ejercicios y exámenes resueltos con él.

## Árbol

```
chemical_equilibrium/
├── README.md                    # este mapa
├── core/                        # MOTOR — no se toca salvo para corregir bugs
│   ├── __init__.py              # rutas del proyecto: ROOT, DATA, SALIDAS
│   ├── entorno.py               # si a python3 le falta numpy/scipy/matplotlib, relanza el script con .venv/bin/python
│   ├── equilibrium.py           # solvers (Método A y Método B), K<->Kc, checks
│   ├── reactions.py             # parser: "C3H8 + 3 H2O = 3 CO + 7 H2" -> nu
│   ├── thermo.py                # Cp/R (C.1), ec. 13.18, Van't Hoff, carga del JSON
│   ├── nonideal.py              # phi de Pitzer (virial), gammas de Wilson y Margules, K_phi (Lewis-Randall)
│   ├── convenciones.py          # ÚNICO lugar que decide cómo entra K_phi: "maestra" (default) | "rigurosa"
│   ├── energia.py               # 2do parcial: balance de energía con reacción, HHV/LHV, T flama adiabática
│   ├── exercise.py              # motor de plantilla: run_exercise(...)
│   ├── gibbs.py                 # S25–S28: minimización de G sin reacciones ni valor inicial (RAND), restricciones extra (ASF), fases condensadas con prueba de fase, datos y validación
│   └── selfcheck.py             # 12 puntos de autovalidación (validate, sweep)
├── data/
│   ├── sva_tables.json          # Smith–Van Ness C.1/C.4, 42 especies verificadas
│   ├── criticas_b1.json         # Tc, Pc, ω (Tabla B.1) para phi de Pitzer
│   ├── janaf_dgf.json           # ΔG°f por T capturados de janaf.nist.gov, con su tabla (gases y fases de referencia)
│   ├── ft_s27_maestra.json      # S27 Fischer–Tropsch: ΔG a 500 K y φ de la maestra (Sesión 27, diapositiva Datos)
│   └── raw/                     # (fotos/pdf de tablas originales; fuera de git)
├── tests/
│   ├── tests.py                 # casos a–x con respuesta conocida (exit 0 = PASS)
│   └── main.py                  # suite original de 4 casos del curso
├── plantilla/
│   ├── plantilla_ejercicio.py   # copia esto para cada ejercicio nuevo (equilibrio; fugacidad="lewis" opcional)
│   ├── plantilla_presion_fugacidad.py  # 1 reacción gas: y a (T, P) o T para y objetivo, gi | Lewis-Randall
│   └── plantilla_balance_energia.py  # Q / T flama adiabática / HHV-LHV (2do parcial)
├── ejercicios/                  # tareas/clase: ej2, sesión 9, sesión 15, ...
│   ├── s26_ej1_reformado_gibbs_barrido.py  # S26: reformado 4:1, 600–1100 K, Método A, prueba de fase C(s), gráfica
│   ├── s26_ej3_h2s_o2_gibbs.py   # S26: H2S + O2, 6 casos a–c, respuesta d) y prueba de fase con S(ref)
│   ├── s27_fischer_tropsch.py    # S27: FT con/sin ASF, con/sin φ, 3 alimentaciones, % en peso, Pregunta A
│   ├── s28_lagrange_reformado.py # S28: sistema de Lagrange explícito (Koretsky 9.45–9.47), 4:1 y 2:3, C(s)
│   ├── s28_si_epitaxial_lagrange.py  # S28: Si epitaxial por Lagrange con prueba de fase vs s24_ej3
│   ├── s26_ej1_reformado.gms, s26_ej3_h2s.gms, s27_fischer_tropsch.gms, s28_si_lagrange.gms  # GAMS (NLP + CNS, IPOPT)
│   └── RESUMEN_s26_s28.md        # resultados, respuestas, supuestos e inconsistencias de S26–S28
├── notebooks/                   # estudio (00_teoria) y versión-reporte (ej2_reporte); NO canónicos
│   └── entrega_s15/             # entregable de la sesión 15
├── examenes/
│   ├── parcial1_otono2026/      # parcial1_p1/p2/p3.py (ejecutables)
│   │   ├── entrega_parcial1/    # CONGELADA: evidencia entregada el 15/09/2026
│   │   └── preliminares/        # versiones previas del mismo parcial (históricas)
│   ├── parcial2_otono2026/practica/  # 5 ejercicios adicionales del 2do parcial + GUIA_planteamientos.md
│   ├── parcial2_otono2026/simulacro/ # simulacro: 1 problema por tipo, copias de las plantillas (correr_simulacro.py)
│   └── primavera2026/           # o-cresol, n-butano→anhídrido maleico, pirólisis
└── salidas/                     # CSV/PNG regenerables (default de run_exercise)
```

## Cómo hacer un ejercicio nuevo

```bash
cp plantilla/plantilla_ejercicio.py ejercicios/ejN_tema.py
# edita SOLO el bloque INPUTS (reacciones como texto, alimentación, fase,
# P, T, K directas o especies_termo para calcular K(T) del JSON)
cd ejercicios && python3 ejN_tema.py
```

Lee el reporte de AUTOVALIDACIÓN: `[FAIL]` = no confíes, `[WARN]` = revisa
a mano, `[info]` 4/5/12 = compara δ_j y n_T(ξ), la expresión simbólica de
cada K y el factor (P/P°)^δ contra tu derivación a mano. CSV y gráfica
quedan en `salidas/`.

## Ejecutar por celdas en VS Code (tablas y gráficas inline)

`plantilla/plantilla_ejercicio.py` y `ejercicios/ej2_mio.py` están divididos
en celdas con marcadores `# %%` (comentarios: no cambian nada al correr
como script). Para verlos en el Interactive Window:

1. Selecciona el intérprete correcto: `Cmd+Shift+P` → "Python: Select
   Interpreter" (el que tenga numpy/scipy/matplotlib; para las tablas
   también `pandas`: `pip3 install pandas`).
2. Abre el `.py`, pon el cursor en una celda y `Shift+Enter` (o el botón
   "Run Cell" que aparece sobre `# %%`). "Run All" corre el archivo de
   arriba a abajo.
3. El Interactive Window muestra la salida de texto, los DataFrames
   (K(T), ξ/n_i/y_i, resumen de selfcheck) y la gráfica y_i vs T inline.
   Las celdas de tabla/gráfica solo tienen efecto en el kernel.

`python3 archivo.py` desde la carpeta del script sigue siendo la forma
canónica y la que valida selfcheck: su salida es exactamente la misma que
antes de dividir el archivo en celdas. Si al correr una celda falla el
`sys.path`, revisa que el directorio de trabajo del kernel sea la carpeta
del script (ajuste `jupyter.notebookFileRoot`, por omisión `${fileDirname}`).

## Notebooks (`notebooks/`)

Los `.py` de `ejercicios/` son la versión canónica y validada; los notebooks
son para estudiar y para versión-reporte. `00_teoria.ipynb` explica ec.
13.18, Cp(T), Van't Hoff y el sistema en ξ con ejemplos calculados por
`core.thermo` (K_WGS(1100 K) ≈ 1.006, K_NH3(298 K) ≈ 762, K_NH3(500 K) ≈
0.316); `ej2_reporte.ipynb` importa el solver desde `core`, agrega la
derivación, tablas y gráfica, y termina en selfcheck (mismos números que
`ejercicios/ej2_mio.py`). Se regeneran con
`cd notebooks && jupyter nbconvert --to notebook --execute --inplace *.ipynb`.

## Cómo correr un ejercicio

**Terminal.** Desde cualquier carpeta y con cualquier `python3`:

```bash
cd "/Users/jpcg/Library/Mobile Documents/com~apple~CloudDocs/IBERO/EQ. Químico/chemical_equilibrium"
python3 ejercicios/s27_fischer_tropsch.py      # o: cd ejercicios && python3 s27_fischer_tropsch.py
```

La ruta va **entre comillas** porque "Mobile Documents" y "EQ. Químico"
llevan espacio; sin comillas, `cd` falla con "too many arguments" o "no such
file or directory".

Cada script de `ejercicios/`, `plantilla/` y `examenes/` importa
`core.entorno` justo después de poner la raíz en `sys.path`. Si el `python3`
de la terminal no tiene numpy, scipy o matplotlib (por ejemplo
`/usr/bin/python3`, el 3.9 de Apple, que no tiene scipy y daba
`ModuleNotFoundError: No module named 'scipy'`), el script se relanza solo
con `.venv/bin/python`. Si no existe el `.venv`, imprime los comandos para
crearlo. Para usar el entorno del proyecto de forma explícita:

```bash
.venv/bin/python ejercicios/s27_fischer_tropsch.py
```

**VS Code.** Abre la carpeta `chemical_equilibrium/` (no `EQ. Químico/`: el
`.venv` quedaría en una subcarpeta y VS Code puede no detectarlo).
`.vscode/settings.json` (local, fuera de git) fija
`python.defaultInterpreterPath` en `${workspaceFolder}/.venv/bin/python`, así
que "Run Python File" usa el `.venv`. Si la barra de estado muestra otro
intérprete: `Python: Select Interpreter` → `.venv`. Para correr por celdas
`# %%`, elige también el kernel `.venv`; un kernel no se puede relanzar solo
y `core.entorno` solo avisa qué falta.

**GAMS.** `gams ejercicios/s27_fischer_tropsch.gms`; si `gams` no está en el
PATH, `/Library/Frameworks/GAMS.framework/Versions/Current/Resources/gams`.
Cada `.gms` de S25–S28 imprime al final una tabla de resultados con los títulos
del `.py` y una línea de estado de los solves, y la guarda en
`ejercicios/<nombre>_resultados.txt`. El `.lst` solo trae los `display` finales
e IPOPT no escribe en la terminal. Con `lo=0` la terminal muestra solo la tabla.

## Cómo correr los tests

```bash
python3 tests/tests.py     # o: cd tests && python3 tests.py
```

Ambas formas funcionan (cada script ejecutable inserta la raíz del
proyecto en `sys.path` con un bloque `pathlib ... parents[N]`). Exit 0 =
todos los casos a–n PASS. `tests/main.py` corre la suite original (4/4).

## Convención de K_phi (fugacidad en fase gas)

`core/convenciones.py` es el ÚNICO lugar donde se decide de qué lado entra
K_phi = Π phi_i^ν (phi_i de especie pura a T, P; Lewis-Randall, Pitzer):

| `CONVENCION_KPHI` | fórmula | Problema 1 d) |
|---|---|---|
| `"maestra"` (**default**, la que pide la profesora) | K = K_y·(P/P°)^δ / K_phi  ⇔  K_y = K·(P/P°)^(−δ)·K_phi | 570.4 K |
| `"rigurosa"` (f̂_i = φ_i y_i P, SVA ec. 13.28) | K = K_phi·K_y·(P/P°)^δ  ⇔  K_y = K·(P/P°)^(−δ)/K_phi | 584.5 K |

Ambas se escriben K = F·K_y·(P/P°)^δ con F = `factor_Kphi(K_phi, convencion)`
y el solver ideal recibe `K_efectiva(K, K_phi, convencion)` = K/F.
`run_exercise(..., fugacidad="lewis", convencion_Kphi=None)` y
`plantilla_presion_fugacidad.py` la usan; cualquier script puede forzar la
otra con `convencion_Kphi="rigurosa"`. Para cambiar el default del proyecto
edita solo `CONVENCION_KPHI`.

## Convenciones

- Los scripts se corren `python3 archivo.py` DESDE la carpeta del script
  (o desde donde sea: no dependen del cwd).
- Imports del motor siempre como `from core.thermo import ...`; dentro de
  `core/` son relativos (`from .equilibrium import ...`).
- Datos: `core` expone `ROOT`, `DATA`, `SALIDAS`; el JSON se carga con
  `core.thermo.cargar_datos([...])` (ruta `DATA / "sva_tables.json"`).
- Salidas regenerables → `salidas/` (default de `run_exercise`) o junto al
  script del examen; NUNCA en la raíz.
- `examenes/.../entrega_parcial1/` está congelada como evidencia (ver su
  LEEME.txt): sus copias de scripts apuntan a la estructura vieja a
  propósito; las versiones ejecutables viven un nivel arriba.

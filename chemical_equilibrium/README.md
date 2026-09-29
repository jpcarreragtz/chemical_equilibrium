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
│   ├── equilibrium.py           # solvers (Método A y Método B), K<->Kc, checks
│   ├── reactions.py             # parser: "C3H8 + 3 H2O = 3 CO + 7 H2" -> nu
│   ├── thermo.py                # Cp/R (C.1), ec. 13.18, Van't Hoff, carga del JSON
│   ├── nonideal.py              # phi de Pitzer (virial), gammas de Wilson
│   ├── exercise.py              # motor de plantilla: run_exercise(...)
│   └── selfcheck.py             # 12 puntos de autovalidación (validate, sweep)
├── data/
│   ├── sva_tables.json          # Smith–Van Ness C.1/C.4, 42 especies verificadas
│   └── raw/                     # (fotos/pdf de tablas originales; fuera de git)
├── tests/
│   ├── tests.py                 # casos a–l con respuesta conocida (exit 0 = PASS)
│   └── main.py                  # suite original de 4 casos del curso
├── plantilla/
│   └── plantilla_ejercicio.py   # copia esto para cada ejercicio nuevo
├── ejercicios/                  # tareas/clase: ej2, sesión 9, sesión 15, ...
│   └── entrega_s15/             # entregable de la sesión 15
├── examenes/
│   ├── parcial1_otono2026/      # parcial1_p1/p2/p3.py (ejecutables)
│   │   ├── entrega_parcial1/    # CONGELADA: evidencia entregada el 15/09/2026
│   │   └── preliminares/        # versiones previas del mismo parcial (históricas)
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

## Cómo correr los tests

```bash
python3 tests/tests.py     # o: cd tests && python3 tests.py
```

Ambas formas funcionan (cada script ejecutable inserta la raíz del
proyecto en `sys.path` con un bloque `pathlib ... parents[N]`). Exit 0 =
todos los casos a–l PASS. `tests/main.py` corre la suite original (4/4).

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

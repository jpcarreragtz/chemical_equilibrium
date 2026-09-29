# Equilibrio Químico — IBERO CDMX (Otoño 2026)

Trabajo del curso de **Equilibrio Químico** (Ingeniería Química, Universidad
Iberoamericana, Grupo B, Otoño 2026): un solver propio en Python con
autovalidación, la base de datos de tablas de Smith–Van Ness transcrita y
verificada, y los ejercicios/parciales resueltos con ambos.

## Estructura

```
chemical_equilibrium/   Solver en Python (Método A: avance de reacción ξ,
                        Smith–Van Ness; Método B: tabla estequiométrica de
                        Fogler) + ejercicios y parciales resueltos
tablas_sva/             Base de datos Smith–Van Ness 7a ed. (Apéndices B y C)
                        en Excel y JSON, verificada con doble pase
excel/                  Libros de Excel de ejercicios resueltos (uno por
                        ejercicio)
```

Dentro de `chemical_equilibrium/` (ver su README para el detalle):

- `core/`: `equilibrium.py` (los dos solvers), `reactions.py` (parser de
  reacciones en texto), `thermo.py` (Cp/R, ΔH°(T), ΔG°(T), ec. 13.18, Van't
  Hoff), `nonideal.py` (phi de Pitzer, Wilson), `exercise.py` (motor de
  plantilla) y `selfcheck.py` (12 puntos de autovalidación);
  `tests/tests.py` (casos de respuesta conocida).
- `plantilla/plantilla_ejercicio.py` — se copia para cada ejercicio nuevo
  (`modo_K = "directo" | "termo"`).
- `data/sva_tables.json` — 42 especies gaseosas (C.1 + C.4) que carga el
  modo termo, verificadas contra el libro.
- `ejercicios/` y `examenes/`: `ej2_mio.py`, `sesion9_ej3_butano.py`,
  `parcial1_otono2026/parcial1_p1/p2/p3.py`, etc.; la entrega del primer
  parcial (congelada) en `examenes/parcial1_otono2026/entrega_parcial1/`.

## Cómo correr

```bash
pip3 install -r requirements.txt
cd chemical_equilibrium
python3 tests/tests.py      # casos a–l de respuesta conocida; exit 0 = todo PASS
cd examenes/parcial1_otono2026 && python3 parcial1_p3.py   # cualquier script igual
```

Cada script se corre `python3 archivo.py` desde su propia carpeta (o desde
donde sea: insertan la raíz del proyecto en `sys.path` y
`data/sva_tables.json` se resuelve vía `core/__init__.py`).

## Flujo de trabajo para un ejercicio nuevo

1. **Resolver a mano** (o al menos plantear: reacciones balanceadas, balances
   n_i(ξ), expresión de cada K).
2. Copiar la plantilla: `cp plantilla/plantilla_ejercicio.py
   ejercicios/ejN_nombre.py` y llenar
   SOLO el bloque INPUTS (reacciones como texto, alimentación, fase, P, T,
   K directas o `especies_termo` para calcular K(T) de las tablas).
3. Correr y **verificar las expresiones simbólicas** que imprime el punto 5
   del selfcheck (cada K en función de y_i y (P/P°)^δ) contra lo planteado a
   mano, igual que δ_j y n_T(ξ) del punto 4.
4. Revisar el **selfcheck**: `[FAIL]` = no confiar; `[WARN]` = revisar a
   mano; los 12 puntos cubren balance de elementos, residuos de K,
   re-solución desde 20 arranques aleatorios, contraste Método A vs B y
   Le Chatelier en barridos de T.
5. Interpretar: tabla comparativa, X del limitante, CSV y gráfica y_i vs T.

## Notas

- Los PDF del libro/apéndices y las fotos de tablas (`data/raw/`) están
  **excluidos** del repo por derechos de autor (ver `.gitignore`); las tablas
  usadas viven transcritas y verificadas en `tablas_sva/`.
- Suite original de referencia: `python3 tests/main.py` (4 tests del curso).

# Equilibrio químico — solver del curso (SVA / Fogler)

## Cómo meter un problema nuevo en 2 minutos

```bash
cp plantilla_ejercicio.py ej5_mi_problema.py   # 1. copia la plantilla
# 2. edita SOLO el bloque INPUTS (abajo el ejemplo)
python3 ej5_mi_problema.py                      # 3. corre
```

Obtienes por cada T: ξ_j, n_i, y_i, X del reactivo limitante, K usadas
(+ ΔCp/R, ΔH°(T), ΔG°(T) si es modo termo), **autovalidación automática**
(balance de elementos, n_i ≥ 0, |K_calc−K|/K < 1e-4, re-solución desde 20
arranques aleatorios, Le Chatelier del barrido), tabla comparativa, **CSV**
y **gráfica y_i vs T** (PNG).

Ejemplo llenado (el caso del examen — pirólisis de propano):

```python
titulo = "Pirólisis de propano (caso examen)"
reacciones = ["C3H8 = C2H4 + CH4",     # escribe TÚ los coeficientes:
              "C3H8 = C3H6 + H2"]      # "C3H8 + 3 H2O = 3 CO + 7 H2"
alimentacion = {"C3H8": 10.0}          # mol o mol/s; inertes aquí también
fase = "gas"                           # "gas" | "liquido"
P, P0, unidad = 1.0, 1.0, "bar"
temperaturas = [650, 700, 750, 800, 850, 900, 950, 1000]
modo_K = "termo"                       # "directo" | "termo"
K = None                               # directo: {700: [K1, K2], 900: [...]}
K_es_Kc = False                        # True si das K_C en mol/L
especies_termo = ["C3H8", "C2H4", "CH4", "C3H6", "H2"]   # termo: se cargan
datos_termo = None                     # solas de data/sva_tables.json
limitante = "C3H8"
```

Reglas de captura:

- Reacciones **balanceadas por ti**, especies por su **fórmula** (el parser
  deduce los átomos y rechaza reacciones desbalanceadas nombrando el
  elemento). Acepta `"1/2 O2"`, `"0.5 N2"`, `"3H2O"`.
- `modo_K="directo"`: K en forma de fracción mol, `K = Π y^ν (P/P0)^δ`,
  P° = 1 bar. `K_es_Kc=True` si tu dato es `K_C = C_T0^δ Π y^ν` en mol/L
  (se convierte solo, con C_T0 = P0/(RT)).
- `fase="liquido"`: `K = Π x_i^ν`, sin factor de presión (ε = 0); K_C solo
  si δ = 0 (ahí K_C = K_x, C_A0 se cancela).
- Lee el reporte de autovalidación: **[FAIL] = no confíes**, [WARN] =
  revisa a mano, [info] 4/5/12 = compara con tu derivación (δ_j y n_T(ξ),
  la expresión simbólica de cada K, el factor (P/P0)^δ usado).

## Verificación

```bash
python3 tests.py     # casos de respuesta conocida (SVA 13.5/13.9, Fogler,
                     # reformado del curso, termo NH3) — salida 0 = todo PASS
python3 main.py      # suite original de 4 tests del curso
```

## Archivos

| Archivo | Qué hace |
|---|---|
| `equilibrium.py` | Núcleo: Método A (`solve_extents`, avance ξ, multirreacción, log-residuales, multi-arranque) y Método B (`StoichiometricTable`, tabla de Fogler, `solve_Xe` con brentq). Conversión K↔K_C, unidades atm/bar/kPa, checks. |
| `reactions.py` | Parser: `"C3H8 + 3 H2O = 3 CO + 7 H2"` → matriz ν + átomos + inertes; rechaza reacciones desbalanceadas. |
| `thermo.py` | Cp/R = A+BT+CT²+DT⁻² (Tabla C.1), ICPH/ICPS, ΔH°(T), ΔG°(T) (ec. 13.18), K(T), Van't Hoff. R = 8.314 J/mol·K, T₀ = 298.15 K. |
| `exercise.py` | `run_exercise(...)`: el motor de la plantilla (resuelve el barrido, imprime todo, autovalida, CSV, PNG). |
| `selfcheck.py` | `validate` / `validate_sweep` / `solve_and_validate`: los 12 puntos de autovalidación. |
| `plantilla_ejercicio.py` | La que copias para cada ejercicio nuevo. |
| `data/sva_tables.json` | Hf298/Gf298 (C.4) y Cp/R (C.1) por especie, con estatus de verificación. `data/raw/` es para las fotos del libro. |
| `tests.py` | Casos de respuesta conocida (a–i). |
| `main.py` | Suite original (N2O4, SO3, reformado). |
| `ej2_mio.py`, `ej2_reformado_propano.py` | Ejercicio 2 resuelto y validado contra la clave del profesor. |

## Las fórmulas que implementa (las del curso)

- `n_i = n_i0 + Σ_j ν_ij ξ_j`, `n_T = n_T0 + Σ_j δ_j ξ_j`, `y_i = n_i/n_T`
- Gas ideal: `Π y_i^ν_i = K (P/P°)^(−δ)` → residual log:
  `Σ ν ln y + δ ln(P/P°) − ln K = 0` (estable para K de 1e-3 a 1e10)
- `K_C = C_T0^δ Π y^ν`, `C_T0 = P/(RT)`, R = 0.08206 L·atm/mol·K
- Método B (Fogler): `Θ_i`, `δ`, `ε = y_A0 δ`, `C_i(X)` a V constante o
  variable, `X_e` de `K_C = Π C_i^ν` con brentq en `[0, X_max]`
- Termo: `ΔH°(T) = ΔH°298 + R·ICPH`;
  `ΔG°/RT = (ΔG°298−ΔH°298)/(R·T₀) + ΔH°298/(RT) + ICPH/T − ICPS`
  (ec. 13.18 SVA, verificada contra integración numérica de Van't Hoff);
  `K = exp(−ΔG°/RT)`

## Datos que hay que traer para modo "termo" (tablas de SVA)

- Los datos viven en **`data/sva_tables.json`** y se cargan con
  `especies_termo=[...]`. Estado actual: **PROVISIONAL** — 6 especies
  confirmadas con datos del curso, el resto de memoria anclada a checks
  físicos y marcadas `"pendiente_foto"`. Pon las **fotos de las Tablas
  C.1–C.4** en `data/raw/` para transcribirlas y cotejarlas (ver
  `data/raw/LEEME.txt`); el test j vigila que el JSON no se corrompa.
- **Tabla C.4**: ΔH°f,298 y ΔG°f,298 [J/mol] de cada especie (gas ideal;
  para elementos en su estado estándar, H2/N2/O2, valen 0 pero la especie
  **sí** necesita su Cp).
- **Tabla C.1**: A, B, C, D de Cp°/R = A + B·T + C·T² + D·T⁻² por especie
  (ojo con el ×10³ de B, ×10⁶ de C y ×10⁻⁵ de D en el encabezado de la
  tabla: en el JSON se guardan ya multiplicados, p. ej. B = 28.785e-3), y
  su Tmax de validez (el motor avisa si el barrido lo excede).
- Si el enunciado da K_C, tráela con sus unidades y usa `K_es_Kc=True`.
- Si el enunciado da las K directamente (como el reformado de la sesión 6),
  usa `modo_K="directo"` — esas K son dato del problema y no tienen por qué
  coincidir con las que salen de C.1/C.4 (las del reformado del curso, de
  hecho, no coinciden: ver test g).

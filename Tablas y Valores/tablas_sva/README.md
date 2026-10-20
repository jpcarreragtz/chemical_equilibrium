# Base de datos Smith–Van Ness (7a ed., Apéndices B y C)

Transcripción propia de las tablas B.1 y C.1–C.4 para usarlas en el solver
del curso y en Excel.

## Archivos

- **`SVA_Tablas.xlsx`** — hojas: `Indice`, `B1_Propiedades` (M, ω, Tc, Pc,
  Zc, Vc, Tn), `C1_Gases`, `C2_Solidos`, `C3_Liquidos` (Cp/R), `C4_Formacion`
  (ΔH°f,298 y ΔG°f,298), `BD_Gases` (C.1+C.4 unidas para gases) y `Buscar`
  (consulta por especie).
- **`SVA_tablas.json`** — mismo contenido para uso programático. Estructura:
  `{"_meta": {...}, "B1": [filas], "C1_gases": [filas], "C2_solidos": [...],
  "C3_liquidos": [...], "C4_formacion": [filas]}`; cada fila trae `grupo`,
  `nombre` (español), `name` (inglés) y `formula`.

(El solver NO lee estos archivos directamente: usa
`chemical_equilibrium/data/sva_tables.json`, un subconjunto de 42 especies
gaseosas ya en el formato que consume `thermo.py`.)

## Convenciones

- **Coeficientes de Cp/R = A + B·T + C·T² + D·T⁻²** (T en K): las columnas
  `B_e3`, `C_e6`, `D_em5` son las cifras TAL COMO las imprime el libro; las
  columnas `B`, `C`, `D` son los valores "reales" ya multiplicados por sus
  factores de encabezado: B = B_e3×10⁻³, C = C_e6×10⁻⁶, D = D_em5×10⁵.
  `Tmax` es el límite de validez del ajuste.
- **C.4 (formación)**: `Hf298` y `Gf298` en J/mol; el estado estándar va en
  el campo `estado` (`g`, `l`, `s`, `ac`) — ojo con especies dobles, p. ej.
  H2O: `g` = −241 818 y `l` = −285 830 J/mol (para gas ideal usa la `g`).
- **Isómeros**: se distinguen en `formula`/`nombre` (p. ej. isobutano =
  `"C4H10 (iso)"`, «iso-Butano»).
- **Claves del JSON del solver** (`chemical_equilibrium/data/sva_tables.json`,
  dict por fórmula): claves especiales `cC6H12` (ciclohexano) y `C2H4Oox`
  (óxido de etileno); `C2H4O` = acetaldehído; `C6H12` = 1-hexeno.

## Verificación

- Doble pase de transcripción: lectura visual + cotejo automático contra la
  capa de texto del PDF del libro — **0 discrepancias en 161 filas** de
  C.1–C.4 (2026-09-22).
- Checks físicos independientes (vía `chemical_equilibrium/tests.py` y la
  ec. 13.18): K_WGS(1100 K) = 1.006 ≈ 1 (SVA ej. 13.5), K_NH3(298 K) = 762
  (= exp(16450/RT₀)), K_NH3(500 K) = 0.316 vs ~0.32 de literatura.

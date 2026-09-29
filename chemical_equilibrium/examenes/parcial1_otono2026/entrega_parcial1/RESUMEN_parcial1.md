# Resumen — Primer examen parcial, Equilibrio Químico (15/09/2026)

Herramienta: solver propio del curso (`chemical_equilibrium/`), Método A
(avance de reacción, Smith–Van Ness) + Método B (tabla de Fogler), con
autovalidación (balance de elementos, residuos de K, 20 arranques
aleatorios, contraste A vs B). Evidencia de ejecución: `EVIDENCIA_ejecucion.txt`.

---

## Problema 1 — A ⇌ R + 2S (gas, 10 atm, 500 K, v0 = 500 dm³/s de A puro)

**Supuestos:** gas ideal; operación continua isotérmica e isobárica →
volumen variable v = v0(1 + εX); A puro (y_A0 = 1).

**Ecuaciones clave:**
- C_A0 = C_T0 = P/(RT) = 10/(0.08206·500) = 0.2437 mol/dm³; F_A0 = C_A0·v0 = 121.86 mol/s
- δ = 2, ε = y_A0·δ = 2; F_A = F_A0(1−X), F_R = F_A0X, F_S = 2F_A0X, F_T = F_A0(1+2X); v = v0(1+2X)
- K_C = C_R·C_S²/C_A = 4·C_A0²·X³/[(1+2X)²(1−X)] = 0.2  →  cúbica **1.0376X³ − 0.6X − 0.2 = 0**
  (las otras dos raíces son complejas conjugadas −0.446 ± 0.133i → se descartan)

**Resultados:**

| Cantidad | Valor |
|---|---|
| X_Ae | **0.8913** |
| X = 0.70·X_Ae | 0.6239 |
| F_A / F_R / F_S / F_T (mol/s) | 45.83 / 76.03 / 152.07 / 273.93 |
| v (dm³/s) | 1123.9 |
| C_A / C_R / C_S (mol/dm³) | 0.0408 / 0.0677 / 0.1353 |

Comprobación: K_C recalculada en X_Ae = 0.200000; C_T = F_T/v = 0.2437 =
P/RT constante; Métodos A y B coinciden (dif < 1e-15). Selfcheck: PASS.

---

## Problema 2 — Acetaldehído por oxidación de etileno (NO es equilibrio)

**Supuestos:** los avances se despejan de las MEDICIONES de salida
(sistema lineal); reacciones balanceadas por mí (el examen las da
esqueléticas): R1: C₂H₄ + ½O₂ → CH₃CHO; R2: C₂H₄ + 3O₂ → 2CO₂ + 2H₂O.

**Ecuaciones clave:** n_CH3CHO = ξ₁ = 0.78; n_H2O = 2ξ₂ = 0.35 → ξ₂ = 0.175;
n_T = 4.5 − ½ξ₁ = 4.11 mol.

**Resultados:**

| Cantidad | Valor |
|---|---|
| ξ₁, ξ₂ | 0.78, 0.175 mol |
| n_i (C₂H₄/O₂/CH₃CHO/CO₂/H₂O) | 0.545 / 2.085 / 0.78 / 0.35 / 0.35 mol |
| (b) X_C2H4 | **0.6367** |
| (c) y_O2 | **0.5073** |
| (d) Rendimiento al deseado = 0.78/(1.5·1) | **0.52** |
| (e) Selectividad global (deseado/no deseado) = 0.78/0.70 | **1.114** (vs CO₂ 2.229; vs H₂O 2.229) |

Comprobación: átomos C 3→3, H 6→6, O 6→6. Selfcheck (adaptado, sin K): PASS.

---

## Problema 3 — CH₄ + H₂O = CO + 3H₂ ; CO + H₂O = CO₂ + H₂ (1:1, equilibrio)

**Supuestos:** gas ideal; **P = 1 bar (supuesta**, el enunciado no la da; =
P° de las K); H₂O como **vapor** (Tabla C.4, gas ideal: ΔH°f = −241 818
J/mol); K(T) por la ec. 13.18 de SVA con Tablas C.1/C.4; R1 balanceada con
3H₂. Balances: n_T = 2 + 2ξ₁.

ΔH°298: R1 = +205 813, R2 = −41 166 J/mol; ΔG°298: R1 = +141 863, R2 = −28 618 J/mol.

| T (K) | K₁ | K₂ | ξ₁ | ξ₂ | X_CH4 | y (CH₄/H₂O/CO/H₂/CO₂) |
|---|---|---|---|---|---|---|
| 700 | 2.685e-4 | 9.516 | **0.1238** | **0.1160** | 0.124 | 0.390/0.338/0.0035/0.217/0.0516 |
| 1000 | 25.81 | 1.461 | **0.7878** | **0.0644** | 0.788 | 0.0593/0.0413/0.202/0.679/0.0180 |

Newton desde (0.1, 0.01): converge en 7 iteraciones (700 K) y 6 (1000 K)
al mismo resultado que el solver multi-arranque (dif < 1e-15). Selfcheck:
PASS en ambas T y en el barrido.

**Discusión:**

(i) R1 (reformado) es fuertemente **endotérmica** (ΔH° ≈ +220 kJ/mol a las
T del problema). Por Van't Hoff, K₁ sube ~10⁵ veces al pasar de 700 a
1000 K (2.7e-4 → 25.8), y con ella ξ₁ (0.124 → 0.788) y el hidrógeno: y_H2
pasa de 0.217 a 0.679. Favorecer T alta es la decisión de diseño correcta
para producir H₂ (y además δ₁ = +2 favorecería P baja).

(ii) R2 (WGS) es **exotérmica** (ΔH° ≈ −36 kJ/mol): K₂ cae de 9.52 a 1.46
y ξ₂ de 0.116 a 0.064, así que el CO₂ baja (y_CO2 0.0516 → 0.0180) y el CO
sobrevive (y_CO 0.0035 → 0.202). A T alta el reformado "produce CO más
rápido de lo que la WGS lo consume": el gas de síntesis se enriquece en CO.

(iii) La alimentación 1:1 deja al **H₂O como limitante efectivo**: R1 y R2
compiten por la misma agua (n_H2O = 1 − ξ₁ − ξ₂ → a 1000 K queda solo
0.148 mol y X_CH4 se estanca en 0.788). Con relación vapor/metano 4:1
(corrida adicional, mismas K): X_CH4 = 0.31 a 700 K y **0.994 a 1000 K**,
con n_H2 = 3.49 mol vs 2.43 mol (aunque y_H2 = 0.499 por la dilución con
vapor). Por eso industrialmente el reformado usa exceso de vapor: empuja
ambas reacciones y evita además la deposición de carbono.

---

*Datos termodinámicos: `data/sva_tables.json` (Hf/Gf de H₂O, CO y CO₂
confirmados contra los valores del curso usados en este parcial; los
coeficientes de Cp de esas tres especies siguen pendientes de cotejo
contra la Tabla C.1 impresa).*

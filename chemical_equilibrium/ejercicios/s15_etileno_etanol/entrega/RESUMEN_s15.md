# Sesión 15 — Hidratación de etileno con fases L+V (SVA Ej. 13.13)

**Condiciones:** T = 200 °C (473.15 K), P = 34.5 bar (fuente de etileno).
**Regla de fases:** F = 2−2+3−1 = 2 → con T y P fijas las composiciones
quedan determinadas; **no aplica ξ** ni tabla estequiométrica.
**Supuestos:** líquido sin C2H4 (SVA); Lewis-Randall (φ̂ᵢ ≈ φᵢ puro,
virial de Pitzer); Poynting ≈ 1 (H2O 1.0087, EtOH 1.0064);
Wilson etanol(1)/agua(2) DECHEMA.

## Constante de equilibrio (ec. 13.18, tablas verificadas)
ΔH°298 = -45792 J/mol, ΔG°298 = -8378 J/mol,
K298 = 29.37; IDCPH = -17.882 K, IDCPS = -0.057675;
**K(473.15) = 0.03104** (ln K = -3.4726); Van't Hoff ΔH cte: 0.0317.
> Nota: el valor SVA/curso 0.0323 equivale a invertir el signo de la
> corrección de Cp y contradice Van't Hoff (K debe ser < 0.0317 por ser
> exotérmica con |ΔH| creciente); se incluye como sensibilidad y ahí se
> reproducen los números del libro.

## Coeficientes de fugacidad (Pitzer)
φ(T,P): C2H4 0.9634, H2O 0.8451,
EtOH 0.7528 (EtOH fuera de la zona virial → extrapolado);
φsat: H2O 0.9270, EtOH 0.7798.
Wilson: Λ12 = 0.2050, Λ21 = 1.1756.

## Solución (base, K = 0.03104)
| | EtOH | H2O | C2H4 |
|---|---|---|---|
| x (líquido) | 0.0798 | 0.9202 | — |
| γ | 2.588 | 1.018 | — |
| y (vapor) | 0.1874 | 0.4630 | 0.3496 |
| f̂ [bar] | 4.868 | 13.499 | 11.619 |

Residual |ln K_calc − ln K| = 8.9e-16. Checks: PASS.

## Sensibilidad (CSV: sensibilidad_s15.csv)
| caso | x_EtOH | x_H2O | y_C2H4 | y_H2O | y_EtOH | g_EtOH | g_H2O |
|---|---|---|---|---|---|---|---|
| base (K = 0.03104) | 0.0798 | 0.9202 | 0.3496 | 0.4630 | 0.1874 | 2.5881 | 1.0178 |
| gamma = 1 (líq. ideal) | 0.2092 | 0.7908 | 0.4193 | 0.3909 | 0.1898 | 1.0000 | 1.0000 |
| phi = 1 en todo (gas ideal) | 0.0791 | 0.9209 | 0.3978 | 0.4223 | 0.1799 | 2.5966 | 1.0175 |
| Psat_EtOH = 32.6 bar (Antoine extrapolada) | 0.0731 | 0.9269 | 0.3476 | 0.4652 | 0.1872 | 2.6682 | 1.0153 |
| con Poynting | 0.0788 | 0.9212 | 0.3456 | 0.4674 | 0.1871 | 2.6002 | 1.0174 |
| K = 0.0323 (valor SVA/curso) | 0.0832 | 0.9168 | 0.3457 | 0.4619 | 0.1924 | 2.5501 | 1.0192 |

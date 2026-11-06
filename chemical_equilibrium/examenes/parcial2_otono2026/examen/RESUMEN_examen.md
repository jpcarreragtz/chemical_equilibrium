# 2do parcial práctico — Equilibrio Químico (Grupo B, 22/10/2026)

Scripts en esta carpeta (`python3 correr_examen.py` corre los tres y guarda
`salida_examen.txt`). Motor: `chemical_equilibrium/core` sin modificar.
Datos: Tablas C.1/C.4 (`data/sva_tables.json`) y Tabla B.1
(`data/criticas_b1.json`) de Smith–Van Ness–Abbott.

## Resultados

### Problema 1 — C8H10(g) → C8H8(g) + H2(g), 2 mol/s de etilbenceno (`p1_etilbenceno_estireno.py`)

| Inciso | Condición | Modelo | Resultado |
|---|---|---|---|
| datos | 298.15 K | — | ΔH°298 = 117 440 J, ΔG°298 = 83 010 J; ΔA = 4.175, ΔB = −4.766×10⁻³, ΔC = 1.814×10⁻⁶, ΔD = 8300 |
| a) | 800 K, 1 bar | gas ideal | K = 0.05701; **X = 0.2322**; y_EB = 0.6231, y_ST = y_H2 = 0.1885; salida: EB 1.5355, ST 0.4645, H2 0.4645 mol/s; Kφ(800 K, 1 bar) = 0.99994 ≈ 1 |
| b) | X = 0.80, 15 bar | solución ideal | **T = 1174.0 K** (900.8 °C); gas ideal 1173.97 K; φ_EB 0.9936, φ_ST 0.9920, φ_H2 1.0017; Kφ = 1.00004 |
| c) | X = 0.80, 1100 K | solución ideal | K(1100) = 10.879; **P = 6.12 bar** (6.1196); gas ideal 6.1197 bar; Kφ = 0.99998 |

Las dos convenciones de Kφ dan lo mismo: |ΔT| = 0.008 K en b) y
|ΔP| = 0.0002 bar en c).

### Problema 2 — combustible 0.6 CH4 / 0.4 H2 con aire estequiométrico (`p2_combustible_CH4_H2.py`)

| Inciso | Resultado |
|---|---|
| a) | O2 = 0.6·2 + 0.4·0.5 = **1.40 mol**; aire = **6.667 mol**; N2 = **5.267 mol** |
| b) | **HHV = 648.72 kJ/mol** (62.18 MJ/kg); **LHV = 578.30 kJ/mol** (55.43 MJ/kg); HHV − LHV = 1.6·44.012 = 70.42 kJ |
| c) | Productos: CO2 0.6, H2O 1.6, N2 5.267 mol → **T_ad = 2340.0 K** (2066.8 °C); ⟨Cp⟩_H = 283.23 J/K; Q por trayectoria = Q por formación = 0 |

### Problema 3 — oxidación de H2S, R1: H2S + 3/2 O2 → SO2 + H2O, R2: SO2 + 1/2 O2 → SO3 (`p3_oxidacion_H2S.py`)

| Inciso | Resultado |
|---|---|
| a) | ξ1 = 0.8, ξ2 = 0.3. Alimentación: H2S 1.000, O2 3.000, N2 11.286 mol (y = 0.0654 / 0.1963 / 0.7383; 15.286 mol). Salida: H2S 0.2, O2 1.65, SO2 0.5, H2O 0.8, **SO3 0.3**, **N2 11.286** mol (14.736 mol) |
| b) | ΔH°298,1 = −518 018 J, ΔH°298,2 = −98 890 J; ΔH°_R = −34 112 J; Σξ·ΔH° = −444 081 J; ΔH°_P = 390 515 J → **Q = −87.68 kJ por mol de H2S alimentado** (−5.95 kJ por mol de gas de salida). Q < 0: el proceso libera calor y hay que retirarlo para mantener la salida a 830 °C |
| c) | Adiabático con la misma salida: **T_ad = 1269.8 K** (996.6 °C), mayor que 1103.15 K, coherente con Q < 0 |

## Consideraciones

- **P1 a): gas ideal.** A 1 bar, Kφ = 0.99994, así que la corrección por no idealidad es despreciable.
- **P1 b) y c): solución ideal de gases (regla de Lewis-Randall).** φ̂_i = φ_i de la especie pura a T y P. Cada φ sale de la **virial truncada generalizada de Pitzer**: ln φ = (Pr/Tr)(B⁰ + ωB¹), con B⁰ = 0.083 − 0.422/Tr^1.6 y B¹ = 0.139 − 0.172/Tr^4.2. Se cumple Vr ≥ 2 para todas las especies, así que la correlación es válida (Fig. 3.14 de SVA).
- **Convención de Kφ:** la "maestra" (K = K_y·(P/P°)^δ / Kφ) es el default de `core/convenciones.py`. También se reporta la "rigurosa" (SVA 13.28, K = Kφ·K_y·(P/P°)^δ). Como Kφ ≈ 1, las dos coinciden.
- **K(T)** se calcula con la ec. 13.18 de SVA, que integra ΔCp°(T); no es Van't Hoff con ΔH° constante. P° = 1 bar.
- **Aire 21 % O2 / 79 % N2** (P2 y P3). El N2 es inerte.
- **P2:** combustión completa (todo el C pasa a CO2 y no se forma CO), con aire teórico, así que no queda O2 en los productos. El H2O de los productos es **vapor**, por eso la T de flama se calcula con el LHV.
- **Cp extrapolado en P2 c):** T_ad = 2340 K está por encima de Tmax = 2000 K de la Tabla C.1 para N2, CO2 y H2O. El resultado es una extrapolación. En P3 c) (1269.8 K) no hay extrapolación.
- **P3:** solo entran H2S y aire, sin agua. Las reacciones se balancearon a R1: H2S + 3/2 O2 → SO2 + H2O y R2: SO2 + 1/2 O2 → SO3.
- **P3 c):** las conversiones (ξ1 = 0.8, ξ2 = 0.3) y la composición de salida se mantienen constantes en el caso adiabático. Solo cambia T_salida.
- **Balances de energía:** Q = ΔH, con W_s = 0 y ΔEc ≈ ΔEp ≈ 0. La trayectoria pasa por 298.15 K y se comprueba con el método de entalpías de formación; los dos métodos coinciden a precisión de máquina.

## Verificación

- `correr_examen.py`: los 3 scripts dan **VEREDICTO: OK**, con 38 checks en PASS y 0 en FAIL.
- `python3 tests/tests.py`: todos los tests en PASS.
- Ningún resultado difiere de los valores esperados en más de 0.1 %.

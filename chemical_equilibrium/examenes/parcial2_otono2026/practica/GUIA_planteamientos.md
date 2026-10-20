# 2do parcial EQ. Químico — planteamientos a mano (ejercicios adicionales)

Temas del examen (Sesión 21): **avance de reacción con efecto de la presión
(con fugacidad)**, **combustión**, **HHV y LHV**, **balance de energía /
temperatura de flama adiabática**.
Datos: Tablas C.1 (Cp/R) y C.4 (ΔH°f, ΔG°f) de Smith–Van Ness (en
`data/sva_tables.json`), Tabla B.1 (Tc, Pc, ω) en `data/criticas_b1.json`.

| Problema | Script | Resultado del código | Maestra | Estado |
|---|---|---|---|---|
| 1 a) | `p1_amoniaco_presion.py` | y_NH3 = 0.9349 | 0.9349 | ✔ |
| 1 b) | | T = 399.50 K | 399.5 K | ✔ |
| 1 c) | | T = 577.64 K | 577.6 K (compañero 551) | ✔ |
| 1 d) | | **T = 570.4 K** (convención "maestra", default) | 568.6 K (compañero 557) | ✔ (dif 1.8 K = término ω·B¹; ver nota) |
| 2 a) | `p2_isomerizacion_margules.py` | x_A = 0.3955 | 0.3955 | ✔ |
| 2 b) | | x_A = 0.4005 (error ideal +1.27 %) | 0.4005 | ✔ |
| 3 a) | `p3_gas_natural_hhv.py` | −921.71 kJ/mol | −921.71 | ✔ |
| 3 b) | | −946.19 kJ/mol | −946.16 | ✔ (redondeo) |
| 3 c) | | **−932.87 kJ/mol** | −923.87 | ⚠ dígitos 3↔2 |
| 4 | `p4_reformado_calor.py` | Q = 54 881 J/mol producto | 54 881 J | ✔ |
| 5 | `p5_calentador_flujo_gas.py` | V̇ = 0.407 m³/s | 0.407 m³/s | ✔ |

Correr todo: `python3 correr_todo.py` (resume veredictos y guarda `salida_practica.txt`).

---

## Fórmulas base (las que se usan en todo)

* **K y presión (gas):** K = F · K_y · (P/P°)^δ, con K_y = Π y_i^ν_i, K_φ = Π φ̂_i^ν_i, δ = Σν_i.
  - Gas ideal: K_φ = 1 (F = 1). Solución ideal de gases (Lewis-Randall): φ̂_i = φ_i (especie pura a T, P).
  - **Convención de la maestra (default en `core/convenciones.py`):** F = 1/K_φ ⇒ K = K_y·(P/P°)^δ / K_φ ⇔ **K_y = K·(P/P°)^(−δ)·K_φ**.
  - Convención rigurosa (f̂_i = φ_i y_i P, SVA 13.28): F = K_φ ⇒ K_y = K·(P/P°)^(−δ)/K_φ. Disponible con `convencion_Kphi="rigurosa"`.
* **φ por la virial de Pitzer (Sesión 13):** ln φ = (P_r/T_r)(B⁰ + ωB¹), B⁰ = 0.083 − 0.422/T_r^1.6, B¹ = 0.139 − 0.172/T_r^4.2.
* **K(T) rigurosa (ec. 13.18 SVA):**
  −ln K = ΔG°/RT = (ΔG°₀ − ΔH°₀)/(RT₀) + ΔH°₀/(RT) + (1/T)·∫(ΔCp/R)dT − ∫(ΔCp/R)dT/T
  con ∫(ΔCp/R)dT = ΔA(T−T₀) + ΔB/2(T²−T₀²) + ΔC/3(T³−T₀³) + ΔD(T−T₀)/(T·T₀).
* **Líquido (Sesión 15/EQ. Físico):** a_i = x_i γ_i ⇒ K = Π (x_i γ_i)^ν_i (P no entra).
  Margules 1 par.: ln γ₁ = A x₂², ln γ₂ = A x₁². 2 par.: ln γ₁ = x₂²[A₁₂ + 2(A₂₁−A₁₂)x₁], ln γ₂ = x₁²[A₂₁ + 2(A₁₂−A₂₁)x₂].
* **Combustión / HHV-LHV:** ΔH°c = Σν ΔH°f; HHV = −ΔH°c con H₂O(l) (ΔH°f = −285 830 J/mol); LHV con H₂O(g) (−241 818); HHV = LHV + n_H₂O·44 012 J/mol; mezcla: Σ y_i HHV_i. Dulong: HHV(MJ/kg) = 33.9C + 144.2(H − O/8) + 9.4S; LHV = HHV − 2.442·(9H).
* **Balance de energía (W_s = 0):** Q = ΔH = ΔH°_R + Σ ξ_j ΔH°298,j + ΔH°_P
  ΔH°_R = −Σ n_i,ent ∫_{298}^{T_ent} Cp dT; ΔH°_P = Σ n_i,sal ∫_{298}^{T_sal} Cp dT. Flama adiabática: Q = 0 ⇒ despejar T_sal (iterando con ⟨Cp⟩_H como en tu Excel s18).

---

## Problema 1 — Síntesis de NH₃ (efecto de P y fugacidad)

½N₂ + 3/2 H₂ → NH₃, alimentación estequiométrica (base n_N₂ = ½, n_H₂ = 3/2).

1. **Tabla en ξ:** n_NH₃ = ξ, n_N₂ = ½(1−ξ), n_H₂ = 3/2(1−ξ), n_T = 2 − ξ, δ = −1.
2. **K_y = y_NH₃ / (y_N₂^½ y_H₂^3/2)** = [ξ/(2−ξ)] / {[½(1−ξ)/(2−ξ)]^½ [3/2(1−ξ)/(2−ξ)]^3/2}.
3. **K(T)** (ec. 13.18): ΔH°₀ = −46 110, ΔG°₀ = −16 450 J/mol; ΔA = −2.9355, ΔB = 2.0905×10⁻³, ΔC = 0, ΔD = −0.3305×10⁵.
4. **a)** K(300 K) = 679.57; a 1 bar K_y = K → ξ = 0.96636 → **y_NH₃ = 0.9349**.
5. **y_NH₃ = 0.5 ⇒ ξ = 2/3**, y_N₂ = 0.125, y_H₂ = 0.375 ⇒ **K_y = 6.1584** (fijo: sólo cambia qué K se necesita).
6. **b)** P = 1 bar: K = K_y = 6.1584 → resolver K(T) = 6.1584 ⇒ **T = 399.5 K**.
7. **c)** P = 100 bar, GI: K = K_y·(100)⁻¹ = 0.061584 ⇒ **T = 577.6 K**.
8. **d)** Solución ideal (convención de la maestra): K_y = K·(P/P°)·K_φ ⇒ K requerida = K_y/(100·K_φ), K_φ = φ_NH₃/(φ_N₂^½ φ_H₂^3/2), φ a (T, 100 bar) → iterar T:
   en T = 570.4 K: φ_NH₃ = 0.9173, φ_N₂ = 1.0334, φ_H₂ = 1.0217 ⇒ K_φ = 0.8737 ⇒ K = 0.06158/0.8737 = 0.07048 ⇒ **T = 570.4 K** (resultado oficial).

> **Nota d):** el 568.6 K de la clave es el mismo cálculo omitiendo el término ω·B¹ en ln φ (ln φ = P_r/T_r·B⁰); la diferencia (1.8 K) es solo ese término. La convención **rigurosa** (f̂_i = φ_i y_i P, SVA 13.28: K_y = K·P/K_φ) da **584.5 K** > T_c: con K_φ < 1 la K requerida es menor y la no idealidad favorece al NH₃ (proceso Haber). Queda disponible con `convencion_Kphi="rigurosa"`; el script y los tests reportan ambas. φ_NH₃ por Lee-Kesler (Tablas E.13/E.14) es 0.927 vs 0.925 de Pitzer, así que el método de φ no explica la diferencia entre convenciones. Para c) tu compañero (551 K) está mal: la clave y el código coinciden en 577.6 K.

## Problema 2 — A(l) → B(l) con Margules

1. ln K = −ΔG°/RT = 1000/(8.314·298.15) = 0.40342 ⇒ K = 1.4969.
2. K = x_Bγ_B/(x_Aγ_A); Margules 1 par. (A = 0.1): ln γ_A = 0.1x_B², ln γ_B = 0.1x_A².
3. ln[(1−x_A)/x_A] + 0.1(x_A² − x_B²) = ln K, y como x_A² − x_B² = 2x_A − 1:
   **ln[(1−x_A)/x_A] + 0.1(2x_A − 1) = 0.40342 ⇒ x_A = 0.3955** (γ_A = 1.0372, γ_B = 1.0158).
4. Ideal: x_A = 1/(1+K) = **0.4005** ⇒ error +1.27 % en x_A (−0.83 % en x_B).

## Problema 3 — Calor de combustión de gases naturales (H₂O líquida)

1. CH₄ + 2O₂ → CO₂ + 2H₂O(l): ΔH°c = −393 509 + 2(−285 830) + 74 520 = **−890 649 J/mol**
   C₂H₆ + 7/2 O₂ → 2CO₂ + 3H₂O(l): **−1 560 688 J/mol**; C₃H₈ + 5O₂ → 3CO₂ + 4H₂O(l): **−2 219 167 J/mol**. N₂ no aporta.
2. **Expresión:** ΔH°c = −890.649 y_CH₄ − 1560.688 y_C₂H₆ − 2219.167 y_C₃H₈ (kJ/mol de gas).
3. a) −921.71; b) **−946.19 (el mayor)**; c) −932.87 kJ/mol.

> **⚠ Nota c):** la clave −923.87 es −932.87 con los dígitos 3 y 2 intercambiados. No cambia la conclusión (b es el mayor).

## Problema 4 — Reformado (2 reacciones), Q por mol de producto

Base 1 mol de producto: CO₂ 0.0275, CO 0.1725, H₂O 0.1725, H₂ 0.6275.
1. **Balance atómico:** C: n_CH₄ = 0.0275 + 0.1725 = 0.2; H: 4(0.2) + 2n_H₂O = 2(0.1725) + 2(0.6275) ⇒ n_H₂O = 0.4; O: 0.4 = 2(0.0275) + 0.1725 + 0.1725 ✔.
2. **Avances:** ξ₁ = 0.2 (CH₄ consumido), ξ₂ = 0.0275 (CO₂ formado).
3. **ΔH°298:** R1 = −110 525 + 74 520 + 241 818 = 205 813 J; R2 = −393 509 + 110 525 + 241 818 = −41 166 J ⇒ Σξ·ΔH = 41 162.6 − 1 132.1 = 40 031 J.
4. **ΔH°_R** (0.2 CH₄ + 0.4 H₂O, 773.15 → 298.15 K) = −11 455 J; **ΔH°_P** (productos 298.15 → 1123.15 K) = 26 305 J.
5. **Q = −11 455 + 40 031 + 26 305 = 54 881 J** por mol de gas producido ✔ (comprobado también por entalpías de formación).

## Problema 5 — Calentador: flujo volumétrico de combustible

Base 1 mol de combustible (0.94 CH₄, 0.06 N₂).
1. ξ₁ (→CO₂) = 0.70·0.94 = 0.658; ξ₂ (→CO) = 0.30·0.94 = 0.282.
2. O₂ teórico (combustión completa) = 1.88 ⇒ alimentado 1.35·1.88 = 2.538; N₂ del aire = 9.5477.
3. Salida: CO₂ 0.658, CO 0.282, H₂O(g) 1.88, O₂ 2.538 − 2(0.658) − 1.5(0.282) = 0.799, N₂ 9.6077.
4. Entradas a 77 °F = 298.15 K ⇒ ΔH°_R = 0.
   ΔH°298: 0.658(−802 625) + 0.282(−519 641) = −674 666 J; ΔH°_P(298.15 → 483.15 K) = 75 412 J ⇒ ΔH = −599 254 J.
   **Calor al agua = 599 254 J por mol de combustible.**
5. Agua: ṁ = 75 lbm/s = 34.019 kg/s; Δh = h_l(95 °C) − h_l(25 °C) = 398.09 − 104.83 = 293.26 kJ/kg ⇒ Q = 9 976.5 kW.
6. ṅ = 9.9765×10⁶/599 254 = 16.648 mol/s ⇒ **V̇ = ṅRT/P (77 °F, 1 atm) = 0.407 m³/s** (14.4 ft³/s) ✔.

---

## Checklist para el examen práctico

1. Escribe reacciones balanceadas y ν; δ = Σν (sólo especies gas).
2. Tabla en ξ (o balance atómico si te dan la salida); comprueba balance de C, H, O, N.
3. ¿Qué te piden? y a T dada → resolver ξ; T para y objetivo → K_y fijo, resolver K(T).
4. ¿Gas ideal, solución ideal (φ puro) o líquido (γ)? Escribe K = F·K_y·(P/P°)^δ (F = 1/K_φ con la convención de la maestra, `core/convenciones.py`) o K = Π(xγ)^ν.
5. Energía: base de cálculo, T de cada corriente, trayectoria por 298.15 K, ΔH°298 por reacción, ΔH°_R (signo −), ΔH°_P.
6. Corre el script/plantilla y verifica que **todos los checks den PASS** y que la tendencia física tenga sentido.

Plantillas: `plantilla/plantilla_ejercicio.py` (equilibrio, barrido de T;
`fugacidad="lewis"` opcional), `plantilla/plantilla_presion_fugacidad.py`
(1 reacción gas: y a (T, P) o T para y objetivo, gi | Lewis-Randall) y
`plantilla/plantilla_balance_energia.py` (Q, T_ad, HHV/LHV).

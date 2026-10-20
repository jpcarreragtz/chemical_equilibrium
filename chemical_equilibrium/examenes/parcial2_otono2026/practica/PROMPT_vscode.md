# Prompt para Claude en VS Code (día del examen práctico)

Pega esto en Claude Code (VS Code) con la carpeta `chemical_equilibrium/`
abierta, y debajo el enunciado del problema:

---
Estoy en el examen práctico del 2do parcial de Equilibrio Químico (IBERO).
Usa SOLO el motor de este repo: `core/thermo.py` (K(T) ec. 13.18 SVA),
`core/nonideal.py` (phi Pitzer, K_phi_lewis, Margules, Wilson),
`core/energia.py` (balance_energia, T_salida_para_Q, combustion_mezcla,
alimentacion_por_atomos, avances_desde_flujos, h_liquido_sat) y los datos
de `data/sva_tables.json` y `data/criticas_b1.json`. Antes de programar,
escribe el planteamiento a mano (reacciones y ν, δ, tabla en ξ o balance
atómico, expresión de K o de Q = ΔH°_R + Σ ξ ΔH°298 + ΔH°_P). Luego crea
un script nuevo en `examenes/parcial2_otono2026/` copiando el estilo de
`practica/p1…p5_*.py` (bloque DATOS editable, impresión paso a paso y
sección AUTOVALIDACIÓN con PASS/FAIL: balance atómico, residual, Q por
trayectoria = Q por entalpías de formación, tendencia física). Para
equilibrio a varias T usa `plantilla/plantilla_ejercicio.py`; para
Q / T de flama / HHV usa `plantilla/plantilla_balance_energia.py`.
Córrelo con `python3`, corre también `python3 tests/tests.py`, y dame el
resultado final con unidades. No modifiques `core/` salvo que encuentres
un bug (y dímelo).

ENUNCIADO:
<pega aquí>
---

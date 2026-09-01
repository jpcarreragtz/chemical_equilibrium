"""PLANTILLA para un ejercicio NUEVO de equilibrio químico (Método A).

Cómo usarla
-----------
1. Copia este archivo:   cp plantilla_ejercicio.py ej<N>_<nombre>.py
2. Rellena SOLO los bloques marcados con  # <-- EDITA
3. Corre:                python3 ej<N>_<nombre>.py
4. Lee el reporte de AUTOVALIDACIÓN (selfcheck), puntos 1-12:
     [FAIL]  => NO confíes en el resultado: revisa datos / planteamiento.
     [WARN]  => sospechoso: revísalo a mano antes de confiar.
     [info]  => compáralo contra tu derivación a mano:
                 punto 4: delta_j y n_T(xi)
                 punto 5: expresión simbólica de cada K
                 punto 12: factor (P/P0)^delta realmente usado
   Si todo es PASS/info, el resultado tiene: reacciones balanceadas,
   residuos de K en cero, átomos conservados, la misma solución desde
   arranques aleatorios y (si es 1 reacción) el mismo resultado por
   Método B — es decir, es confiable aunque no tengas clave de respuestas.

Los datos precargados son un EJEMPLO COMPLETO Y FUNCIONAL (síntesis de
amoniaco con argón inerte a 10 bar) para que la plantilla corra tal cual;
reemplázalos por los de tu ejercicio.

Convención de K: forma de fracción mol,
    K_j = prod_i y_i^nu_ij * (P/P0)^delta_j .
(Si tu dato es Kc en mol/L, conviértelo antes con equilibrium.Kc_to_Ky.)
"""

from selfcheck import Problem, solve_and_validate, validate_sweep

# ===== 1) ESPECIES — el ORDEN define las columnas de nu ====================
especies = ["N2", "H2", "NH3", "Ar"]                           # <-- EDITA

# ===== 2) MOLES INICIALES (mol; especie ausente en el dict = 0) ===========
n0 = {"N2": 1.0, "H2": 3.0, "NH3": 0.0, "Ar": 0.1}             # <-- EDITA

# ===== 3) MATRIZ ESTEQUIOMÉTRICA — una FILA por reacción ==================
# Reactivos < 0, productos > 0, no participa = 0. MISMO orden que `especies`.
#          N2   H2  NH3   Ar
nu = [
    [   -1,  -3,   2,   0 ],   # R1: N2 + 3 H2 <-> 2 NH3      # <-- EDITA
    # [ ...                ],  # R2: ...   (agrega una fila por reacción)
]

# ===== 4) INERTES — especies presentes en n0 que NO reaccionan ============
inertes = ["Ar"]                                               # <-- EDITA

# ===== 5) COMPOSICIÓN ATÓMICA — para el checker de balance ================
atomos = {                                                     # <-- EDITA
    "N2":  {"N": 2},
    "H2":  {"H": 2},
    "NH3": {"N": 1, "H": 3},
    "Ar":  {"Ar": 1},
}

# ===== 6) CONSTANTES DE EQUILIBRIO — una lista [K1, K2, ...] por T ========
casos = {                                                      # <-- EDITA
    600.0: [1.66e-3],
    800.0: [6.00e-6],
}

# ===== 7) PRESIÓN — P y P0 en la MISMA unidad =============================
P      = 10.0     # presión del sistema                        # <-- EDITA
P0     = 1.0      # presión de referencia de las K             # <-- EDITA
unidad = "bar"    # "bar" | "atm" | "kPa"                      # <-- EDITA

# ==========================================================================
# De aquí para abajo NO necesitas editar nada:
# resuelve cada T, autovalida (puntos 1-10 y 12) e imprime resultados;
# al final corre el chequeo de Le Chatelier (punto 11) sobre el barrido.
# ==========================================================================
xi_por_T = {}
for T, K in casos.items():
    print(f"\n{'=' * 64}\n  T = {T:g} K\n{'=' * 64}")
    prob = Problem(species=especies, n0=n0, nu=nu, K=K, P=P, P0=P0,
                   atoms=atomos, inerts=inertes, T=T, pressure_unit=unidad)
    res = solve_and_validate(prob)   # resuelve y AUTOVALIDA en un paso
    xi_por_T[T] = tuple(res.xi)

    print("\n  Resultados:")
    for j, xi in enumerate(res.xi):
        print(f"    xi_{j + 1} = {xi:.6f} mol")
    print(f"    n_T  = {res.n_total:.6f} mol")
    print(f"    {'Especie':8s} {'n_i (mol)':>12s} {'y_i':>10s}")
    for sp in especies:
        print(f"    {sp:8s} {res.moles[sp]:12.6g} "
              f"{res.mole_fractions[sp]:10.6f}")

# ---- Punto 11: Le Chatelier sobre el barrido (>= 2 temperaturas) ---------
if len(casos) >= 2:
    validate_sweep(casos, xi_por_T)

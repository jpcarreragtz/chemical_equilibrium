"""P3 examen: descomposición térmica de propano.

    R1: C3H8 = C2H4 + CH4       (craqueo)
    R2: C3H8 = C3H6 + H2        (deshidrogenación)

10 mol/s de propano a un reactor de flujo continuo; salida en equilibrio.
Composición de salida para 8 temperaturas entre 650 y 1000 K a 1 bar,
con K(T) calculadas de las Tablas C.4 y C.1 de SVA (ec. 13.18), gráfica
y_i vs T y CSV.

Esperado: dH°_R1 y dH°_R2 > 0 (endotérmicas) -> K1 y K2 crecen con T;
X_C3H8 creciente (~0.86 a 650 K -> ~1.0 a 1000 K); selfcheck PASS en
todas las T.
"""

from exercise import run_exercise

# ══════════════════════════ INPUTS ════════════════════════════════════════
titulo = "P3 examen: descomposicion termica de propano"

reacciones = [
    "C3H8 = C2H4 + CH4",             # R1: craqueo
    "C3H8 = C3H6 + H2",              # R2: deshidrogenación
]

alimentacion = {"C3H8": 10.0}        # mol/s

fase = "gas"
P, P0, unidad = 1.0, 1.0, "bar"
temperaturas = [650, 700, 750, 800, 850, 900, 950, 1000]

modo_K = "termo"
K = None
K_es_Kc = False

# Smith-Van Ness Tablas C.4 (Hf, Gf en J/mol) y C.1 (Cp/R = A+BT+CT²+DT⁻²)
datos_termo = {
    "C3H8": {"Hf298": -104680.0, "Gf298": -24290.0,
             "A": 1.213, "B": 28.785e-3, "C": -8.824e-6, "D": 0.0},
    "C2H4": {"Hf298": 52510.0, "Gf298": 68460.0,
             "A": 1.424, "B": 14.394e-3, "C": -4.392e-6, "D": 0.0},
    "CH4":  {"Hf298": -74520.0, "Gf298": -50460.0,
             "A": 1.702, "B": 9.081e-3, "C": -2.164e-6, "D": 0.0},
    "C3H6": {"Hf298": 19710.0, "Gf298": 62205.0,
             "A": 1.637, "B": 22.706e-3, "C": -6.915e-6, "D": 0.0},
    "H2":   {"Hf298": 0.0, "Gf298": 0.0,
             "A": 3.249, "B": 0.422e-3, "C": 0.0, "D": 0.083e5},
}

limitante = "C3H8"
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

resultado = run_exercise(
    titulo=titulo, reacciones=reacciones, alimentacion=alimentacion,
    fase=fase, P=P, P0=P0, unidad=unidad, temperaturas=temperaturas,
    modo_K=modo_K, K=K, datos_termo=datos_termo, K_es_Kc=K_es_Kc,
    limitante=limitante,
)

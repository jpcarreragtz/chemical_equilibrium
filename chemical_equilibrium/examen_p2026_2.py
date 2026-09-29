"""P2 examen: oxidación de n-butano a anhídrido maleico.

    2 C4H10 + 7 O2 = 2 C4H2O3 + 8 H2O      (delta = +1)

Entran 10 mol de n-C4H10 con aire en exceso (relación butano:O2 = 1:5,
N2 = 50·0.79/0.21 = 188.095 mol como inerte). P = 1 bar, T = 350 K,
K(350) = 490.
(a) avance de reacción, (b) fracciones mol en el equilibrio,
(c) conversión del limitante (n-butano).

Esperado: xi = 4.392, y_C4H2O3 = 0.0348, y_H2O = 0.139, X_C4H10 = 0.878.
"""

from exercise import run_exercise

# ══════════════════════════ INPUTS ════════════════════════════════════════
titulo = "P2 examen: n-butano a anhidrido maleico"

reacciones = ["2 C4H10 + 7 O2 = 2 C4H2O3 + 8 H2O"]

alimentacion = {"C4H10": 10.0, "O2": 50.0, "N2": 188.095}   # N2 inerte (aire)

fase = "gas"
P, P0, unidad = 1.0, 1.0, "bar"
temperaturas = [350.0]

modo_K = "directo"
K = {350.0: [490.0]}
K_es_Kc = False

datos_termo = None
limitante = "C4H10"
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

resultado = run_exercise(
    titulo=titulo, reacciones=reacciones, alimentacion=alimentacion,
    fase=fase, P=P, P0=P0, unidad=unidad, temperaturas=temperaturas,
    modo_K=modo_K, K=K, datos_termo=datos_termo, K_es_Kc=K_es_Kc,
    limitante=limitante,
)

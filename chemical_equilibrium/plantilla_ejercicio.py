"""PLANTILLA — llena UN solo bloque de inputs y corre. Todo lo demás
(matriz nu, átomos, inertes, K(T) si es modo termo, solución, X, tablas,
autovalidación, CSV y gráfica y_i vs T) se calcula solo.

Cómo usarla
-----------
1. cp plantilla_ejercicio.py ej<N>_<nombre>.py
2. Edita SOLO el bloque INPUTS de abajo.
3. python3 ej<N>_<nombre>.py
4. Lee el reporte de AUTOVALIDACIÓN al final de cada T:
   [FAIL] = no confíes; [WARN] = revísalo a mano; [info] 4/5/12 = compara
   con tu derivación (delta_j y n_T, expresión de cada K, factor (P/P0)^δ).

Reglas de captura
-----------------
* Reacciones como texto CON coeficientes: "C3H8 + 3 H2O = 3 CO + 7 H2".
  Especies por su FÓRMULA (C3H8, H2O...), no por nombre. Acepta "1/2 O2".
* modo_K = "directo": das K = {T: [K1, K2, ...]} (fracción mol, con el
  factor (P/P0)^delta incluido en la definición, P0 = 1 bar del curso).
* modo_K = "termo": lista las especies en especies_termo = [...] y los
  datos (Hf298/Gf298 de la Tabla C.4, Cp/R de la C.1) se cargan solos de
  data/sva_tables.json; K(T) se calcula con la ec. 13.18 para cada T.
  (También puedes pasar datos_termo = {...} a mano; gana sobre el JSON.)
* K_es_Kc = True si tus K son K_C en mol/L (se convierten con
  K = K_C/C_T0^delta, C_T0 = P0/(R·T)); en fase líquida solo si delta = 0.
* fase = "liquido": se usa K = Π x_i^ν y el factor de presión no aplica.

El ejemplo precargado es el CASO DEL EXAMEN: pirólisis de propano,
C3H8 = C2H4 + CH4 y C3H8 = C3H6 + H2, 10 mol/s de C3H8, 1 bar,
T = 650–1000 K, con K calculadas de las Tablas C.1/C.4 de SVA.
"""

from exercise import run_exercise

# ══════════════════════════ INPUTS (edita SOLO esto) ══════════════════════
titulo = "Pirólisis de propano (caso examen)"

reacciones = [                       # texto con coeficientes, ya balanceadas
    "C3H8 = C2H4 + CH4",             # R1: craqueo
    "C3H8 = C3H6 + H2",              # R2: deshidrogenación
]

alimentacion = {"C3H8": 10.0}        # mol (o mol/s); inertes van aquí también

fase = "gas"                         # "gas" | "liquido"
P, P0, unidad = 1.0, 1.0, "bar"      # P y P0 en la MISMA unidad
temperaturas = [650, 700, 750, 800, 850, 900, 950, 1000]   # K

modo_K = "termo"                     # "directo" | "termo"

# --- modo "directo": K por temperatura, una lista por T (orden = reacciones)
K = None
# K = {700: [7.2285e-3, 7.5034], 900: [1.3407e6, 1.5591]}
K_es_Kc = False                      # True si las K de arriba son K_C (mol/L)

# --- modo "termo": se cargan Hf298/Gf298/Cp de data/sva_tables.json
#     (lista TODAS las especies de tus reacciones; error claro si falta una)
especies_termo = ["C3H8", "C2H4", "CH4", "C3H6", "H2"]
datos_termo = None                   # o pásalos a mano: {especie: {Hf298,
                                     #  Gf298, A, B, C, D}} (gana sobre
                                     #  especies_termo)

limitante = "C3H8"                   # None = elige el 1er reactivo de R1
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

resultado = run_exercise(
    titulo=titulo, reacciones=reacciones, alimentacion=alimentacion,
    fase=fase, P=P, P0=P0, unidad=unidad, temperaturas=temperaturas,
    modo_K=modo_K, K=K, datos_termo=datos_termo,
    especies_termo=especies_termo, K_es_Kc=K_es_Kc,
    limitante=limitante,
)

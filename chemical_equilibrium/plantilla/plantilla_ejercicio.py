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

El archivo está dividido en celdas `# %%` para el Interactive Window de
VS Code (Shift+Enter sobre una celda): las celdas de tabla y gráfica solo
muestran cosas ahí; como script (`python3 archivo.py`) no imprimen nada y
la salida es exactamente la de siempre.

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

# %% Setup — shim de sys.path e imports
import sys, pathlib
try:
    _RAIZ = pathlib.Path(__file__).resolve().parents[1]          # raíz del proyecto
except NameError:                                                # celda/notebook sin __file__
    _RAIZ = next((p for p in [pathlib.Path.cwd(), *pathlib.Path.cwd().parents]
                  if (p / "core" / "exercise.py").exists()), None)
    if _RAIZ is None:
        raise RuntimeError("no encuentro core/ desde el cwd: corre el kernel "
                           "con la carpeta del script como directorio de trabajo")
sys.path.insert(0, str(_RAIZ))

from core.exercise import run_exercise

try:
    import pandas as pd              # solo para las tablas inline (pip3 install pandas)
except ImportError:                  # sin pandas el script corre igual; las
    pd = None                        # celdas de tabla no muestran nada

# %% [markdown]
# ## Enunciado
#
# _(Escribe aquí el enunciado tal como viene en la tarea/examen.)_
#
# **Reacciones**
#
# - R1: $\mathrm{C_3H_8 \rightleftharpoons C_2H_4 + CH_4}$ &nbsp; (δ₁ = +1)
# - R2: $\mathrm{C_3H_8 \rightleftharpoons C_3H_6 + H_2}$ &nbsp; (δ₂ = +1)
#
# **Alimentación y condiciones:** 10 mol/s de C₃H₈, P = 1 bar, T = 650–1000 K.
#
# **Balances con ξ:** $n_T = n_0 + \delta_1\,\xi_1 + \delta_2\,\xi_2$
# (revisa el `[info] 4` del selfcheck contra tu derivación a mano).

# %% Reacciones y alimentación — INPUTS (edita SOLO esto)
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

# %% Resolver — K(T) (ec. 13.18 si modo termo) + solver + selfcheck + CSV + PNG
# run_exercise hace todo en una llamada e imprime el reporte completo
# (K por T, ξ, n_i, y_i, X, los 12 puntos de autovalidación, resumen del
# barrido). Las celdas de abajo solo RE-LEEN lo que devuelve; no recalculan.
resultado = run_exercise(
    titulo=titulo, reacciones=reacciones, alimentacion=alimentacion,
    fase=fase, P=P, P0=P0, unidad=unidad, temperaturas=temperaturas,
    modo_K=modo_K, K=K, datos_termo=datos_termo,
    especies_termo=especies_termo, K_es_Kc=K_es_Kc,
    limitante=limitante,
)

# %% K(T) — constantes usadas en cada T (en modo termo, con ΔH° y ΔG° de la ec. 13.18)
por_T = resultado["resultados"]          # {T: {"res", "K", "X", "termo", "reporte"}}
filas_K = []
for T in sorted(por_T):
    fila = {"T (K)": T}
    fila.update({f"K{j + 1}": k for j, k in enumerate(por_T[T]["K"])})
    if por_T[T]["termo"] is not None:    # modo termo: ΔH°(T), ΔG°(T) por reacción
        for j, t in enumerate(por_T[T]["termo"]):
            fila[f"dH°_R{j + 1} (kJ/mol)"] = t["dH"] / 1000.0
            fila[f"dG°_R{j + 1} (kJ/mol)"] = t["dG"] / 1000.0
    filas_K.append(fila)
tabla_K = pd.DataFrame(filas_K).set_index("T (K)") if pd else filas_K
tabla_K

# %% Tabla de resultados — ξ_j, X del limitante, n_i y fracciones mol por T
especies = resultado["especies"]
frac = "x" if fase != "gas" else "y"
filas = []
for T in sorted(por_T):
    r = por_T[T]["res"]                  # ExtentResult del solver
    fila = {"T (K)": T}
    fila.update({f"xi_{j + 1}": x for j, x in enumerate(r.xi)})
    fila[f"X_{resultado['limitante']}"] = por_T[T]["X"]
    fila.update({f"n_{sp}": r.moles[sp] for sp in especies})
    fila.update({f"{frac}_{sp}": r.mole_fractions[sp] for sp in especies})
    filas.append(fila)
tabla_resultados = pd.DataFrame(filas).set_index("T (K)") if pd else filas
tabla_resultados

# %% Gráfica — y_i vs T (la misma PNG que run_exercise ya guardó en salidas/)
# core dibuja con el backend Agg y cierra la figura, así que aquí se muestra
# el PNG guardado tal cual (idéntico al archivo) en vez de redibujar.
try:
    from IPython.display import Image
    grafica = Image(filename=resultado["png"]) if resultado["png"] else None
except ImportError:
    grafica = None
grafica

# %% Selfcheck — resumen de los 12 puntos por T (el detalle ya se imprimió en "Resolver")
filas_sc = [{"T (K)": T,
             "PASS": por_T[T]["reporte"].n_pass,
             "WARN": por_T[T]["reporte"].n_warn,
             "FAIL": por_T[T]["reporte"].n_fail,
             "ok": por_T[T]["reporte"].ok,
             "veredicto global": resultado["veredicto"]}
            for T in sorted(por_T)]
tabla_selfcheck = pd.DataFrame(filas_sc).set_index("T (K)") if pd else filas_sc
tabla_selfcheck

"""P1 examen: hidrogenación de o-cresol a 2-metilciclohexanona.

    C7H8O + 2 H2 = C7H12O      (gas, flujo continuo, isotérmico e isobárico)

Alimentación equimolar (1:1), 2 atm, 300 °C = 573.15 K.
(a) Tabla estequiométrica, C_i y P_i a X = 90 % (base: H2, el limitante).
(b) K_C = 325 a 300 °C -> conversión de equilibrio X_H2,e.

NOTA de unidades: el enunciado da K_C = 325 "L/mol", pero para A + 2B -> C
(delta = -2) las unidades consistentes son L²/mol²; se usa 325 con C_i en
mol/L tal cual (supuesto reportado — el X_e esperado 0.206 confirma esta
lectura).

Esperado: C_T0 = 0.04252 mol/L, C_H2,0 = 0.02126 mol/L, eps = -0.5,
X_H2,e ~ 0.206.
"""

from exercise import run_exercise

# ══════════════════════════ INPUTS ════════════════════════════════════════
titulo = "P1 examen: hidrogenacion de o-cresol"

reacciones = ["C7H8O + 2 H2 = C7H12O"]

alimentacion = {"C7H8O": 1.0, "H2": 1.0}     # equimolar; H2 limita (pide 2)

fase = "gas"
P, P0, unidad = 2 * 1.01325, 1.0, "bar"      # 2 atm expresadas en bar
temperaturas = [573.15]                      # 300 °C

modo_K = "directo"
K = {573.15: [325.0]}                        # K_C tal como la da el enunciado
K_es_Kc = True                               # se convierte a K de fracción mol

datos_termo = None
limitante = "H2"
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

resultado = run_exercise(
    titulo=titulo, reacciones=reacciones, alimentacion=alimentacion,
    fase=fase, P=P, P0=P0, unidad=unidad, temperaturas=temperaturas,
    modo_K=modo_K, K=K, datos_termo=datos_termo, K_es_Kc=K_es_Kc,
    limitante=limitante,
)

# ═══════════ Inciso (a): tabla estequiométrica y C_i, P_i a X = 0.90 ══════
# Método B (Fogler) normalizado por mol del limitante H2 (nu_B = -1):
#   v = v0(1 + eps·X) a T,P constantes;  C_i = F_i/v;  P_i = C_i·R·T
from equilibrium import R_ATM, StoichiometricTable, format_table

T0 = 573.15
tbl = StoichiometricTable(
    nu={"H2": -1.0, "C7H8O": -0.5, "C7H12O": 0.5},   # por mol de H2
    theta={"H2": 1.0, "C7H8O": 1.0},                  # alimentación equimolar
    limiting="H2",
    T0=T0, P0=2.0, pressure_unit="atm",               # 2 atm, 300 °C
    variable_volume=True,                             # flujo a T,P constantes
)
X = 0.90
print(f"\n{'=' * 70}\n  Inciso (a): tabla estequiométrica (base H2) y "
      f"C_i, P_i a X = {X:g}\n{'=' * 70}")
print(f"  y_H2,0 = {tbl.y_A0:g}   delta' = {tbl.delta:+g} (por mol H2)   "
      f"eps = y_H2,0·delta' = {tbl.epsilon:+g}")
print(f"  C_T0 = P/(R·T) = {tbl.C_T0:.5f} mol/L   "
      f"C_H2,0 = y_H2,0·C_T0 = {tbl.C_A0:.5f} mol/L")
print(f"  v/v0 = 1 + eps·X = {1 + tbl.epsilon * X:.4f}\n")
print("  Tabla simbólica (por mol de H2):")
print("  " + format_table(tbl.symbolic_rows(),
      ["Especie", "Inicial", "Cambio", "Remanente", "Concentración"]
      ).replace("\n", "\n  "))
n, C = tbl.moles(X), tbl.concentrations(X)
filas = [[sp, n[sp], C[sp], C[sp] * R_ATM * T0] for sp in tbl.species]
filas.append(["TOTAL", sum(n.values()), sum(C.values()),
              sum(C.values()) * R_ATM * T0])
print(f"\n  Evaluada en X = {X:g}:")
print("  " + format_table(
    filas, ["Especie", "F_i (mol)", "C_i (mol/L)", "P_i (atm)"]
).replace("\n", "\n  "))

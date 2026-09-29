"""Parcial 1 — Problema 3 (4.5 pts): reformado de metano con vapor.

Balanceadas (en el examen R1 viene esquelética):
    R1: CH4 + H2O = CO + 3 H2      (delta = +2, endotérmica)
    R2: CO + H2O = CO2 + H2        (WGS, delta = 0, exotérmica)

Alimentación 1 mol CH4 + 1 mol H2O; salida en equilibrio a 700 K y
1000 K. El examen no da P: se asume P = 1 bar (= P° de las K).
K(T) de las Tablas C.4/C.1 de SVA vía ec. 13.18 (data/sva_tables.json).
"""

from exercise import run_exercise

# ══════════════════════════ INPUTS ════════════════════════════════════════
titulo = "Parcial P3: reformado de metano"

reacciones = ["CH4 + H2O = CO + 3 H2",     # R1 (balanceada con 3 H2)
              "CO + H2O = CO2 + H2"]       # R2: WGS

alimentacion = {"CH4": 1.0, "H2O": 1.0}

fase = "gas"
P, P0, unidad = 1.0, 1.0, "bar"            # supuesto: P = 1 bar (no dada)
temperaturas = [700.0, 1000.0]

modo_K = "termo"
especies_termo = ["CH4", "H2O", "CO", "H2", "CO2"]
datos_termo = None
K = None
K_es_Kc = False

limitante = "CH4"
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

resultado = run_exercise(
    titulo=titulo, reacciones=reacciones, alimentacion=alimentacion,
    fase=fase, P=P, P0=P0, unidad=unidad, temperaturas=temperaturas,
    modo_K=modo_K, K=K, datos_termo=datos_termo,
    especies_termo=especies_termo, K_es_Kc=K_es_Kc, limitante=limitante,
)

print("\n  Discusión (inciso final):")
for T in temperaturas:
    t1, t2 = resultado["resultados"][T]["termo"]
    xi1, xi2 = resultado["resultados"][T]["res"].xi
    print(f"    T = {T:g} K: dH_R1 = {t1['dH'] / 1000:+.1f} kJ/mol "
          f"(endotérmica), dH_R2 = {t2['dH'] / 1000:+.1f} kJ/mol "
          f"(exotérmica);  K1 = {t1['K']:.4g}, K2 = {t2['K']:.4g};  "
          f"xi = ({xi1:.4f}, {xi2:.4f})")
print("    R1 endotérmica -> K1 y xi_1 crecen con T (Van't Hoff); "
      "R2 exotérmica -> K2 cae con T.")

"""Sesión 9, Ej. 3 (Koretsky, Fig. 14-5): craqueo térmico de n-butano.

    C4H10 = C3H6 + CH4        (gas, delta = +1)

Reactor de estado estacionario a T = 500 K y P constante; alimentación
10 mol/s de n-butano; salida en equilibrio. Flujo de propileno para:
  A) P = 1 bar      B) P = 25 bar   (gas ideal)

K(500) se calcula de data/sva_tables.json (Tablas C.4/C.1 de SVA,
C4H10 cotejado contra el libro por el usuario) con la ec. 13.18.

Esperado (cálculo a mano): K(500) ~ 1.13;
  P = 1 bar : xi ~ 7.28 mol/s, y_C3H6 ~ 0.421, X ~ 0.728
  P = 25 bar: xi ~ 2.08 mol/s, y_C3H6 ~ 0.172, X ~ 0.208
A mayor P baja la conversión porque delta = +1 (Le Chatelier).
"""

from exercise import run_exercise

# ══════════════════════════ INPUTS ════════════════════════════════════════
titulo_base = "Craqueo de n-butano (Koretsky 14-5)"

reacciones = ["C4H10 = C3H6 + CH4"]
alimentacion = {"C4H10": 10.0}       # mol/s

fase = "gas"
P0, unidad = 1.0, "bar"
presiones = [1.0, 25.0]              # los dos incisos del enunciado
temperaturas = [500.0]

modo_K = "termo"
especies_termo = ["C4H10", "C3H6", "CH4"]   # se cargan de data/sva_tables.json
datos_termo = None

limitante = "C4H10"

esperado = {1.0: {"xi": 7.28, "y_C3H6": 0.421, "X": 0.728},
            25.0: {"xi": 2.08, "y_C3H6": 0.172, "X": 0.208}}
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

for P in presiones:
    r = run_exercise(
        titulo=f"{titulo_base} - P = {P:g} bar",
        reacciones=reacciones, alimentacion=alimentacion, fase=fase,
        P=P, P0=P0, unidad=unidad, temperaturas=temperaturas,
        modo_K=modo_K, datos_termo=datos_termo,
        especies_termo=especies_termo, limitante=limitante,
    )

    # ---- la ecuación que resolvió el solver, con números, para el examen ----
    caso = r["resultados"][500.0]
    K = caso["K"][0]
    xi = float(caso["res"].xi[0])
    y = caso["res"].mole_fractions
    delta = 1
    lhs = y["C3H6"] * y["CH4"] / y["C4H10"]
    rhs = K * (P0 / P) ** delta
    print(f"\n  Ecuación de equilibrio resuelta (P = {P:g} bar), "
          f"para reproducir a mano:")
    print("    y_C3H6·y_CH4 / y_C4H10 = K·(P°/P)^δ")
    print("    en función de ξ (n_T = 10 + ξ):  "
          "[ξ/(10+ξ)]·[ξ/(10+ξ)] / [(10−ξ)/(10+ξ)] = ξ²/(100 − ξ²)")
    print(f"    lado izq.:  {y['C3H6']:.5f}·{y['CH4']:.5f} / {y['C4H10']:.5f}"
          f" = {lhs:.5f}")
    print(f"    lado der.:  {K:.5f}·(1/{P:g})^{delta} = {rhs:.5f}"
          f"      (|dif| = {abs(lhs - rhs):.1e})")
    print(f"    => ξ = 10·sqrt(c/(1+c)) con c = K·(P°/P)^δ = {rhs:.5f}"
          f"  ->  ξ = {10.0 * (rhs / (1 + rhs)) ** 0.5:.4f} mol/s")

    e = esperado[P]
    print(f"  Cotejo vs cálculo a mano:  ξ = {xi:.4f} (≈{e['xi']:g})   "
          f"F_C3H6 = {xi:.4f} mol/s   y_C3H6 = {y['C3H6']:.4f} "
          f"(≈{e['y_C3H6']:g})   X = {caso['X']:.4f} (≈{e['X']:g})")

print("\n  Conclusión (Le Chatelier): delta = +1, así que subir P de 1 a "
      "25 bar empuja el equilibrio hacia el reactivo:\n  el flujo de "
      "propileno cae de ~7.3 a ~2.1 mol/s.")

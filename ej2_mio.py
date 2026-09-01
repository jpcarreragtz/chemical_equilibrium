"""Ejercicio 2 (sesión 6): Reformado de propano con vapor — método de avance de reacción.

Reacciones (como en mis notas, pág. 91 del PDF):
    R1: C3H8 + 3 H2O <-> 3 CO  + 7 H2      delta_1 = +6
    R2: CO   +   H2O <-> CO2   +   H2      delta_2 =  0
    R3: C3H8 + 2 H2  <-> 3 CH4             delta_3 =  0

Balances de moles en función de los avances xi_1, xi_2, xi_3:
    n_C3H8 = 1   - xi_1        - xi_3
    n_H2O  = 4   - 3 xi_1 - xi_2
    n_CO   =       3 xi_1 - xi_2
    n_H2   =       7 xi_1 + xi_2 - 2 xi_3
    n_CO2  =               xi_2
    n_CH4  =                      3 xi_3
    n_N2   = 0.5  (inerte)
    n_T    = 5.5 + 6 xi_1

Alimentación: 1 mol C3H8, 4 mol H2O, 0.5 mol N2. P = 1 bar.

Rendimientos de H2:
    Y_H2 (teórico)      = n_H2 / 10   (10 mol H2 por mol C3H8 si R1 y R2 completas)
    Y_H2 (alimentación) = n_H2 / 8    (máximo con 4 mol H2O disponibles)
"""

from equilibrium import solve_extents

# ---------- Especies (el orden define las columnas de nu) ----------
especies = ["C3H8", "H2O", "CO", "H2", "CO2", "CH4", "N2"]

# ---------- Moles iniciales ----------
n0 = {"C3H8": 1.0, "H2O": 4.0, "CO": 0.0, "H2": 0.0,
      "CO2": 0.0, "CH4": 0.0, "N2": 0.5}

# ---------- Matriz estequiométrica nu (filas = R1, R2, R3) ----------
#          C3H8  H2O   CO   H2  CO2  CH4   N2
nu = [
    [  -1,  -3,   3,   7,   0,   0,   0 ],   # R1: reformado
    [   0,  -1,  -1,   1,   1,   0,   0 ],   # R2: WGS
    [  -1,   0,   0,  -2,   0,   3,   0 ],   # R3: metanación
]

# ---------- Composición atómica (checker de balance) ----------
atomos = {
    "C3H8": {"C": 3, "H": 8},
    "H2O":  {"H": 2, "O": 1},
    "CO":   {"C": 1, "O": 1},
    "H2":   {"H": 2},
    "CO2":  {"C": 1, "O": 2},
    "CH4":  {"C": 1, "H": 4},
    "N2":   {"N": 2},
}

# ---------- Constantes de equilibrio por temperatura (orden K1, K2, K3) ----------
casos = {
    700:  [7.2285e-3, 7.5034,  2.4203e10],
    900:  [1.3407e6,  1.5591,  2.4864e8 ],
    1000: [1.0501e9,  0.89956, 5.0080e7 ],
}
P = 1.0   # bar

# ---------- (OPCIONAL) Mis valores resueltos a mano, para comparar ----------
# Llena con tus xi calculados a mano (o deja None). Ej: 700: (0.24, 0.69, 0.76)
xi_a_mano = {
    700:  (None, None, None),
    900:  (None, None, None),
    1000: (None, None, None),
}

H2_MAX_TEORICO = 10.0   # mol H2 por mol C3H8 si R1 y R2 son completas

# Máximo permitido por la alimentación (1 mol C3H8, 4 mol H2O):
#   R1 completa: 1 C3H8 + 3 H2O -> 3 CO + 7 H2   (quedan 1 mol H2O, 3 mol CO)
#   R2 con el agua sobrante: 1 CO + 1 H2O -> 1 CO2 + 1 H2
#   Total = 7 + 1 = 8 mol H2; propano y agua se agotan a la vez.
# (NO es 10/6*4 = 6.667: eso reparte el agua proporcionalmente entre
#  R1 y R2, pero el reformado rinde 7/3 H2 por H2O y la WGS solo 1,
#  así que el máximo se logra reformando primero y usando el sobrante en WGS.)
H2_MAX_ALIMENTACION = 8.0


def balances(xi1, xi2, xi3):
    """Balances de moles — exactamente los de mis notas."""
    n = {
        "C3H8": 1 - xi1 - xi3,
        "H2O":  4 - 3*xi1 - xi2,
        "CO":   3*xi1 - xi2,
        "H2":   7*xi1 + xi2 - 2*xi3,
        "CO2":  xi2,
        "CH4":  3*xi3,
        "N2":   0.5,
    }
    n_T = 5.5 + 6*xi1
    return n, n_T


resumen = []
xi_resueltos = {}   # xi de cada T, para la verificación final contra la clave
for T, K in casos.items():
    print(f"\n{'='*64}\n  T = {T} K\n{'='*64}")

    # ---------- Resolver (Método A) ----------
    res = solve_extents(especies, n0, nu, K, P, atoms=atomos)  # <<< CLAUDE: usa la firma real
    xi1, xi2, xi3 = res.xi                                     # <<< CLAUDE: extrae los xi como los devuelva el módulo
    xi_resueltos[T] = (xi1, xi2, xi3)

    # ---------- a) Avances de reacción ----------
    print("  a) Avances de reacción:")
    print(f"     xi_1 (reformado)  = {xi1:.5f} mol")
    print(f"     xi_2 (WGS)        = {xi2:.5f} mol")
    print(f"     xi_3 (metanación) = {xi3:.5f} mol")

    if all(v is not None for v in xi_a_mano[T]):
        m1, m2, m3 = xi_a_mano[T]
        print(f"     a mano:  xi_1 = {m1:.4f}  xi_2 = {m2:.4f}  xi_3 = {m3:.4f}   "
              f"(dif: {abs(xi1-m1):.4f}, {abs(xi2-m2):.4f}, {abs(xi3-m3):.4f})")

    # ---------- b) Moles y fracciones mol ----------
    n, n_T = balances(xi1, xi2, xi3)
    print(f"\n  b) Moles y fracciones mol   (n_T = 5.5 + 6 xi_1 = {n_T:.5f} mol)")
    print(f"     {'Especie':7s} {'balance':>26s} {'n_i (mol)':>12s} {'y_i':>10s}")
    expr = {
        "C3H8": "1 - xi1 - xi3", "H2O": "4 - 3 xi1 - xi2", "CO": "3 xi1 - xi2",
        "H2": "7 xi1 + xi2 - 2 xi3", "CO2": "xi2", "CH4": "3 xi3", "N2": "0.5",
    }
    for esp in especies:
        print(f"     {esp:7s} {expr[esp]:>26s} {n[esp]:12.5g} {n[esp]/n_T:10.5f}")
    print(f"     {'':7s} {'':>26s} {'sum y_i =':>12s} {sum(n.values())/n_T:10.6f}")

    # ---------- c) Conversión del propano ----------
    X = (n0["C3H8"] - n["C3H8"]) / n0["C3H8"]     # = xi1 + xi3
    S_ref = xi1 / (xi1 + xi3)                     # fracción vía reformado
    print(f"\n  c) Conversión del propano  X = (1 - n_C3H8)/1 = xi1 + xi3 = {X:.5f}")
    print(f"     Fracción vía reformado  S_reformado = xi1/(xi1+xi3) = {S_ref:.5f}")

    # ---------- d) Selectividad respecto a TODOS los productos con carbono ----------
    C_total = n["CO"] + n["CO2"] + n["CH4"]
    S_CO  = n["CO"]  / C_total
    S_CO2 = n["CO2"] / C_total
    S_CH4 = n["CH4"] / C_total
    print(f"\n  d) Selectividades  S_i = n_i / (n_CO + n_CO2 + n_CH4),  "
          f"denominador = {C_total:.5f} mol")
    print(f"     S_CO = {S_CO:.5f}   S_CO2 = {S_CO2:.5f}   S_CH4 = {S_CH4:.5f}   "
          f"(suma = {S_CO+S_CO2+S_CH4:.5f})")

    # ---------- e) Rendimiento de H2 y razones del gas de síntesis ----------
    Y_H2   = n["H2"] / H2_MAX_TEORICO
    Y_H2_alim = n["H2"] / H2_MAX_ALIMENTACION
    H2_CO  = n["H2"] / n["CO"]
    H2_COx = n["H2"] / (n["CO"] + n["CO2"])
    print(f"\n  e) Rendimiento  Y_H2 = n_H2/10 = {Y_H2:.5f}")
    print(f"     Y_H2 (máx. por alimentación) = n_H2/8 = {Y_H2_alim:.5f}")
    print(f"     H2/CO = {H2_CO:.4f}    H2/COx = n_H2/(n_CO+n_CO2) = {H2_COx:.4f}")

    resumen.append((T, X, S_ref, S_CO, S_CO2, S_CH4, Y_H2, H2_CO, H2_COx, Y_H2_alim))

# ---------- Tabla final (mismo formato que la clave de respuestas) ----------
print(f"\n{'='*64}\n  Conversión, selectividades y rendimiento\n{'='*64}")
print(f"  {'T_K':>5s} {'X_C3H8':>7s} {'S_reformado':>12s} {'S_CO':>9s} {'S_CO2':>9s} "
      f"{'S_CH4':>9s} {'Y_H2':>9s} {'Y_H2_alim':>9s} {'H2_CO':>8s} {'H2_COx':>8s}")
for T, X, S_ref, S_CO, S_CO2, S_CH4, Y_H2, H2_CO, H2_COx, Y_H2_alim in resumen:
    print(f"  {T:5d} {X:7.4g} {S_ref:12.5f} {S_CO:9.5f} {S_CO2:9.5f} "
          f"{S_CH4:9.5f} {Y_H2:9.5f} {Y_H2_alim:9.5f} {H2_CO:8.4f} {H2_COx:8.4f}")

# ---------- Verificación contra la clave de respuestas del profesor ----------
TOL_CLAVE = 1e-4   # tolerancia relativa por valor

clave_xi = {   # T: (xi_1, xi_2, xi_3)
    700:  (0.24044, 0.69089, 0.75956),
    900:  (0.51886, 0.69625, 0.48114),
    1000: (0.75718, 0.41774, 0.24282),
}
clave_metricas = {   # T: (S_CO, S_CO2, S_CH4, Y_H2, H2_CO, H2_COx)
    700:  (0.010139, 0.2303,  0.75956, 0.085482, 28.105, 1.1851),
    900:  (0.28677,  0.23208, 0.48114, 0.33659,  3.9124, 2.1624),
    1000: (0.61794,  0.13925, 0.24282, 0.52324,  2.8225, 2.3034),
}

print(f"\n{'='*64}\n  Verificación contra la clave del profesor "
      f"(tol. rel. 1e-4)\n{'='*64}")
todo_pass = True
for T, X, S_ref, S_CO, S_CO2, S_CH4, Y_H2, H2_CO, H2_COx, Y_H2_alim in resumen:
    print(f"  T = {T} K")
    comparaciones = list(zip(("xi_1", "xi_2", "xi_3"),
                             xi_resueltos[T], clave_xi[T]))
    comparaciones += list(zip(("S_CO", "S_CO2", "S_CH4", "Y_H2", "H2_CO", "H2_COx"),
                              (S_CO, S_CO2, S_CH4, Y_H2, H2_CO, H2_COx),
                              clave_metricas[T]))
    for nombre, calc, ref in comparaciones:
        dif = abs(calc - ref)
        ok = dif / abs(ref) < TOL_CLAVE
        todo_pass &= ok
        print(f"    {nombre:6s}  calc = {calc:9.5f}   clave = {ref:9.5f}   "
              f"|dif| = {dif:.1e}   [{'PASS' if ok else 'FAIL'}]")
print(f"\n  Resultado: "
      f"{'TODOS los valores PASS' if todo_pass else '*** HAY VALORES FAIL ***'}")

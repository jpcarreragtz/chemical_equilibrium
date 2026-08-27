"""Ejercicio 2 - Sesión 6: Reformado de propano con vapor.
Alimentación: 1 mol C3H8, 4 mol H2O, 0.5 mol N2 (inerte). P = 1 bar.
Método A (avance de reacción), 3 reacciones simultáneas."""

from equilibrium import (conversion, format_table, run_checks_method_a,
                         selectivity, solve_extents, yield_fraction)

# ---------- Paso 1: especies (el ORDEN aquí manda en todo lo demás) ----------
especies = ["C3H8", "H2O", "CO", "H2", "CO2", "CH4", "N2"]

# ---------- Paso 2: moles iniciales ----------
n0 = {"C3H8": 1.0, "H2O": 4.0, "CO": 0.0, "H2": 0.0,
      "CO2": 0.0, "CH4": 0.0, "N2": 0.5}

# ---------- Paso 3: matriz estequiométrica ----------
# Filas = reacciones R1, R2, R3. Columnas = especies EN EL ORDEN de la lista.
# Reactivos negativos, productos positivos, inerte = 0 en todas.
#        C3H8  H2O   CO   H2  CO2  CH4   N2
nu = [
    [  -1,  -3,   3,   7,   0,   0,   0 ],   # R1: reformado
    [   0,  -1,  -1,   1,   1,   0,   0 ],   # R2: WGS
    [  -1,   0,   0,  -2,   0,   3,   0 ],   # R3: metanación
]

# ---------- Paso 4: composición atómica (red de seguridad) ----------
atomos = {
    "C3H8": {"C": 3, "H": 8},
    "H2O":  {"H": 2, "O": 1},
    "CO":   {"C": 1, "O": 1},
    "H2":   {"H": 2},
    "CO2":  {"C": 1, "O": 2},
    "CH4":  {"C": 1, "H": 4},
    "N2":   {"N": 2},
}

# ---------- Paso 5: K's por temperatura (mismo orden que las filas de nu) ----------
casos = {
    700:  [7.2285e-3, 7.5034,  2.4203e10],
    900:  [1.3407e6,  1.5591,  2.4864e8 ],
    1000: [1.0501e9,  0.89956, 5.0080e7 ],
}
P = 1.0   # bar

# H2 máximo teórico: R1+R2 completas equivalen a C3H8 + 6 H2O -> 3 CO2 + 10 H2,
# o sea 10 mol H2 por mol de C3H8; pero esa ruta pide 6 mol H2O por mol de
# propano y solo se alimentan 4, así que el agua limita a (10/6)*4 = 6.667 mol.
n_A0 = n0["C3H8"]
H2_MAX_PROPANO = 10.0 * n_A0
H2_MAX_AGUA = 10.0 / 6.0 * n0["H2O"]

# ---------- Resolver para cada T ----------
resultados = {}   # T -> ExtentResult, para las tablas comparativas del final
for T, K in casos.items():
    print(f"\n{'='*50}\n  T = {T} K\n{'='*50}")
    res = solve_extents(especies, n0, nu, K, P=P, P0=1.0, atoms=atomos)
    resultados[T] = res

    print("  Avances de reacción:  "
          + "  ".join(f"xi_{j + 1} = {x:.5f} mol"
                      for j, x in enumerate(res.xi)))
    print("  " + format_table(
        [[sp, res.moles[sp], res.mole_fractions[sp]] for sp in especies],
        ["Especie", "n_i (mol)", "y_i"]).replace("\n", "\n  "))
    run_checks_method_a(res, n0, atoms=atomos)

# ---------- Salidas del ejercicio ----------
print(f"\n{'='*50}\n  RESULTADOS DEL EJERCICIO\n{'='*50}")
print(f"  H2 máximo teórico: {H2_MAX_PROPANO:g} mol (10 por mol de C3H8, "
      f"R1+R2 completas)")
print(f"  H2 máximo por el agua alimentada: {H2_MAX_AGUA:.4g} mol "
      f"(= 10/6 * 4 mol H2O)")

filas = []
for T, res in sorted(resultados.items()):
    n = res.moles
    xi1, xi2, xi3 = res.xi

    # a) conversión de propano: se consume en R1 y R3 -> X = (xi1+xi3)/n_A0
    X = (xi1 + xi3) / n_A0
    assert abs(X - conversion(n_A0, n["C3H8"])) < 1e-9   # misma X por balance

    # b) selectividades entre productos con carbono (cocientes de moles)
    S_CO_CH4 = selectivity(n["CO"], n["CH4"])
    S_CO_CO2 = selectivity(n["CO"], n["CO2"])
    S_CH4_CO2 = selectivity(n["CH4"], n["CO2"])

    # c) rendimiento de H2 de las dos formas
    Y_H2_teo = yield_fraction(n["H2"], n_A0, 10.0)   # n_H2 / 10
    Y_H2_agua = n["H2"] / H2_MAX_AGUA                # n_H2 / 6.667

    print(f"\n  T = {T} K:")
    print(f"    a) X_C3H8 = xi_1 + xi_3 = {xi1:.5f} + {xi3:.5f} = {X:.5f}")
    print(f"    b) S CO/CH4 = {S_CO_CH4:.5g}   S CO/CO2 = {S_CO_CO2:.5g}   "
          f"S CH4/CO2 = {S_CH4_CO2:.5g}")
    print(f"    c) rendimiento H2: {Y_H2_teo:.4f} del máximo teórico (n_H2/10)"
          f" | {Y_H2_agua:.4f} del máximo por agua (n_H2/{H2_MAX_AGUA:.4g})")

    filas.append([T, xi1, xi2, xi3, X, n["H2"], S_CO_CH4, S_CO_CO2,
                  S_CH4_CO2, Y_H2_teo, Y_H2_agua])

# d) tabla comparativa entre las tres temperaturas
print("\n  d) Comparación entre temperaturas:")
print("  " + format_table(
    filas,
    ["T (K)", "xi_1", "xi_2", "xi_3", "X_C3H8", "n_H2 (mol)",
     "S CO/CH4", "S CO/CO2", "S CH4/CO2", "Y_H2 (n/10)",
     "Y_H2 (agua)"]).replace("\n", "\n  "))

print("\n  Fracciones molares y_i en el equilibrio:")
temps = sorted(resultados)
print("  " + format_table(
    [[sp] + [resultados[T].mole_fractions[sp] for T in temps]
     for sp in especies],
    ["Especie"] + [f"T = {T} K" for T in temps]).replace("\n", "\n  "))

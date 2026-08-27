"""Ejercicio 2 (mi versión): Reformado de propano con vapor.
R1: C3H8 + 3 H2O <-> 3 CO + 7 H2
R2: CO + H2O <-> CO2 + H2
R3: C3H8 + 2 H2 <-> 3 CH4
Alimentación: 1 mol C3H8, 4 mol H2O, 0.5 mol N2 (inerte). P = 1 bar.
"""

import inspect
from equilibrium import solve_extents

# ---------- Especies (el orden manda en las columnas de nu) ----------
especies = ["C3H8", "H2O", "CO", "H2", "CO2", "CH4", "N2"]

# ---------- Moles iniciales ----------
n0 = {"C3H8": 1.0, "H2O": 4.0, "CO": 0.0, "H2": 0.0,
      "CO2": 0.0, "CH4": 0.0, "N2": 0.5}

# ---------- Matriz estequiométrica (filas=R1,R2,R3; reactivos negativos) ----------
#          C3H8  H2O   CO   H2  CO2  CH4   N2
nu = [
    [  -1,  -3,   3,   7,   0,   0,   0 ],   # R1: reformado
    [   0,  -1,  -1,   1,   1,   0,   0 ],   # R2: WGS
    [  -1,   0,   0,  -2,   0,   3,   0 ],   # R3: metanación
]

# ---------- Composición atómica (para el checker de balance) ----------
atomos = {
    "C3H8": {"C": 3, "H": 8},
    "H2O":  {"H": 2, "O": 1},
    "CO":   {"C": 1, "O": 1},
    "H2":   {"H": 2},
    "CO2":  {"C": 1, "O": 2},
    "CH4":  {"C": 1, "H": 4},
    "N2":   {"N": 2},
}

# ---------- K's por temperatura (orden = filas de nu: K1, K2, K3) ----------
casos = {
    700:  [7.2285e-3, 7.5034,  2.4203e10],
    900:  [1.3407e6,  1.5591,  2.4864e8 ],
    1000: [1.0501e9,  0.89956, 5.0080e7 ],
}
P = 1.0   # bar

# ---------- Máximos de H2 para el inciso e ----------
H2_MAX_TEORICO = 10.0            # 10 mol H2 por mol C3H8 (R1 + R2 completas)
H2_MAX_AGUA = 10.0 / 6.0 * 4.0   # 6.667 mol: límite por alimentar solo 4 mol H2O


# ---------- Llamada adaptativa: lee la firma real de solve_extents ----------
def llamar_solver(K_lista):
    """Empata nuestros datos con los nombres de parámetros reales del módulo."""
    valores = {
        "species": especies, "especies": especies, "names": especies,
        "n0": n0, "moles0": n0, "initial_moles": n0, "n_initial": n0,
        "nu": nu, "stoich": nu, "stoich_matrix": nu, "nu_matrix": nu,
        "K": K_lista, "Ks": K_lista, "K_values": K_lista, "Keq": K_lista,
        "P": P, "pressure": P,
        "atoms": atomos, "atomos": atomos, "atomic": atomos,
        "atomic_composition": atomos, "composition": atomos, "elements": atomos,
    }
    firma = inspect.signature(solve_extents)
    kwargs = {}
    for nombre, par in firma.parameters.items():
        if nombre in valores:
            kwargs[nombre] = valores[nombre]
        elif par.default is inspect.Parameter.empty:
            raise TypeError(
                f"No sé qué pasar al parámetro requerido '{nombre}'. "
                f"Firma real: {firma}. Agrega '{nombre}' al dict 'valores'."
            )
    return solve_extents(**kwargs)


def extraer_xi(res):
    """Saca (xi1, xi2, xi3) sin importar si res es dict, objeto o arreglo."""
    if isinstance(res, dict):
        for llave in ("xi", "extents", "xis", "extent"):
            if llave in res:
                return tuple(float(v) for v in res[llave])
        raise KeyError(f"No encontré los xi en el dict. Llaves: {list(res)}")
    for attr in ("xi", "extents", "x"):
        if hasattr(res, attr):
            return tuple(float(v) for v in getattr(res, attr))
    return tuple(float(v) for v in res)   # último recurso: res ES el arreglo


# ---------- Resolver para cada temperatura ----------
resumen = []
for T, K in casos.items():
    print(f"\n{'='*50}\n  T = {T} K\n{'='*50}")

    res = llamar_solver(K)
    xi1, xi2, xi3 = extraer_xi(res)

    # Balances de moles (los de las anotaciones de la diapositiva)
    n_C3H8 = 1 - xi1 - xi3
    n_H2O  = 4 - 3*xi1 - xi2
    n_CO   = 3*xi1 - xi2
    n_H2   = 7*xi1 + xi2 - 2*xi3
    n_CO2  = xi2
    n_CH4  = 3*xi3
    n_N2   = 0.5
    n_T = 5.5 + 6*xi1

    print(f"  xi_1 = {xi1:.5f}   xi_2 = {xi2:.5f}   xi_3 = {xi3:.5f}")
    print(f"  {'Especie':8s} {'n_i (mol)':>12s} {'y_i':>10s}")
    tabla = [("C3H8", n_C3H8), ("H2O", n_H2O), ("CO", n_CO), ("H2", n_H2),
             ("CO2", n_CO2), ("CH4", n_CH4), ("N2", n_N2)]
    for nombre, n_i in tabla:
        print(f"  {nombre:8s} {n_i:12.5g} {n_i/n_T:10.5f}")

    # Checks básicos hechos aquí mismo (independientes del módulo)
    suma_y = sum(n_i for _, n_i in tabla) / n_T
    negativos = [nom for nom, n_i in tabla if n_i < -1e-9]
    print(f"  Checks: sum(y_i) = {suma_y:.10f} | "
          f"{'todos n_i >= 0' if not negativos else 'NEGATIVOS: ' + str(negativos)}")

    # ---------- Incisos c, d, e ----------
    X = xi1 + xi3                                  # c) conversión del propano
    S_CO_CH4  = n_CO  / n_CH4                      # d) selectividades
    S_CO_CO2  = n_CO  / n_CO2
    S_CH4_CO2 = n_CH4 / n_CO2
    Y_teorico = n_H2 / H2_MAX_TEORICO              # e) rendimientos
    Y_agua    = n_H2 / H2_MAX_AGUA

    print(f"\n  c) X_C3H8 = xi1 + xi3 = {X:.5f}")
    print(f"  d) S CO/CH4 = {S_CO_CH4:.4f}   S CO/CO2 = {S_CO_CO2:.4f}   "
          f"S CH4/CO2 = {S_CH4_CO2:.4f}")
    print(f"  e) Y_H2 teorico = {Y_teorico:.4f}   Y_H2 por agua = {Y_agua:.4f}")

    resumen.append((T, xi1, xi2, xi3, X, n_H2, S_CO_CH4, Y_teorico, Y_agua))

# ---------- Tabla comparativa final ----------
print(f"\n{'='*50}\n  COMPARACIÓN ENTRE TEMPERATURAS\n{'='*50}")
print(f"  {'T(K)':>5s} {'xi1':>8s} {'xi2':>8s} {'xi3':>8s} {'X':>6s} "
      f"{'n_H2':>8s} {'S CO/CH4':>9s} {'Y_teo':>7s} {'Y_agua':>7s}")
for T, xi1, xi2, xi3, X, nH2, S1, Yt, Ya in resumen:
    print(f"  {T:5d} {xi1:8.4f} {xi2:8.4f} {xi3:8.4f} {X:6.3f} "
          f"{nH2:8.4f} {S1:9.4f} {Yt:7.4f} {Ya:7.4f}")
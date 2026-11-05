"""Sesiones 23–24 — Si epitaxial: Si(s) depositado en el equilibrio desde
SiCl4 o SiCl3H con H2, a 1300 K y 1 bar.

Especies: Si(s), SiCl2, SiCl4, SiCl3H, SiCl2H2, SiClH3, SiH4, H2, HCl
(en las matrices Si(s) se escribe "Si"). Elementos Si, Cl, H. Gas ideal;
Si(s) puro con actividad 1. Base: 1 mol de clorosilano + R mol de H2,
R = 1 ... 150.

Componentes con el criterio del proyecto (CLAUDE.md, elegir_componentes):
especies de alimentación en el orden dado, SiCl4, SiCl3H, H2.

Método:
  1. Matrices (core/independencia.py): A, Aᵀ, rango, r, [I | Z], F → −F
     para [−I | −Z], y ν completa (−Z sobre la identidad, una columna por
     reacción); las reacciones se leen de los signos de ν. Cada columna
     forma UNA no componente a partir de SiCl4, SiCl3H y H2:
         R1: 4 SiCl3H -> Si + 3 SiCl4 + 2 H2         (único con Si(s))
         R2–R5: SiCl2, SiCl2H2, SiClH3, SiH4 desde SiCl3H/SiCl4
         R6: SiCl4 + H2 -> SiCl3H + HCl
  2. Con u = ln y de las componentes y K_j = Π a_i^ν (gas: a_i = y_i·P/P°):
         ln a_k = ln K_j − δ_j ln(P/P°) − Σ_c ν_cj u_c   (k = no componente)
     Con Si(s) (R1 a R6, a_Si = 1). Incógnitas u_SiCl4, u_SiCl3H, u_H2;
     residuales logarítmicos
         f1 = ln Σy_i = 0
         f2 = ln(Σy_i·Cl_i / Σy_i·H_i) − ln(Cl_0/H_0) = 0
         f3 = ln a_Si (de R1) = 0
     n_T = Cl_0 / Σy_i·Cl_i ,   n_Si(s) = 1 − n_T·Σy_i·Si_i
     Si n_Si(s) < 0 no hay sólido: régimen SOLO GAS con R2 a R6, y f3 se
     cambia por el balance
         f3 = ln(Σy_i·Si_i / Σy_i·Cl_i) − ln(Si_0/Cl_0) = 0
     Se reporta n_Si(s) = 0, y R1 da a_Si, que debe ser < 1.
  3. Verificación: minimización directa de G (SLSQP sobre n_i con
     balances de elementos lineales), balances, Σy = 1, las K de esta base
     y de la base anterior (Si, HCl, H2), y que SiCl4 con R = 10 y SiCl3H
     con R = 7 (mismo Cl:H) den el mismo gas.

core/equilibrium.solve_extents trata todas las especies como gas (todas
entran en Σy) y no maneja un sólido puro ni el cambio de régimen, por eso
ambos se resuelven aquí.

Datos: ΔfG°(1300 K) de NIST-JANAF (kJ/mol); NO se usa data/sva_tables.json.
R = 8.314462618 J/(mol·K) (CODATA).
"""

# %% Setup
import csv
import pathlib
import sys
from math import exp, log

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto

import numpy as np
from scipy.optimize import brentq, least_squares, minimize

from core import SALIDAS
from core.independencia import (comprobar_conservacion, elegir_componentes,
                                escalar_columnas, formatear_matriz,
                                generan_mismo_espacio, grados_libertad,
                                matriz_atomica, matriz_estequiometrica,
                                negar_renglones, nu_completa, nu_por_especie,
                                num_reacciones_independientes, rango,
                                reaccion_desde_columna, reaccion_texto,
                                reducir_por_componentes, transpuesta)

GUARDAR = __name__ == "__main__"     # importado desde tests/: no escribe archivos

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
ESPECIES = ["Si", "SiCl2", "SiCl4", "SiCl3H", "SiCl2H2", "SiClH3", "SiH4",
            "H2", "HCl"]                       # "Si" = Si(s), sólido puro
ELEMENTOS = ["Si", "Cl", "H"]
ALIMENTACION = ["SiCl4", "SiCl3H", "H2"]       # orden para elegir componentes
COMPLEMENTO = None                             # None = demás especies en orden
COMPONENTES_ENUNCIADO = None                   # el enunciado no las da
T = 1300.0                                     # K
P, P0 = 1.0, 1.0                               # bar
R_GAS = 8.314462618                            # J/(mol·K), CODATA
CLOROSILANOS = ["SiCl4", "SiCl3H"]             # 1 mol de clorosilano + R mol H2
RELACIONES = [1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 75, 100, 125, 150]

# ΔfG°(1300 K), kJ/mol — NIST-JANAF Thermochemical Tables
DGF = {"SiCl2": -216.012, "SiCl4": -492.536, "SiCl3H": -356.537,
       "SiCl2H2": -199.368, "SiClH3": -28.482, "SiH4": 151.897,
       "HCl": -102.644, "H2": 0.0, "Si": 0.0}

GAS = [sp for sp in ESPECIES if sp != "Si"]


def titulo(t: str) -> None:
    print(f"\n{'=' * 72}\n{t}\n{'=' * 72}")


# %% 1. Matriz atómica A y su transpuesta
A, elementos = matriz_atomica(ESPECIES, ELEMENTOS)
At = transpuesta(A)
N = len(ESPECIES)

titulo("1. Matriz atómica A (especies × elementos); Si = Si(s)")
print(formatear_matriz(A, filas=ESPECIES, columnas=elementos))
print("\nAᵀ (elementos × especies)")
print(formatear_matriz(At, filas=elementos, columnas=ESPECIES))

# %% 2. Rango y número de reacciones independientes
rg = rango(At)
r = num_reacciones_independientes(A)

titulo("2. Rango y reacciones independientes")
print(f"N = {N} especies, rango(A) = rango(Aᵀ) = {rg}")
print(f"r = N − rango = {N} − {rg} = {r}")

# %% 3. Elección de componentes y reducción de Aᵀ a [I | Z]
COMPONENTES, CRITERIO = elegir_componentes(ESPECIES, A, ALIMENTACION,
                                           COMPLEMENTO, COMPONENTES_ENUNCIADO)
red = reducir_por_componentes(At, ESPECIES, COMPONENTES)
k = len(COMPONENTES)
filas_f = [f"f{i + 1}" for i in range(len(At))]

titulo("3. Componentes y operaciones de renglón sobre Aᵀ")
print(CRITERIO)
print(f"\nComponentes (pivotes): {', '.join(COMPONENTES)}; "
      f"no componentes: {', '.join(red['no_componentes'])}")
print("\nAᵀ con columnas reordenadas [componentes | no componentes]:")
print(formatear_matriz(red["inicial"], filas=filas_f, columnas=red["columnas"],
                       corte=k))
for n_paso, (desc, M) in enumerate(red["pasos"], 1):
    print(f"\nPaso {n_paso}: {desc}")
    print(formatear_matriz(M, filas=filas_f, columnas=red["columnas"], corte=k))

print("\nMatriz reducida [I | Z]:")
print(formatear_matriz(red["reducida"], filas=COMPONENTES,
                       columnas=red["columnas"], corte=k))

# %% 4. F → −F: [−I | −Z]
menos_IZ = negar_renglones(red["reducida"])

titulo("4. F → −F en todos los renglones: [−I | −Z]")
print("−Z contiene directamente los coeficientes ν de las componentes.")
print(formatear_matriz(menos_IZ, filas=COMPONENTES, columnas=red["columnas"],
                       corte=k))

# %% 5. Matriz ν completa (una columna por reacción)
nu_frac, FILAS_NU = nu_completa(red)
NU, factores = escalar_columnas(nu_frac)
nombres_rx = [f"R{i + 1}" for i in range(len(NU[0]))]
nu = nu_por_especie(NU, FILAS_NU, ESPECIES)     # una lista por reacción (K, solver)

titulo("5. Matriz ν completa: −Z sobre la identidad de las no componentes")
if any(f != 1 for f in factores):
    print(formatear_matriz(nu_frac, filas=FILAS_NU, columnas=nombres_rx))
    print("\nCada columna se escala al mínimo entero positivo: "
          + ", ".join(f"{n} ×{f}" for n, f in zip(nombres_rx, factores)) + "\n")
print(formatear_matriz(NU, filas=FILAS_NU, columnas=nombres_rx))

# %% 6. Reacciones leídas de ν
REACCIONES = [reaccion_desde_columna(col, FILAS_NU, ESPECIES)
              for col in transpuesta(NU)]

# base anterior (componentes Si, HCl, H2): solo para verificar
COMPONENTES_ANT = ["Si", "HCl", "H2"]
nu_ant = matriz_estequiometrica(reducir_por_componentes(At, ESPECIES,
                                                        COMPONENTES_ANT), ESPECIES)
MISMO_ESPACIO = generan_mismo_espacio(
    [dict(zip(ESPECIES, f)) for f in nu], [dict(zip(ESPECIES, f)) for f in nu_ant])

titulo("6. Reacciones leídas de ν (ν < 0 reactivo, ν > 0 producto)")
for nombre, txt in zip(nombres_rx, REACCIONES):
    print(f"{nombre}: {txt}")
print(f"\nMismo espacio que la base anterior ({', '.join(COMPONENTES_ANT)}): "
      f"{'SÍ' if MISMO_ESPACIO else 'NO'}")

# %% 7. Comprobación Aᵀν = 0
prod, conserva = comprobar_conservacion(red["inicial"], transpuesta(NU))

titulo("7. Comprobación Aᵀν = 0 (columnas de Aᵀ en el orden de los renglones de ν)")
print(formatear_matriz(prod, filas=elementos, columnas=nombres_rx))
print(f"\nAᵀν = 0: {'SÍ' if conserva else 'NO'} — cada reacción conserva Si, Cl y H")

# %% 8. Regla de las fases: con y sin Si(s)
A_gas, _ = matriz_atomica(GAS, ELEMENTOS)
r_gas = num_reacciones_independientes(A_gas)
FASES = {"con Si(s)": dict(pi=2, N=N, r=r, s=1),
         "sin Si(s)": dict(pi=1, N=len(GAS), r=r_gas, s=2)}
for d in FASES.values():
    d["F"] = grados_libertad(d["pi"], d["N"], d["r"], d["s"])

titulo("8. Regla de las fases F = 2 − π + N − r − s")
c, g = FASES["con Si(s)"], FASES["sin Si(s)"]
print(f"Con Si(s): π = {c['pi']}, N = {c['N']}, r = {c['r']}, s = {c['s']} "
      f"(la alimentación solo fija Cl:H en el gas) -> F = {c['F']}")
print(f"Sin Si(s): π = {g['pi']}, N = {g['N']}, r = {g['r']} (rango 3 con "
      f"{g['N']} especies), s = {g['s']} (fija Si:Cl y Cl:H) -> F = {g['F']}")
print("En ambos casos, con T y P fijas el estado de equilibrio queda determinado.")

# %% 9. ΔG° y K de cada reacción a 1300 K
i_si = ESPECIES.index("Si")
i_gas = [ESPECIES.index(s) for s in GAS]
i_comp = [ESPECIES.index(s) for s in COMPONENTES]
DGF_vec = np.array([DGF[s] for s in ESPECIES])


def termo(nu_mat):
    """(ν float, ΔG° kJ/mol, ln K, δ_gas) de un conjunto de reacciones."""
    nu_m = np.array(nu_mat, dtype=float)
    dG = nu_m @ DGF_vec
    return nu_m, dG, -dG * 1000.0 / (R_GAS * T), nu_m[:, i_gas].sum(axis=1)


nu_f, DG, lnK, delta_gas = termo(nu)
K = np.exp(lnK)
nu_ant_f, DG_ant, lnK_ant, delta_ant = termo(nu_ant)
K_ANT = np.exp(lnK_ant)

titulo(f"9. ΔG° y K a {T:g} K (ΔfG° de NIST-JANAF, R = {R_GAS} J/mol/K)")
print(f"{'':4}{'reacción':38}{'ΔG° (kJ/mol)':>14}{'K':>14}{'δ_gas':>7}")
for j, fila in enumerate(nu):
    print(f"{nombres_rx[j]:4}{REACCIONES[j]:38}"
          f"{DG[j]:14.3f}{K[j]:14.5g}{delta_gas[j]:7.0f}")
print(f"\nBase anterior ({', '.join(COMPONENTES_ANT)}), solo para verificar:")
for j, fila in enumerate(nu_ant):
    print(f"{'R' + str(j + 1) + 'ant':6}{reaccion_texto(fila, ESPECIES):36}"
          f"{DG_ant[j]:14.3f}{K_ANT[j]:14.5g}{delta_ant[j]:7.0f}")

# %% 10. Solver con residuales logarítmicos y cambio de régimen
A_gas_f = np.array([[float(x) for x in A[ESPECIES.index(s)]] for s in GAS])
NO_COMP = red["no_componentes"]          # especie formada por cada reacción
J_SI = NO_COMP.index("Si")               # R1: la única reacción con Si(s)


def ln_actividades(u) -> dict:
    """ln a de cada no componente (ln y si es gas, ln a_Si para Si(s))
    dadas u = ln y de las componentes."""
    return {sp: lnK[j] - delta_gas[j] * log(P / P0) - nu_f[j, i_comp] @ u
            for j, sp in enumerate(NO_COMP)}


def ln_y_gas(u) -> np.ndarray:
    """ln y de cada gas (orden GAS): componentes = u; resto por R2–R6."""
    ln = dict(zip(COMPONENTES, u))
    ln.update({sp: v for sp, v in ln_actividades(u).items() if sp != "Si"})
    return np.array([ln[s] for s in GAS])


def ln_a_si(u) -> float:
    return ln_actividades(u)["Si"]                        # por R1


def moles_alim(clorosilano: str, R: float) -> dict:
    return {clorosilano: 1.0, "H2": float(R)}


def elementos_alim(alim: dict) -> np.ndarray:
    """Moles de Si, Cl, H en la alimentación."""
    return sum(n * np.array([float(x) for x in A[ESPECIES.index(sp)]])
               for sp, n in alim.items())


def residuales(u, b, con_solido: bool):
    y = np.exp(ln_y_gas(u))
    e = y @ A_gas_f                                       # Σy·(Si, Cl, H)
    f = [log(y.sum()), log(e[1] / e[2]) - log(b[1] / b[2])]
    f.append(ln_a_si(u) if con_solido                     # R1 con a_Si = 1
             else log(e[0] / e[1]) - log(b[0] / b[1]))    # balance Si:Cl
    return f


def resolver(b, con_solido: bool):
    """Multi-arranque; devuelve (u, n_soluciones, dispersión entre ellas)."""
    sols = []
    for u1 in (-1.0, -3.0, -6.0, -9.0):                   # ln y_SiCl4
        for u2 in (-2.0, -4.0, -7.0):                     # ln y_SiCl3H
            for u3 in (-0.05, -0.5, -1.5):                # ln y_H2
                s = least_squares(residuales, [u1, u2, u3], args=(b, con_solido),
                                  bounds=(-60.0, 0.0), xtol=1e-15, ftol=1e-15,
                                  gtol=1e-15, max_nfev=500)
                if np.max(np.abs(s.fun)) < 1e-12:
                    sols.append(s.x)
    if not sols:
        raise RuntimeError("el solver no convergió en ningún arranque")
    sols = np.array(sols)
    return sols[0], len(sols), float(np.ptp(sols, axis=0).max())


def n_si_con_solido(b) -> float:
    """n_Si(s) suponiendo Si(s) presente (puede salir < 0)."""
    u, _, _ = resolver(b, True)
    e = np.exp(ln_y_gas(u)) @ A_gas_f
    return b[0] - b[1] / e[1] * e[0]


def equilibrio(clorosilano: str, R: float) -> dict:
    b = elementos_alim(moles_alim(clorosilano, R))
    u, n_sol, disp = resolver(b, True)
    y = np.exp(ln_y_gas(u))
    e = y @ A_gas_f
    n_T = b[1] / e[1]
    n_si_supuesto = b[0] - n_T * e[0]
    if n_si_supuesto >= 0:
        return dict(regimen="Si(s)", y=y, n_T=n_T, n_Si=n_si_supuesto,
                    a_Si=1.0, n_Si_supuesto=n_si_supuesto, b=b,
                    n_sol=n_sol, disp=disp)
    u, n_sol, disp = resolver(b, False)
    y = np.exp(ln_y_gas(u))
    n_T = b[1] / (y @ A_gas_f)[1]
    return dict(regimen="solo gas", y=y, n_T=n_T, n_Si=0.0, a_Si=exp(ln_a_si(u)),
                n_Si_supuesto=n_si_supuesto, b=b, n_sol=n_sol, disp=disp)


RESULTADOS = {(cl, R): equilibrio(cl, R) for cl in CLOROSILANOS
              for R in RELACIONES}

titulo(f"10. Equilibrio a {T:g} K y {P:g} bar, por mol de clorosilano")
print("Incógnitas: ln y de SiCl4, SiCl3H y H2. Con Si(s): R1–R6 con a_Si = 1.")
print("n_Si(s)* = valor si se supone Si(s) presente; si es < 0 se resuelve "
      "solo gas\n(R2–R6 + balance Si:Cl; n_Si(s) = 0 y a_Si sale de R1).")
for cl in CLOROSILANOS:
    print(f"\n{cl} + R H2")
    print(f"{'R':>5} {'régimen':>9} {'n_T':>9} {'n_Si(s)':>9} {'n_Si(s)*':>10} "
          f"{'a_Si':>8} {'arranques':>9}")
    for R in RELACIONES:
        e = RESULTADOS[(cl, R)]
        print(f"{R:5g} {e['regimen']:>9} {e['n_T']:9.5f} {e['n_Si']:9.5f} "
              f"{e['n_Si_supuesto']:10.5f} {e['a_Si']:8.5f} {e['n_sol']:9d}")
for cl in CLOROSILANOS:
    print(f"\nFracciones mol del gas, {cl} + R H2")
    print(f"{'R':>5} " + " ".join(f"{s:>10}" for s in GAS))
    for R in RELACIONES:
        y = RESULTADOS[(cl, R)]["y"]
        print(f"{R:5g} " + " ".join(f"{v:10.3e}" for v in y))
disp_max = max(e["disp"] for e in RESULTADOS.values())
print(f"\nTodos los arranques convergidos llegan a la misma solución "
      f"(dispersión máx. en ln y: {disp_max:.1e}).")

# %% 11. Verificación: minimización directa de G
MU_RT = np.array([DGF[s] * 1000.0 / (R_GAS * T) for s in GAS])
A_todas = np.vstack([[1.0, 0.0, 0.0], A_gas_f])           # Si(s) + gases


def G_RT(n):
    ng = np.maximum(n[1:], 1e-300)                        # n[0] = Si(s), μ° = 0
    return float(ng @ (MU_RT + np.log(ng / ng.sum()) + log(P / P0)))


def grad_G_RT(n):
    ng = np.maximum(n[1:], 1e-300)
    return np.append(0.0, MU_RT + np.log(ng / ng.sum()) + log(P / P0))


def minimizar_G(clorosilano: str, R: float):
    """min G/RT sujeto a balances de Si, Cl, H y n ≥ 0 (6 arranques)."""
    alim = moles_alim(clorosilano, R)
    b = elementos_alim(alim)
    rng = np.random.default_rng(0)
    mejor = None
    for intento in range(6):
        n0 = np.full(9, 0.05) if intento == 0 else 0.5 * rng.random(9)
        for sp, n in alim.items():
            n0[1 + GAS.index(sp)] += n
        s = minimize(G_RT, n0, jac=grad_G_RT, method="SLSQP",
                     bounds=[(0.0, None)] + [(1e-25, None)] * len(GAS),
                     constraints=[{"type": "eq", "fun": lambda n: n @ A_todas - b,
                                   "jac": lambda n: A_todas.T}],
                     options={"ftol": 1e-16, "maxiter": 5000})
        if mejor is None or s.fun < mejor.fun:
            mejor = s
    n = mejor.x
    return n[1:] / n[1:].sum(), n[0]


def max_res_K(nu_m, lnK_m, delta_m, ln_y, a_si):
    """máx |ln Q − ln K|; Q incluye a_Si^ν_Si (ν_Si = 0 si no hay Si)."""
    ln_Q = (nu_m[:, i_gas] @ ln_y + delta_m * log(P / P0)
            + nu_m[:, i_si] * log(a_si))
    return float(np.max(np.abs(ln_Q - lnK_m)))


VERIF = {}
for (cl, R), e in RESULTADOS.items():
    y_min, n_si_min = minimizar_G(cl, R)
    n_gas = e["y"] * e["n_T"]
    elem = n_gas @ A_gas_f + np.array([e["n_Si"], 0.0, 0.0])
    ln_y = np.log(e["y"])
    VERIF[(cl, R)] = dict(
        dy_minG=float(np.max(np.abs(y_min - e["y"]))),
        dnSi_minG=abs(n_si_min - e["n_Si"]),
        balance=float(np.max(np.abs(elem - e["b"]) / e["b"])),
        suma_y=abs(e["y"].sum() - 1.0),
        res_K=max_res_K(nu_f, lnK, delta_gas, ln_y, e["a_Si"]),
        res_K_ant=max_res_K(nu_ant_f, lnK_ant, delta_ant, ln_y, e["a_Si"]))

titulo("11. Verificación: minimización de G, balances, Σy = 1 y las K")
print("|lnQ−lnK| R1–R6: esta base; base ant.: R1–R6 con Si, HCl, H2 "
      "(mismas ΔfG°).")
print(f"{'':7}{'R':>5} {'|Δy| minG':>10} {'|ΔnSi| minG':>12} "
      f"{'balance rel':>12} {'|Σy−1|':>9} {'|lnQ−lnK|':>10} {'base ant.':>10}")
for (cl, R), v in VERIF.items():
    print(f"{cl:7}{R:5g} {v['dy_minG']:10.1e} {v['dnSi_minG']:12.1e} "
          f"{v['balance']:12.1e} {v['suma_y']:9.1e} {v['res_K']:10.1e} "
          f"{v['res_K_ant']:10.1e}")
peor = {c: max(v[c] for v in VERIF.values()) for c in next(iter(VERIF.values()))}
gas_ok = all(e["a_Si"] < 1 for e in RESULTADOS.values() if e["regimen"] == "solo gas")
print(f"\nPeor caso: |Δy| = {peor['dy_minG']:.1e}, |Δn_Si(s)| = "
      f"{peor['dnSi_minG']:.1e}, balance = {peor['balance']:.1e}, "
      f"|Σy−1| = {peor['suma_y']:.1e},\n|ln Q − ln K| = {peor['res_K']:.1e} "
      f"(esta base), {peor['res_K_ant']:.1e} (base anterior)")
print("En el régimen solo gas, a_Si sale de R1, así que R1 se cumple por "
      "construcción;\nla comprobación real son R2–R6 y la base anterior.")
print(f"Régimen solo gas con a_Si < 1 en todos los puntos: {'SÍ' if gas_ok else 'NO'}")

# %% 12. Mismo Cl:H -> mismo gas: SiCl4 con R = 10 vs SiCl3H con R = 7
e4, e3 = RESULTADOS[("SiCl4", 10)], RESULTADOS[("SiCl3H", 7)]
DY_EQUIV = float(np.max(np.abs(e4["y"] - e3["y"])))

titulo("12. SiCl4 + 10 H2 vs SiCl3H + 7 H2 (Cl:H = 4:20 = 3:15 = 1:5)")
print(f"{'':10}{'SiCl4, R=10':>14}{'SiCl3H, R=7':>14}")
for s, a, b_ in zip(GAS, e4["y"], e3["y"]):
    print(f"y_{s:8}{a:14.6e}{b_:14.6e}")
print(f"{'n_T':10}{e4['n_T']:14.6f}{e3['n_T']:14.6f}")
print(f"{'n_Si(s)':10}{e4['n_Si']:14.6f}{e3['n_Si']:14.6f}")
print(f"\nmáx |Δy| = {DY_EQUIV:.1e}: con Si(s) presente (F = 2, s = 1) el gas "
      "depende solo de T, P y Cl:H;\nn_T y n_Si(s) sí cambian.")

# %% 13. R mínimo para depositar Si desde SiCl4
R_MIN = brentq(lambda R: n_si_con_solido(elementos_alim(moles_alim("SiCl4", R))),
               5.0, 10.0, xtol=1e-12)

titulo("13. R mínimo para depositar Si(s) desde SiCl4")
print(f"n_Si(s) = 0 en R_min = {R_MIN:.4f}; por debajo, todo el Si queda en "
      "el gas (a_Si < 1).")
print("Desde SiCl3H hay depósito en todo el intervalo estudiado (R ≥ 1).")

# %% 14. Resumen: n_Si(s) por mol de clorosilano
R_RESUMEN = [1, 5, 10, 20, 50, 100, 150]
titulo("14. Resumen: n_Si(s) por mol de clorosilano")
print(f"{'R':8}" + "".join(f"{R:>9g}" for R in R_RESUMEN))
for cl in CLOROSILANOS:
    print(f"{cl + ':':8}" + "".join(f"{RESULTADOS[(cl, R)]['n_Si']:9.5f}"
                                    for R in R_RESUMEN))
print(f"R_min (SiCl4): {R_MIN:.4f}")

# %% 15. Salidas: CSV y gráfica n_Si(s) vs R
CSV_PATH = SALIDAS / "resultados_s24_ej3_si_epitaxial.csv"
PNG_PATH = SALIDAS / "grafica_s24_ej3_si_epitaxial.png"

if GUARDAR:
    SALIDAS.mkdir(exist_ok=True)
    with open(CSV_PATH, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["clorosilano", "R_H2", "regimen", "n_T", "n_Si(s)", "a_Si"]
                   + [f"y_{s}" for s in GAS])
        for (cl, R), e in RESULTADOS.items():
            w.writerow([cl, R, e["regimen"], f"{e['n_T']:.8g}",
                        f"{e['n_Si']:.8g}", f"{e['a_Si']:.8g}"]
                       + [f"{v:.8e}" for v in e["y"]])

    from matplotlib.figure import Figure
    from matplotlib.ticker import NullFormatter, ScalarFormatter

    TINTA, TINTA_2, REJILLA = "#0b0b0b", "#52514e", "#e4e3df"
    ESTILO = {"SiCl4": dict(color="#2a78d6", marker="o"),
              "SiCl3H": dict(color="#eb6834", marker="s")}
    fig = Figure(figsize=(7.5, 4.8), dpi=150, facecolor="#fcfcfb")
    ax = fig.add_subplot(facecolor="#fcfcfb")
    for cl in CLOROSILANOS:
        ys = [RESULTADOS[(cl, R)]["n_Si"] for R in RELACIONES]
        xs_l, ys_l = list(RELACIONES), list(ys)
        if cl == "SiCl4":                    # quiebre real en (R_min, 0)
            i = int(np.searchsorted(xs_l, R_MIN))
            xs_l.insert(i, R_MIN)
            ys_l.insert(i, 0.0)
        ax.plot(xs_l, ys_l, lw=2, color=ESTILO[cl]["color"])
        ax.plot(RELACIONES, ys, lw=0, ms=6, markeredgecolor="#fcfcfb",
                markeredgewidth=1, **ESTILO[cl])
        ax.plot([], [], lw=2, ms=6, label=f"{cl} + R H₂",
                markeredgecolor="#fcfcfb", markeredgewidth=1, **ESTILO[cl])
        ax.annotate(cl, (RELACIONES[-1], ys[-1]), xytext=(6, 0),
                    textcoords="offset points", va="center", color=TINTA,
                    fontsize=9)
    ax.axvline(R_MIN, color=TINTA_2, lw=1, ls="--")
    ax.annotate(f"R mín. (SiCl₄) = {R_MIN:.2f}", (R_MIN, 0.72), xytext=(5, 0),
                textcoords="offset points", color=TINTA_2, fontsize=9)
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100, 150])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(0.8, 260)
    ax.set_ylim(-0.02, 0.85)
    ax.set_xlabel("R = mol H₂ / mol clorosilano (escala log)", color=TINTA)
    ax.set_ylabel("n_Si(s) / mol de clorosilano", color=TINTA)
    ax.set_title(f"Si(s) depositado en el equilibrio, {T:g} K y {P:g} bar",
                 color=TINTA, fontsize=11, loc="left")
    ax.grid(True, which="major", color=REJILLA, lw=0.8)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color(TINTA_2)
    ax.tick_params(colors=TINTA_2)
    ax.legend(frameon=False, loc="upper left", labelcolor=TINTA)
    fig.tight_layout()
    fig.savefig(PNG_PATH, facecolor=fig.get_facecolor())

    titulo("15. Salidas")
    print(f"CSV:     {CSV_PATH.relative_to(SALIDAS.parent)}")
    print(f"Gráfica: {PNG_PATH.relative_to(SALIDAS.parent)}")

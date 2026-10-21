"""
tests.py — Casos de respuesta CONOCIDA para auditar el solver antes del parcial.

    python3 tests.py        (código de salida 0 = todo PASS)

Casos:
  a) Reformado de propano (3 rxn, K de 7e-3 a 2.4e10) a 700/900/1000 K
  b) Craqueo de n-butano a 750 K (SVA ej. 13.9): 2 reacciones paralelas
  c) Esterificación en fase líquida (Fogler): A+B<->C+D, Kc = 9, X_e = 0.9397
  d) Gas A <-> B + 2C con Kc = 0.7 mol²/dm⁶ (Fogler): X_e = 0.4425
  e) Water-gas shift con K = 1 (SVA ej. 13.5): xi = 1/2 y 2/3
  f) Termo NH3: K298 = exp(16450/RT0) ~ 764; Van't Hoff K500 ~ 0.4147
  g) K(T) por ec. 13.18 (datos SVA) vs las K del curso para el reformado
  h) Independencia de la semilla inicial del solver
  i) Parser de reacciones (fórmulas, coeficientes, balance)
  ... j–l) ver cada función
  m) Energía/combustión: HHV/LHV, T flama CH4 (SVA 4.7), Problemas 3–5 del 2do parcial, Dulong
  n) Presión/no idealidad: Problema 1 (NH3, Lewis-Randall con las dos
     convenciones de K_phi de core/convenciones.py, run_exercise con
     fugacidad="lewis") y 2 (Margules)
"""

from __future__ import annotations

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto

import sys
from math import exp

import numpy as np

from core.equilibrium import R_ATM, StoichiometricTable, solve_extents
from core.reactions import build_system, parse_formula, parse_reaction
from core.thermo import K_vant_hoff, R_J, T_REF, reaction_thermo

FALLAS = []


def check(label: str, got: float, exp_val: float,
          atol: float = None, rtol: float = None) -> None:
    if atol is not None:
        ok = abs(got - exp_val) <= atol
        crit = f"atol {atol:g}"
    else:
        ok = abs(got - exp_val) <= rtol * abs(exp_val)
        crit = f"rtol {rtol:g}"
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: obtenido {got:.6g}, "
          f"esperado {exp_val:.6g} ({crit})")
    if not ok:
        FALLAS.append(label)


def head(t: str) -> None:
    print(f"\n{'=' * 72}\n{t}\n{'=' * 72}")


# Datos SVA (Tablas C.1 y C.4): Hf298/Gf298 en J/mol; Cp/R = A+BT+CT²+DT⁻²
DATOS_SVA = {
    "CH4":  {"Hf298": -74520.0, "Gf298": -50460.0,
             "A": 1.702, "B": 9.081e-3, "C": -2.164e-6, "D": 0.0},
    "C2H4": {"Hf298": 52510.0, "Gf298": 68460.0,
             "A": 1.424, "B": 14.394e-3, "C": -4.392e-6, "D": 0.0},
    "C2H6": {"Hf298": -83820.0, "Gf298": -31855.0,
             "A": 1.131, "B": 19.225e-3, "C": -5.561e-6, "D": 0.0},
    "C3H6": {"Hf298": 19710.0, "Gf298": 62205.0,
             "A": 1.637, "B": 22.706e-3, "C": -6.915e-6, "D": 0.0},
    "C3H8": {"Hf298": -104680.0, "Gf298": -24290.0,
             "A": 1.213, "B": 28.785e-3, "C": -8.824e-6, "D": 0.0},
    "H2":   {"Hf298": 0.0, "Gf298": 0.0,
             "A": 3.249, "B": 0.422e-3, "C": 0.0, "D": 0.083e5},
    "N2":   {"Hf298": 0.0, "Gf298": 0.0,
             "A": 3.280, "B": 0.593e-3, "C": 0.0, "D": 0.040e5},
    "CO":   {"Hf298": -110525.0, "Gf298": -137169.0,
             "A": 3.376, "B": 0.557e-3, "C": 0.0, "D": -0.031e5},
    "CO2":  {"Hf298": -393509.0, "Gf298": -394359.0,
             "A": 5.457, "B": 1.045e-3, "C": 0.0, "D": -1.157e5},
    "H2O":  {"Hf298": -241818.0, "Gf298": -228572.0,
             "A": 3.470, "B": 1.450e-3, "C": 0.0, "D": 0.121e5},
    "NH3":  {"Hf298": -46110.0, "Gf298": -16450.0,
             "A": 3.578, "B": 3.020e-3, "C": 0.0, "D": -0.186e5},
}


# ---------------------------------------------------------------- caso a
def caso_a():
    head("a) Reformado de propano — 3 reacciones, 1 C3H8 + 4 H2O + 0.5 N2, "
         "1 bar")
    reacciones = ["C3H8 + 3 H2O = 3 CO + 7 H2",
                  "CO + H2O = CO2 + H2",
                  "C3H8 + 2 H2 = 3 CH4"]
    alim = {"C3H8": 1.0, "H2O": 4.0, "N2": 0.5}
    especies, nu, atomos, _ = build_system(reacciones, alim)
    clave = {700: (7.2285e-3, 7.5034, 2.4203e10, (0.24044, 0.69089, 0.75956)),
             900: (1.3407e6, 1.5591, 2.4864e8, (0.51886, 0.69625, 0.48114)),
             1000: (1.0501e9, 0.89956, 5.0080e7, (0.75718, 0.41774, 0.24282))}
    for T, (K1, K2, K3, xis) in clave.items():
        res = solve_extents(especies, alim, nu, [K1, K2, K3], P=1.0,
                            atoms=atomos)
        for j, xe in enumerate(xis):
            check(f"{T} K  xi_{j + 1}", res.xi[j], xe, atol=2e-4)
        if T == 1000:
            check("1000 K  y_H2", res.mole_fractions["H2"], 0.52099,
                  atol=1e-4)


# ---------------------------------------------------------------- caso b
def caso_b():
    head("b) Craqueo de n-butano a 750 K, 1 bar (SVA ej. 13.9)")
    reacciones = ["C4H10 = C2H4 + C2H6", "C4H10 = C3H6 + CH4"]
    alim = {"C4H10": 1.0}
    especies, nu, atomos, _ = build_system(reacciones, alim)
    res = solve_extents(especies, alim, nu, [3.856, 268.4], P=1.0,
                        atoms=atomos)
    check("xi_1 (C2H4+C2H6)", res.xi[0], 0.1068, atol=2e-4)
    check("xi_2 (C3H6+CH4)", res.xi[1], 0.8916, atol=2e-4)


# ---------------------------------------------------------------- caso c
def caso_c():
    head("c) Esterificación líquida A+B<->C+D, 23%/77% masa, Kc = 9 "
         "(delta = 0)")
    M_A, M_B = 46.068, 60.052
    nA0, nB0 = 23.0 / M_A, 77.0 / M_B            # base: 100 g de mezcla
    # Método A: líquido => K = Π x^ν, sin factor de presión (P = P0)
    res = solve_extents(["EtOH", "AcOH", "Ester", "H2O"],
                        {"EtOH": nA0, "AcOH": nB0},
                        [[-1, -1, 1, 1]], [9.0], P=1.0, P0=1.0)
    check("Método A: X_e = xi/n_A0", res.xi[0] / nA0, 0.9397, atol=2e-4)
    # Método B: volumen constante; con delta = 0 C_A0 se cancela en Kc
    tbl = StoichiometricTable(nu={"EtOH": -1, "AcOH": -1, "Ester": 1,
                                  "H2O": 1},
                              theta={"EtOH": 1.0, "AcOH": nB0 / nA0},
                              limiting="EtOH", T0=300.0, P0=1.0,
                              variable_volume=False)
    check("Método B: X_e de brentq", tbl.solve_Xe(9.0), 0.9397, atol=2e-4)


# ---------------------------------------------------------------- caso d
def caso_d():
    head("d) Gas A<->B+2C, C_A0 = F_A0/v0 = 2 mol/dm³, Kc = 0.7 mol²/dm⁶ "
         "(Fogler)")
    T0 = 300.0
    P0 = 2.0 * R_ATM * T0            # atm tales que C_A0 = P0/(R T0) = 2
    tbl = StoichiometricTable(nu={"A": -1.0, "B": 1.0, "C": 2.0},
                              theta={"A": 1.0}, limiting="A",
                              T0=T0, P0=P0, pressure_unit="atm",
                              variable_volume=True)
    check("Método B: X_e", tbl.solve_Xe(0.7), 0.4425, atol=5e-4)
    # Método A: K_y = Kc/C_T^delta = 0.7/2² = 0.175 con P = P0
    res = solve_extents(["A", "B", "C"], {"A": 1.0}, [[-1, 1, 2]],
                        [0.7 / 2.0 ** 2], P=1.0, P0=1.0)
    check("Método A: X = xi", res.xi[0], 0.4425, atol=5e-4)


# ---------------------------------------------------------------- caso e
def caso_e():
    head("e) Water-gas shift a 1100 K con K = 1 (SVA ej. 13.5)")
    reacciones = ["CO + H2O = CO2 + H2"]
    for alim, xe in ({"CO": 1.0, "H2O": 1.0}, 0.5), \
                    ({"CO": 1.0, "H2O": 2.0}, 2.0 / 3.0):
        especies, nu, atomos, _ = build_system(reacciones, alim)
        res = solve_extents(especies, alim, nu, [1.0], P=1.0, atoms=atomos)
        check(f"xi con {alim['H2O']:g} mol H2O", res.xi[0], xe, atol=1e-6)


# ---------------------------------------------------------------- caso f
def caso_f():
    head("f) Amoniaco 1/2 N2 + 3/2 H2 <-> NH3 — módulo termo")
    nu = {"N2": -0.5, "H2": -1.5, "NH3": 1.0}
    t298 = reaction_thermo(nu, DATOS_SVA, T_REF)
    check("K(298.15) = exp(16450/RT0)", t298["K"],
          exp(16450.0 / (R_J * 298.15)), rtol=1e-9)
    check("K(298) vs ~764 (redondeo del curso)", t298["K"], 764.0, rtol=5e-3)
    K500_vh = K_vant_hoff(500.0, t298["K"], -46110.0)
    check("Van't Hoff K(500) vs ~0.4147", K500_vh, 0.4147, rtol=1e-2)
    t500 = reaction_thermo(nu, DATOS_SVA, 500.0)
    print(f"  [info] con Cp(T) completo (ec. 13.18): K(500) = "
          f"{t500['K']:.4f}, dH(500) = {t500['dH'] / 1000:.2f} kJ/mol")


# ---------------------------------------------------------------- caso g
def caso_g():
    head("g) Módulo termo: ec. 13.18 vs integración numérica de Van't Hoff, "
         "y ancla física")
    from scipy.integrate import quad
    from core.thermo import icph

    rxns = [("R1 reformado", parse_reaction("C3H8 + 3 H2O = 3 CO + 7 H2")),
            ("WGS", parse_reaction("CO + H2O = CO2 + H2")),
            ("metanación", parse_reaction("C3H8 + 2 H2 = 3 CH4"))]

    # 1) la forma cerrada 13.18 debe coincidir con integrar
    #    d(lnK)/dT = dH(T)/(R T^2) numéricamente
    for nombre, nud in rxns:
        dA, dB, dC, dD = (sum(c * DATOS_SVA[s][k] for s, c in nud.items())
                          for k in ("A", "B", "C", "D"))
        dH0 = sum(c * DATOS_SVA[s]["Hf298"] for s, c in nud.items())
        dG0 = sum(c * DATOS_SVA[s]["Gf298"] for s, c in nud.items())
        for T in (700.0, 1000.0):
            lnK_num = (-dG0 / (R_J * T_REF)
                       + quad(lambda t: (dH0 + R_J * icph(T_REF, t, dA, dB,
                                                          dC, dD))
                              / (R_J * t * t), T_REF, T, limit=200)[0])
            K_forma = reaction_thermo(nud, DATOS_SVA, T)["K"]
            check(f"13.18 vs numérica: {nombre} a {T:g} K",
                  K_forma, exp(lnK_num), rtol=1e-10)

    # 2) ancla física: la MISMA WGS del caso (e): K = 1 a 1100 K
    K_wgs_1100 = reaction_thermo(rxns[1][1], DATOS_SVA, 1100.0)["K"]
    check("K_WGS(1100 K) desde C.1/C.4 vs 1 (SVA ej. 13.5)",
          K_wgs_1100, 1.0, rtol=0.10)

    # 3) informativo: las K del EJERCICIO de reformado del curso NO salen
    #    de las tablas C.1/C.4 de SVA (difieren hasta ~80x en K1); son
    #    datos DADOS de ese enunciado. La WGS de esa tabla daría K = 1 a
    #    ~970 K, contradiciendo el propio caso (e). En el examen: usa las
    #    K que te den (modo 'directo'); usa modo 'termo' cuando el
    #    enunciado pida CALCULAR K desde Tablas C.1/C.4.
    curso = {700.0: (7.2285e-3, 7.5034, 2.4203e10),
             900.0: (1.3407e6, 1.5591, 2.4864e8),
             1000.0: (1.0501e9, 0.89956, 5.0080e7)}
    print("  [info] K termo (SVA) vs K de la tabla del curso (dato dado):")
    for T, Ks in curso.items():
        vals = [reaction_thermo(nud, DATOS_SVA, T)["K"]
                for _, nud in rxns]
        print(f"         {T:g} K: " + "   ".join(
            f"K{j + 1} {v:.4g}/{k:.4g}" for j, (v, k) in
            enumerate(zip(vals, Ks))))


# ---------------------------------------------------------------- caso h
def caso_h():
    head("h) Independencia de la semilla inicial (reformado a 1000 K)")
    reacciones = ["C3H8 + 3 H2O = 3 CO + 7 H2", "CO + H2O = CO2 + H2",
                  "C3H8 + 2 H2 = 3 CH4"]
    alim = {"C3H8": 1.0, "H2O": 4.0, "N2": 0.5}
    especies, nu, atomos, _ = build_system(reacciones, alim)
    xis = [solve_extents(especies, alim, nu,
                         [1.0501e9, 0.89956, 5.0080e7], P=1.0,
                         atoms=atomos, seed=s).xi for s in range(5)]
    spread = max(float(np.max(np.abs(a - b))) for a in xis for b in xis)
    check("dispersión máx de xi entre semillas 0..4", spread, 0.0, atol=1e-8)


# ---------------------------------------------------------------- caso i
def caso_i():
    head("i) Parser de reacciones y fórmulas")
    f = parse_formula("Ca(OH)2")
    ok = f == {"Ca": 1, "O": 2, "H": 2}
    print(f"  [{'PASS' if ok else 'FAIL'}] parse_formula('Ca(OH)2') = {f}")
    if not ok:
        FALLAS.append("parse_formula")
    nu = parse_reaction("C3H8 + 3 H2O = 3 CO + 7 H2")
    ok = nu == {"C3H8": -1.0, "H2O": -3.0, "CO": 3.0, "H2": 7.0}
    print(f"  [{'PASS' if ok else 'FAIL'}] parse_reaction reformado = {nu}")
    if not ok:
        FALLAS.append("parse_reaction")
    nu2 = parse_reaction("SO2 + 1/2 O2 = SO3")
    ok = abs(nu2["O2"] + 0.5) < 1e-12
    print(f"  [{'PASS' if ok else 'FAIL'}] coeficiente '1/2 O2' -> "
          f"{nu2['O2']}")
    if not ok:
        FALLAS.append("coef 1/2")
    try:
        build_system(["C3H8 + H2O = CO + H2"], {"C3H8": 1})   # sin balancear
        print("  [FAIL] reacción desbalanceada NO fue rechazada")
        FALLAS.append("balance parser")
    except ValueError as e:
        print(f"  [PASS] desbalanceada rechazada: {str(e)[:70]}...")


# ---------------------------------------------------------------- caso j
def caso_j():
    head("j) data/sva_tables.json vs los valores del curso que me dictaste "
         "(bloque P3)")
    from core.thermo import cargar_datos
    # Referencia INDEPENDIENTE del JSON: los números que el usuario dictó
    # en el enunciado de P3 (examen_p2026_3).
    ref = {
        "C3H8": {"Hf298": -104680.0, "Gf298": -24290.0,
                 "A": 1.213, "B": 28.785e-3, "C": -8.824e-6, "D": 0.0},
        "CH4":  {"Hf298": -74520.0, "Gf298": -50460.0,
                 "A": 1.702, "B": 9.081e-3, "C": -2.164e-6, "D": 0.0},
        "C3H6": {"Hf298": 19710.0, "Gf298": 62205.0,
                 "A": 1.637, "B": 22.706e-3, "C": -6.915e-6, "D": 0.0},
    }
    datos = cargar_datos(list(ref))
    for sp, campos in ref.items():
        for campo, v in campos.items():
            check(f"{sp}.{campo}", datos[sp][campo], v,
                  atol=abs(v) * 1e-12 + 1e-15)
    try:
        cargar_datos(["C3H8", "XYZ99"])
        print("  [FAIL] especie inexistente NO fue rechazada")
        FALLAS.append("cargar_datos faltante")
    except KeyError as e:
        ok = "XYZ99" in str(e)
        print(f"  [{'PASS' if ok else 'FAIL'}] especie faltante rechazada "
              f"nombrándola: {str(e)[:75]}...")
        if not ok:
            FALLAS.append("cargar_datos mensaje")


# ---------------------------------------------------------------- caso k
def caso_k():
    head("k) Craqueo de n-butano a 500 K (Koretsky 14-5): K termo del JSON "
         "y xi a 1 y 25 bar")
    from core.thermo import cargar_datos
    nud = parse_reaction("C4H10 = C3H6 + CH4")
    datos = cargar_datos(["C4H10", "C3H6", "CH4"])
    t = reaction_thermo(nud, datos, 500.0)
    print(f"  [info] dH298 = {t['dH298']:.0f} J/mol, dG298 = {t['dG298']:.0f}"
          f" J/mol, dH(500) = {t['dH']:.0f} J/mol, dG(500) = {t['dG']:.0f}"
          f" J/mol")
    check("K(500 K) desde el JSON", t["K"], 1.13, rtol=0.02)
    especies, nu, atomos, _ = build_system(["C4H10 = C3H6 + CH4"],
                                           {"C4H10": 10.0})
    for P, xi_exp in ((1.0, 7.28), (25.0, 2.08)):
        res = solve_extents(especies, {"C4H10": 10.0}, nu, [t["K"]],
                            P=P, P0=1.0, atoms=atomos)
        check(f"xi a {P:g} bar", res.xi[0], xi_exp, atol=0.05)


# ---------------------------------------------------------------- caso l
def caso_l():
    head("l) nonideal.py: phi de Pitzer y gammas de Wilson (Sesiones 13/15)")
    import warnings as _w
    from core.nonideal import (phi_pitzer, wilson_gamma, wilson_gamma_inf,
                          wilson_lambdas)
    # phi -> 1 en el límite de gas ideal (Tr alto, Pr bajo)
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        r = phi_pitzer(1400.0, 1.0, 282.3, 50.40, 0.087)   # C2H4, Tr ~ 5
        check("phi(Tr~5, Pr~0.02) ~ 1", r["phi"], 1.0, atol=0.01)
        # valor cotejado a mano (H2O a 473.15 K y 34.5 bar, slides S13):
        # Tr=0.7312, Pr=0.1564, B0=-0.6134, B1=-0.5016 -> phi = 0.845
        r2 = phi_pitzer(473.15, 34.5, 647.1, 220.55, 0.345)
        check("phi H2O(473.15 K, 34.5 bar) vs mano", r2["phi"], 0.845,
              rtol=2e-3)
        check("  Tr reportada", r2["Tr"], 0.7312, rtol=1e-3)
        check("  Pr reportada", r2["Pr"], 0.1564, rtol=1e-3)
        # el warning de zona virial debe dispararse (EtOH a 34.5 bar)
        with _w.catch_warnings(record=True) as avisos:
            _w.simplefilter("always")
            r3 = phi_pitzer(473.15, 34.5, 513.9, 61.48, 0.645)
        ok = (not r3["valido"]) and len(avisos) == 1
        print(f"  [{'PASS' if ok else 'FAIL'}] warning fuera de zona "
              f"virial (EtOH: Vr = {r3['Vr']:.2f} < 2, "
              f"{len(avisos)} aviso)")
        if not ok:
            FALLAS.append("warning virial")
    # Wilson etanol(1)/agua(2), DECHEMA, 473.15 K
    W = dict(V1=58.69, V2=18.07, a12=382.30, a21=955.45)
    L12, L21 = wilson_lambdas(473.15, **W)
    check("Lambda12", L12, 0.205, rtol=2e-3)
    check("Lambda21", L21, 1.176, rtol=2e-3)
    g1, _ = wilson_gamma(1.0, 473.15, **W)
    _, g2 = wilson_gamma(0.0, 473.15, **W)
    check("gamma1(x1=1) = 1", g1, 1.0, atol=1e-12)
    check("gamma2(x2=1) = 1", g2, 1.0, atol=1e-12)
    g1i_num, _ = wilson_gamma(1e-12, 473.15, **W)
    _, g2i_num = wilson_gamma(1.0 - 1e-12, 473.15, **W)
    g1i, g2i = wilson_gamma_inf(473.15, **W)
    check("gamma1_inf numérica vs cerrada", g1i_num, g1i, rtol=1e-6)
    check("gamma2_inf numérica vs cerrada", g2i_num, g2i, rtol=1e-6)


# ---------------------------------------------------------------- caso m
def caso_m():
    head("m) Energía/combustión (core.energia): SVA ej. 4.7 y Problemas 3–5")
    import warnings as _w
    from core.energia import (T_salida_para_Q, balance_energia,
                              calor_combustion, combustion_mezcla,
                              alimentacion_por_atomos, avances_desde_flujos,
                              h_liquido_sat)
    # ΔH°c por ley de Hess (H2O líquida) contra valores de SVA/slides
    check("HHV CH4 (kJ/mol)", calor_combustion("CH4")["HHV"] / 1000,
          890.649, rtol=1e-6)
    check("LHV CH4 (kJ/mol)", calor_combustion("CH4")["LHV"] / 1000,
          802.625, rtol=1e-6)
    check("HHV H2 (slide: 285.8)", calor_combustion("H2")["HHV"] / 1000,
          285.83, rtol=1e-4)
    check("ΔH°c CO (slide: −283.0)", calor_combustion("CO")["dHc_l"] / 1000,
          -282.984, rtol=1e-5)
    # Problema 3 a) y b)
    check("P3 a) ΔH°c (kJ/mol)", combustion_mezcla(
        {"CH4": .95, "C2H6": .02, "C3H8": .02, "N2": .01})["dHc_l"] / 1000,
        -921.71, atol=0.01)
    check("P3 b) ΔH°c (kJ/mol)", combustion_mezcla(
        {"CH4": .90, "C2H6": .05, "C3H8": .03, "N2": .02})["dHc_l"] / 1000,
        -946.16, atol=0.05)
    # SVA ej. 4.7 / Sesión 18 ej. 1: T flama CH4, 20 % exceso → 2066 K
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        N2 = 2.4 * 79 / 21
        T = T_salida_para_Q([({"CH4": 1, "O2": 2.4, "N2": N2}, 298.15)],
                            {"CO2": 1, "H2O": 2, "O2": 0.4, "N2": N2},
                            ["CH4 + 2 O2 = CO2 + 2 H2O"], [1.0])
    check("T flama adiabática CH4 20 % exceso (K)", T, 2066.0, atol=1.0)
    # Problema 4: Q = 54 881 J por mol de producto
    sal = {"CO2": .0275, "CO": .1725, "H2O": .1725, "H2": .6275}
    rx = ["CH4 + H2O = CO + 3 H2", "CO + H2O = CO2 + H2"]
    alim = alimentacion_por_atomos(sal, ["CH4", "H2O"])
    xi, _ = avances_desde_flujos(alim, sal, rx)
    b = balance_energia([(alim, 773.15)], sal, 1123.15, rx, xi)
    check("P4 Q (J)", b["Q"], 54881.0, rtol=1e-4)
    check("P4 Q trayectoria vs formación", b["Q"] - b["Q_formacion"], 0.0,
          atol=1e-6)
    # Problema 5: V = 0.407 m³/s
    xi5 = [0.7 * 0.94, 0.3 * 0.94]
    O2 = 1.35 * 2 * 0.94
    aire = {"O2": O2, "N2": O2 * 79 / 21}
    sal5 = {"CO2": xi5[0], "CO": xi5[1], "H2O": 1.88,
            "O2": O2 - 2 * xi5[0] - 1.5 * xi5[1], "N2": 0.06 + aire["N2"]}
    b5 = balance_energia([({"CH4": .94, "N2": .06}, 298.15), (aire, 298.15)],
                         sal5, 483.15, ["CH4 + 2 O2 = CO2 + 2 H2O",
                                        "CH4 + 3/2 O2 = CO + 2 H2O"], xi5)
    dh = h_liquido_sat(95.0) - h_liquido_sat(25.0)
    V = 75 * 0.45359237 * dh * 1000 / (-b5["Q"]) * R_J * 298.15 / 101325
    check("P5 V (m³/s)", V, 0.407, atol=6e-4)
    # Ecuación de Dulong — columna "Coal" de la tabla de la Sesión 17 (daf):
    # C 0.8294, H 0.0262, O 0.0933, S 0.01. A mano:
    #   HHV = 33.9(0.8294) + 144.2(0.0262 − 0.0933/8) + 9.4(0.01)
    #       = 28.11666 + 2.0963075 + 0.094 = 30.3069675 MJ/kg
    #   LHV = 30.3069675 − 2.442·9·0.0262 = 29.7311439 MJ/kg
    from core.energia import hhv_dulong
    dl = hhv_dulong(0.8294, 0.0262, 0.0933, 0.01)
    check("Dulong Coal HHV (MJ/kg)", dl["HHV"], 30.3069675, rtol=1e-9)
    check("Dulong Coal LHV (MJ/kg)", dl["LHV"], 29.7311439, rtol=1e-9)
    check("Dulong Coal HHV (kJ/kg)", hhv_dulong(0.8294, 0.0262, 0.0933, 0.01,
                                                base="kJ/kg")["HHV"],
          30306.9675, rtol=1e-9)
    try:
        hhv_dulong(82.94, 2.62, 9.33, 1.0)
        print("  [FAIL] Dulong con porcentajes NO fue rechazado")
        FALLAS.append("dulong %")
    except ValueError as e:
        print(f"  [PASS] Dulong con % rechazado: {str(e)[:55]}...")


# ---------------------------------------------------------------- caso n
def caso_n():
    head("n) Presión/no idealidad 2do parcial: Problemas 1 y 2")
    from math import log
    from scipy.optimize import brentq
    from core.nonideal import K_phi_lewis, cargar_criticas, margules_gamma
    nu = {"NH3": 1.0, "N2": -0.5, "H2": -1.5}
    d = {k: DATOS_SVA[k] for k in nu}
    lnK = lambda T: log(reaction_thermo(nu, d, T)["K"])
    Ky = lambda x: (x / (2 - x)) / ((0.5 * (1 - x) / (2 - x)) ** 0.5
                                    * (1.5 * (1 - x) / (2 - x)) ** 1.5)
    xi = brentq(lambda x: log(Ky(x)) - lnK(300.0), 1e-12, 1 - 1e-12)
    check("P1 a) y_NH3 (300 K, 1 bar)", xi / (2 - xi), 0.9349, atol=1e-4)
    ky = Ky(2 / 3)
    check("P1 b) T (K)", brentq(lambda T: lnK(T) - log(ky), 300, 1500),
          399.5, atol=0.1)
    check("P1 c) T (K)", brentq(lambda T: lnK(T) - log(ky / 100), 300, 1500),
          577.6, atol=0.1)
    # --- core/convenciones.py: ÚNICO lugar donde se decide cómo entra K_phi
    from core.convenciones import CONVENCION_KPHI, K_efectiva, factor_Kphi
    print(f"  [info] convención default del proyecto: '{CONVENCION_KPHI}'")
    check("factor_Kphi(0.8, 'rigurosa') = K_phi", factor_Kphi(0.8, "rigurosa"),
          0.8, rtol=1e-12)
    check("factor_Kphi(0.8, 'maestra') = 1/K_phi", factor_Kphi(0.8, "maestra"),
          1.25, rtol=1e-12)
    check("K_efectiva(2, 0.8, 'rigurosa') = K/K_phi",
          K_efectiva(2.0, 0.8, "rigurosa"), 2.5, rtol=1e-12)
    check("K_efectiva(2, 0.8, 'maestra') = K·K_phi",
          K_efectiva(2.0, 0.8, "maestra"), 1.6, rtol=1e-12)
    check("factor_Kphi sin convención usa el default",
          factor_Kphi(0.8), factor_Kphi(0.8, CONVENCION_KPHI), rtol=1e-15)
    try:
        factor_Kphi(0.8, "otra")
        print("  [FAIL] convención inválida NO fue rechazada")
        FALLAS.append("convencion invalida")
    except ValueError as e:
        print(f"  [PASS] convención inválida rechazada: {str(e)[:60]}...")
    # --- P1 d) con las dos convenciones (K_phi Lewis-Randall, Pitzer)
    crit = cargar_criticas(list(nu))

    def T_d(conv):
        return brentq(lambda T: lnK(T) - log(
            factor_Kphi(K_phi_lewis(nu, T, 100.0, crit, avisar=False)[0], conv)
            * ky / 100), 300, 1500)

    check("P1 d) 'maestra' (K_y = K·P·K_phi) T (K)", T_d("maestra"), 570.4,
          atol=0.2)
    check("P1 d) 'rigurosa' (K_y = K·P/K_phi) T (K)", T_d("rigurosa"), 584.5,
          atol=0.2)
    # --- run_exercise con fugacidad="lewis" (K_ef = K/F) y sin fugacidad
    import contextlib, io, os
    from core.exercise import run_exercise

    def y_NH3_run(T, fug):
        with contextlib.redirect_stdout(io.StringIO()):
            r = run_exercise(titulo="test n NH3", reacciones=["1/2 N2 + 3/2 H2 = NH3"],
                             alimentacion={"N2": 0.5, "H2": 1.5}, fase="gas",
                             P=100.0, P0=1.0, unidad="bar", temperaturas=[T],
                             modo_K="termo", especies_termo=["N2", "H2", "NH3"],
                             fugacidad=fug, csv_out=os.devnull)
        return r["resultados"][T]["res"].mole_fractions["NH3"], r["veredicto"]

    y, ver = y_NH3_run(570.4, "lewis")
    check(f"run_exercise fugacidad='lewis' ({CONVENCION_KPHI}) 570.4 K, 100 bar: y_NH3",
          y, 0.500, atol=1e-3)
    ok = ver != "FAIL"
    print(f"  [{'PASS' if ok else 'FAIL'}] selfcheck de run_exercise con "
          f"fugacidad valida contra K_ef: {ver}")
    if not ok:
        FALLAS.append("selfcheck lewis")
    y, ver = y_NH3_run(577.64, None)
    check("run_exercise fugacidad=None 577.64 K, 100 bar: y_NH3 (caso ideal intacto)",
          y, 0.500, atol=1e-3)
    # Problema 2: Margules 1 parámetro
    K = exp(1000.0 / (R_J * 298.15))
    def F(x):
        gA, gB = margules_gamma(x, 0.1)
        return log((1 - x) * gB / (x * gA)) - log(K)
    check("P2 a) x_A (Margules)", brentq(F, 1e-9, 1 - 1e-9), 0.3955,
          atol=6e-5)
    check("P2 b) x_A (ideal)", 1 / (1 + K), 0.4005, atol=6e-5)
    g1, g2 = margules_gamma(0.0, 0.5, 1.2)
    check("Margules 2 par.: ln γ1∞ = A12", log(g1), 0.5, rtol=1e-12)
    g1, g2 = margules_gamma(1.0, 0.5, 1.2)
    check("Margules 2 par.: ln γ2∞ = A21", log(g2), 1.2, rtol=1e-12)


if __name__ == "__main__":
    for f in (caso_a, caso_b, caso_c, caso_d, caso_e, caso_f, caso_g,
              caso_h, caso_i, caso_j, caso_k, caso_l, caso_m, caso_n):
        f()
    print(f"\n{'=' * 72}")
    if FALLAS:
        print(f"RESULTADO: {len(FALLAS)} FALLAS -> " + "; ".join(FALLAS))
        sys.exit(1)
    print("RESULTADO: TODOS los tests PASS")
    sys.exit(0)

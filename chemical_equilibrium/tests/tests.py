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
  o) Reacciones independientes por el método matricial (S23–24,
     core/independencia.py): reformado de CH4, propileno, metanol,
     perclorato de amonio, dependencia y equivalencia de conjuntos
  p) Método matricial con sólido puro: oxidación de H2S con S8(s), π = 2
  q) Si epitaxial desde SiCl4/SiCl3H a 1300 K (ejercicios/s24_ej3): matrices,
     K, n_Si(s) con cambio de régimen, R mínimo y verificación contra min G
  r) Criterio de componentes (CLAUDE.md, elegir_componentes): Si epitaxial
     con SiCl4, SiCl3H, H2; reformado; H2S; perclorato; dependencias y errores
  s) ν leída de [−I | −Z] (CLAUDE.md): matrices ν por columnas, Aᵀν = 0,
     escalado a enteros, signo -> lado de la reacción, mismas reacciones
  t) Perclorato de amonio (ejercicios/s24_ej4): componentes del enunciado,
     [I | Z] con fracciones, ν por columnas y reacciones
  u) Minimización de Gibbs (S25, core/gibbs.py): SVA Ej. 13.14, contraste con
     Método A, steam cracker, trazas, presión, φ, sin arranque y errores
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


def _txt(v) -> str:   # Fraction(-3, 1) -> -3, también dentro de listas
    return "[" + ", ".join(map(_txt, v)) + "]" if isinstance(v, list) else str(v)


def igual(label: str, got, exp_val) -> None:
    """Comparación EXACTA (enteros, Fraction, listas, texto)."""
    ok = got == exp_val
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: obtenido {_txt(got)}, "
          f"esperado {_txt(exp_val)}")
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


# ---------------------------------------------------------------- caso o
def caso_o():
    head("o) Reacciones independientes, método matricial (S23–24)")
    from core.independencia import (
        comprobar_conservacion, generan_mismo_espacio, grados_libertad,
        matriz_atomica, matriz_estequiometrica, num_reacciones_independientes,
        rango, rango_reacciones, reaccion_texto, reducir_por_componentes,
        transpuesta)

    # --- Reformado de metano (resultados del ejercicio s24_ej1)
    sp = ["CH4", "H2O", "CO", "CO2", "H2"]
    A, el = matriz_atomica(sp, ["C", "H", "O"])
    At = transpuesta(A)
    r = num_reacciones_independientes(A)
    igual("reformado: rango", rango(At), 3)
    igual("reformado: r", r, 2)
    igual("reformado: F = 2 − π + N − r − s", grados_libertad(1, len(sp), r, 2), 2)
    red = reducir_por_componentes(At, sp, ["CH4", "H2O", "H2"])
    igual("reformado: columnas [I | Z]", red["columnas"],
          ["CH4", "H2O", "H2", "CO", "CO2"])
    igual("reformado: matriz reducida", red["reducida"],
          [[1, 0, 0, 1, 1], [0, 1, 0, 1, 2], [0, 0, 1, -3, -4]])
    nu = matriz_estequiometrica(red, sp)
    igual("reformado: ν R1", nu[0], [-1, -1, 1, 0, 3])
    igual("reformado: ν R2", nu[1], [-1, -2, 0, 1, 4])
    igual("reformado: texto R1", reaccion_texto(nu[0], sp),
          "CH4 + H2O -> CO + 3 H2")
    igual("reformado: texto R2", reaccion_texto(nu[1], sp),
          "CH4 + 2 H2O -> CO2 + 4 H2")
    igual("reformado: Aᵀν = 0", comprobar_conservacion(At, nu)[1], True)

    # --- r = N − rango para otros sistemas
    for nombre, especies, r_esp in (
            ("propileno", ["C3H6", "O3", "H2O", "CH2O", "C2H4O", "H2O2"], 3),
            ("metanol", ["CH3OH", "C", "CO", "CO2", "CH4", "H2O"], 3)):
        A2, _ = matriz_atomica(especies)
        igual(f"{nombre}: r", num_reacciones_independientes(A2), r_esp)
    sp3 = ["NH4ClO4", "Cl2", "N2O", "NOCl", "HCl", "H2O", "N2", "O2", "ClO2"]
    A3, _ = matriz_atomica(sp3)
    At3 = transpuesta(A3)
    igual("perclorato de amonio: r", num_reacciones_independientes(A3), 5)
    red3 = reducir_por_componentes(At3, sp3, ["NH4ClO4", "Cl2", "N2O", "NOCl"])
    nu3 = matriz_estequiometrica(red3, sp3)
    igual("perclorato de amonio: reacciones generadas", len(nu3), 5)
    igual("perclorato de amonio: rango de ν", rango(nu3), 5)
    igual("perclorato de amonio: Aᵀν = 0", comprobar_conservacion(At3, nu3)[1],
          True)
    try:
        spm = ["CH3OH", "C", "CO", "CO2", "CH4", "H2O"]
        Am, _ = matriz_atomica(spm)
        reducir_por_componentes(transpuesta(Am), spm, ["C", "CO", "CO2"])
        print("  [FAIL] componentes dependientes (C, CO, CO2) NO rechazados")
        FALLAS.append("componentes dependientes")
    except ValueError as e:
        print(f"  [PASS] componentes dependientes rechazados: {str(e)[:60]}...")

    # --- dependencia y equivalencia de conjuntos de reacciones
    R1 = "CH4 + H2O -> CO + 3 H2"
    R2 = "CH4 + 2 H2O -> CO2 + 4 H2"
    igual("dependencia: rango{R1, R2, 2 CH4 + 3 H2O -> CO + CO2 + 7 H2}",
          rango_reacciones([R1, R2, "2 CH4 + 3 H2O -> CO + CO2 + 7 H2"]), 2)
    igual("equivalencia: {R1, WGS} ~ {R1, R2}",
          generan_mismo_espacio([R1, "CO + H2O -> CO2 + H2"], [R1, R2]), True)
    igual("equivalencia: {CO2 + H2 -> CO + H2O, R2} ~ {R1, R2}",
          generan_mismo_espacio(["CO2 + H2 -> CO + H2O", R2], [R1, R2]), True)


# ---------------------------------------------------------------- caso p
def caso_p():
    head("p) Método matricial con sólido puro: H2S + O2 con S8(s)")
    from core.independencia import (
        comprobar_conservacion, grados_libertad, matriz_atomica,
        matriz_estequiometrica, matriz_reacciones, num_reacciones_independientes,
        rango, reaccion_texto, reducir_por_componentes, transpuesta)
    sp = ["H2S", "O2", "SO2", "SO3", "H2O", "S8"]
    A, _ = matriz_atomica(sp, ["H", "O", "S"])
    At = transpuesta(A)
    r = num_reacciones_independientes(A)
    igual("H2S/S8: rango", rango(At), 3)
    igual("H2S/S8: r", r, 3)
    igual("H2S/S8: F = 2 − π + N − r − s (π = 2, s = 1)",
          grados_libertad(2, len(sp), r, 1), 2)
    red = reducir_por_componentes(At, sp, ["H2S", "O2", "S8"])
    nu = matriz_estequiometrica(red, sp)
    igual("H2S/S8: reacciones generadas", len(nu), 3)
    igual("H2S/S8: Aᵀν = 0", comprobar_conservacion(At, nu)[1], True)
    # cada reacción esperada (texto del enunciado) debe estar en ν; el módulo
    # numera por no componente en el orden de ESPECIES (SO2, SO3, H2O)
    esperadas = {"R1": "8 H2S + 4 O2 -> 8 H2O + S8",
                 "R2": "S8 + 8 O2 -> 8 SO2",
                 "R3": "S8 + 12 O2 -> 8 SO3"}
    filas, _ = matriz_reacciones(list(esperadas.values()), sp)
    for (nombre, txt), fila in zip(esperadas.items(), filas):
        pos = [i for i, rx in enumerate(nu) if rx == fila]
        ok = len(pos) == 1
        donde = (f"generada como R{pos[0] + 1}: {reaccion_texto(nu[pos[0]], sp)}"
                 if ok else "NO está en ν")
        print(f"  [{'PASS' if ok else 'FAIL'}] H2S/S8: {nombre} esperada "
              f"'{txt}' -> {donde}")
        if not ok:
            FALLAS.append(f"H2S/S8 {nombre}")


# ---------------------------------------------------------------- caso q
_EJERCICIOS = {}


def cargar_ejercicio(nombre: str):
    """Importa ejercicios/<nombre>.py una sola vez, sin salida (importados,
    los ejercicios no escriben en salidas/)."""
    if nombre not in _EJERCICIOS:
        import contextlib, importlib.util, io
        ruta = (pathlib.Path(__file__).resolve().parents[1] / "ejercicios"
                / f"{nombre}.py")
        spec = importlib.util.spec_from_file_location(nombre, ruta)
        mod = importlib.util.module_from_spec(spec)
        with contextlib.redirect_stdout(io.StringIO()):
            spec.loader.exec_module(mod)
        _EJERCICIOS[nombre] = mod
    return _EJERCICIOS[nombre]


def cargar_ej3():
    return cargar_ejercicio("s24_ej3_si_epitaxial")


def caso_q():
    head("q) Si epitaxial a 1300 K: matrices, K, n_Si(s) y R mínimo (s24_ej3)")
    from core.independencia import (matriz_reacciones, reaccion_texto,
                                    reducir_por_componentes)
    ej = cargar_ej3()
    # --- matrices con la base Si, HCl, H2 (el ejercicio ahora usa las de
    #     alimentación; esta base se conserva como verificación)
    red_ant = reducir_por_componentes(ej.At, ej.ESPECIES, ["Si", "HCl", "H2"])
    igual("Si: rango", ej.rg, 3)
    igual("Si: r", ej.r, 6)
    igual("Si: columnas [I | Z] (base Si, HCl, H2)", red_ant["columnas"],
          ["Si", "HCl", "H2", "SiCl2", "SiCl4", "SiCl3H", "SiCl2H2",
           "SiClH3", "SiH4"])
    igual("Si: matriz reducida (base Si, HCl, H2)", red_ant["reducida"],
          [[1, 0, 0, 1, 1, 1, 1, 1, 1],
           [0, 1, 0, 2, 4, 3, 2, 1, 0],
           [0, 0, 1, -1, -2, -1, 0, 1, 2]])
    esperadas = ["Si + 2 HCl -> SiCl2 + H2", "Si + 4 HCl -> SiCl4 + 2 H2",
                 "Si + 3 HCl -> SiCl3H + H2", "Si + 2 HCl -> SiCl2H2",
                 "Si + HCl + H2 -> SiClH3", "Si + 2 H2 -> SiH4"]
    filas, _ = matriz_reacciones(esperadas, ej.ESPECIES)
    for j, (txt, fila) in enumerate(zip(esperadas, filas)):
        ok = ej.nu_ant[j] == fila
        print(f"  [{'PASS' if ok else 'FAIL'}] Si: R{j + 1} esperada '{txt}' -> "
              f"ν generada: {reaccion_texto(ej.nu_ant[j], ej.ESPECIES)}")
        if not ok:
            FALLAS.append(f"Si R{j + 1}")
    igual("Si: Aᵀν = 0", ej.conserva, True)
    for caso, (pi, N, r, s, F) in (("con Si(s)", (2, 9, 6, 1, 2)),
                                   ("sin Si(s)", (1, 8, 5, 2, 2))):
        d = ej.FASES[caso]
        igual(f"Si: {caso} (π, N, r, s, F)",
              [d["pi"], d["N"], d["r"], d["s"], d["F"]], [pi, N, r, s, F])
    # --- K a 1300 K (dadas a 5 cifras significativas)
    for j, K_esp in enumerate([2.6970, 1.9639e3, 8.9729e1, 5.7828e-1,
                               1.0476e-3, 7.8854e-7]):
        igual(f"Si: K{j + 1} a 5 cifras (base Si, HCl, H2)",
              float(f"{ej.K_ANT[j]:.4e}"), K_esp)
    # --- n_Si(s) por mol de clorosilano (dados a 5 decimales)
    tabla = {"SiCl4":  [0, 0, 0.01537, 0.09780, 0.29284, 0.50228, 0.62914],
             "SiCl3H": [0.21309, 0.24412, 0.28708, 0.36416, 0.53374, 0.69611,
                        0.78325]}
    for cl, valores in tabla.items():
        for R, n_esp in zip([1, 5, 10, 20, 50, 100, 150], valores):
            check(f"Si: n_Si(s) {cl}, R = {R}", ej.RESULTADOS[(cl, R)]["n_Si"],
                  n_esp, atol=5e-6)
    gas = [e for e in ej.RESULTADOS.values() if e["regimen"] == "solo gas"]
    ok = bool(gas) and all(e["n_Si"] == 0 and e["a_Si"] < 1 for e in gas)
    print(f"  [{'PASS' if ok else 'FAIL'}] Si: régimen solo gas en "
          f"{len(gas)} puntos, todos con n_Si(s) = 0 y a_Si < 1")
    if not ok:
        FALLAS.append("Si régimen solo gas")
    check("Si: R mínimo para depositar desde SiCl4", ej.R_MIN, 8.2327, atol=5e-5)
    # --- verificación
    peor = {c: max(v[c] for v in ej.VERIF.values())
            for c in ("dy_minG", "dnSi_minG", "balance", "suma_y", "res_K",
                      "res_K_ant")}
    check("Si: máx |Δy| vs minimización de G", peor["dy_minG"], 0, atol=1e-7)
    check("Si: máx |Δn_Si(s)| vs minimización de G", peor["dnSi_minG"], 0,
          atol=1e-7)
    check("Si: máx error relativo de balances Si, Cl, H", peor["balance"], 0,
          atol=1e-12)
    check("Si: máx |Σy − 1|", peor["suma_y"], 0, atol=1e-12)
    check("Si: máx |ln Q − ln K| (6 reacciones)", peor["res_K"], 0, atol=1e-9)
    check("Si: máx |ln Q − ln K| (base Si, HCl, H2)", peor["res_K_ant"], 0,
          atol=1e-9)
    check("Si: mismo gas SiCl4 R = 10 vs SiCl3H R = 7 (máx |Δy|)",
          ej.DY_EQUIV, 0, atol=1e-12)


# ---------------------------------------------------------------- caso r
def caso_r():
    head("r) Criterio de componentes: elegir_componentes (CLAUDE.md)")
    from core.independencia import (elegir_componentes, generan_mismo_espacio,
                                    matriz_atomica, matriz_reacciones,
                                    matriz_estequiometrica, reaccion_texto,
                                    reducir_por_componentes)
    ej = cargar_ej3()
    # --- Si epitaxial con alimentación SiCl4, SiCl3H, H2
    igual("Si: componentes (alimentación)", ej.COMPONENTES,
          ["SiCl4", "SiCl3H", "H2"])
    ok = "Componentes finales: SiCl4, SiCl3H, H2" in ej.CRITERIO
    print(f"  [{'PASS' if ok else 'FAIL'}] Si: el texto del criterio muestra "
          "las componentes finales")
    if not ok:
        FALLAS.append("Si criterio texto")
    igual("Si: rango", ej.rg, 3)
    igual("Si: r", ej.r, 6)
    igual("Si: columnas [I | Z]", ej.red["columnas"],
          ["SiCl4", "SiCl3H", "H2", "Si", "SiCl2", "SiCl2H2", "SiClH3",
           "SiH4", "HCl"])
    igual("Si: matriz reducida", ej.red["reducida"],
          [[1, 0, 0, -3, -1, -1, -2, -3, 1],
           [0, 1, 0, 4, 2, 2, 3, 4, -1],
           [0, 0, 1, -2, -1, 0, 0, 0, 1]])
    esperadas = [("4 SiCl3H -> Si + 3 SiCl4 + 2 H2", 1.1686e2),
                 ("2 SiCl3H -> SiCl2 + SiCl4 + H2", 6.5788e-1),
                 ("2 SiCl3H -> SiCl2H2 + SiCl4", 1.4106e-1),
                 ("3 SiCl3H -> SiClH3 + 2 SiCl4", 5.5930e-3),
                 ("4 SiCl3H -> SiH4 + 3 SiCl4", 9.2145e-5),
                 ("SiCl4 + H2 -> SiCl3H + HCl", 4.5688e-2)]
    filas, _ = matriz_reacciones([t for t, _ in esperadas], ej.ESPECIES)
    for j, ((txt, K_esp), fila) in enumerate(zip(esperadas, filas)):
        ok = ej.nu[j] == fila
        print(f"  [{'PASS' if ok else 'FAIL'}] Si: R{j + 1} esperada '{txt}' -> "
              f"ν generada: {reaccion_texto(ej.nu[j], ej.ESPECIES)}")
        if not ok:
            FALLAS.append(f"Si nueva R{j + 1}")
        igual(f"Si: K{j + 1} a 5 cifras", float(f"{ej.K[j]:.4e}"), K_esp)
    # --- equivalencia con la base anterior Si, HCl, H2 (calculada aquí)
    nu_ant = matriz_estequiometrica(
        reducir_por_componentes(ej.At, ej.ESPECIES, ["Si", "HCl", "H2"]),
        ej.ESPECIES)
    igual("Si: {SiCl4, SiCl3H, H2} y {Si, HCl, H2} generan el mismo espacio",
          generan_mismo_espacio([dict(zip(ej.ESPECIES, f)) for f in ej.nu],
                                [dict(zip(ej.ESPECIES, f)) for f in nu_ant]),
          True)

    def elegir(especies, alim, compl=None, enun=None):
        A, _ = matriz_atomica(especies)
        return elegir_componentes(especies, A, alim, compl, enun)

    ref = ["CH4", "H2O", "CO", "CO2", "H2"]
    igual("reformado: alimentación + complemento",
          elegir(ref, ["CH4", "H2O"], ["H2"])[0], ["CH4", "H2O", "H2"])
    igual("H2S: alimentación + complemento",
          elegir(["H2S", "O2", "SO2", "SO3", "H2O", "S8"], ["H2S", "O2"],
                 ["S8"])[0], ["H2S", "O2", "S8"])
    per = ["NH4ClO4", "Cl2", "N2O", "NOCl", "HCl", "H2O", "N2", "O2", "ClO2"]
    comp, txt = elegir(per, ["O2", "N2", "H2O", "HCl"], None,
                       ["NH4ClO4", "Cl2", "N2O", "NOCl"])
    igual("perclorato: las del enunciado ganan a la alimentación", comp,
          ["NH4ClO4", "Cl2", "N2O", "NOCl"])
    # metanol: CO2 = 2 CO − C, dependiente de la alimentación previa
    comp, txt = elegir(["CH3OH", "C", "CO", "CO2", "CH4", "H2O"],
                       ["C", "CO", "CO2", "H2O"])
    igual("alimentación dependiente (CO2) no es pivote", comp, ["C", "CO", "H2O"])
    ok = "AVISO (criterio 3): CO2" in txt
    print(f"  [{'PASS' if ok else 'FAIL'}] el texto avisa que CO2 es dependiente")
    if not ok:
        FALLAS.append("aviso dependiente")
    for etiqueta, args in (("alimentación", (ref, ["CH4", "XYZ"])),
                           ("complemento", (ref, ["CH4"], ["Q2"])),
                           ("enunciado", (ref, ["CH4"], None, ["CH4", "Fe"]))):
        try:
            elegir(*args)
            print(f"  [FAIL] especie inexistente en {etiqueta} NO fue rechazada")
            FALLAS.append(f"inexistente {etiqueta}")
        except ValueError as e:
            print(f"  [PASS] especie inexistente en {etiqueta} rechazada: "
                  f"{str(e)[:50]}...")
    try:
        elegir(ref, ["CH4"], ["CO2"])          # CH4 + CO2: solo 2 de rango 3
        print("  [FAIL] complemento insuficiente NO fue rechazado")
        FALLAS.append("rango no alcanzado")
    except ValueError as e:
        print(f"  [PASS] rango no alcanzado rechazado: {str(e)[:50]}...")
    comp, txt = elegir(ref, ["CH4"])
    igual("sin complemento: completa en el orden dado", comp, ["CH4", "H2O", "CO"])
    ok = "Agregadas: H2O, CO." in txt
    print(f"  [{'PASS' if ok else 'FAIL'}] el texto muestra las agregadas: "
          f"'{[l for l in txt.splitlines() if 'Agregadas' in l][0].strip()}'")
    if not ok:
        FALLAS.append("agregadas texto")


# ---------------------------------------------------------------- caso s
def caso_s():
    head("s) ν leída de [−I | −Z]: columnas = reacciones, signos directos")
    from core.independencia import (
        comprobar_conservacion, elegir_componentes, escalar_columnas,
        matriz_atomica, matriz_estequiometrica, matriz_reacciones,
        negar_renglones, nu_completa, nu_por_especie, reaccion_desde_columna,
        reaccion_texto, reducir_por_componentes, transpuesta)
    ej = cargar_ej3()

    def sistema(especies, alim, compl=None, enun=None):
        A, _ = matriz_atomica(especies)
        comp, _ = elegir_componentes(especies, A, alim, compl, enun)
        red = reducir_por_componentes(transpuesta(A), especies, comp)
        nu_frac, filas = nu_completa(red)
        nu, fac = escalar_columnas(nu_frac)
        return red, nu_frac, filas, nu, fac

    sp_ref = ["CH4", "H2O", "CO", "CO2", "H2"]
    sp_h2s = ["H2S", "O2", "SO2", "SO3", "H2O", "S8"]
    sp_per = ["NH4ClO4", "Cl2", "N2O", "NOCl", "HCl", "H2O", "N2", "O2", "ClO2"]
    casos = {
        "reformado": (sp_ref, sistema(sp_ref, ["CH4", "H2O"], ["H2"]),
                      ["CH4", "H2O", "H2", "CO", "CO2"],
                      [[-1, -1, 3, 1, 0], [-1, -2, 4, 0, 1]]),
        "H2S": (sp_h2s, sistema(sp_h2s, ["H2S", "O2"], ["S8"]),
                ["H2S", "O2", "S8", "H2O", "SO2", "SO3"],
                [[-8, -4, 1, 8, 0, 0], [0, -8, -1, 0, 8, 0],
                 [0, -12, -1, 0, 0, 8]]),
        "Si (ejercicio)": (ej.ESPECIES,
                           (ej.red, ej.nu_frac, ej.FILAS_NU, ej.NU, ej.factores),
                           ["SiCl4", "SiCl3H", "H2", "Si", "SiCl2", "SiCl2H2",
                            "SiClH3", "SiH4", "HCl"],
                           [[3, -4, 2, 1, 0, 0, 0, 0, 0],
                            [1, -2, 1, 0, 1, 0, 0, 0, 0],
                            [1, -2, 0, 0, 0, 1, 0, 0, 0],
                            [2, -3, 0, 0, 0, 0, 1, 0, 0],
                            [3, -4, 0, 0, 0, 0, 0, 1, 0],
                            [-1, 1, -1, 0, 0, 0, 0, 0, 1]]),
        "perclorato": (sp_per, sistema(sp_per, [], None,
                                       ["NH4ClO4", "Cl2", "N2O", "NOCl"]),
                       None, None),
    }
    for nombre, (especies, (red, nu_frac, filas, nu, fac), filas_esp,
                 cols_esp) in casos.items():
        cols = [list(c) for c in transpuesta(nu)]
        # --- ν de la sección 4 (renglones y columnas R1, R2, ...)
        if cols_esp is not None and filas == filas_esp:
            igual(f"{nombre}: renglones de ν", filas, filas_esp)
            for j, (c, ce) in enumerate(zip(cols, cols_esp)):
                igual(f"{nombre}: ν R{j + 1}", c, ce)
            igual(f"{nombre}: número de reacciones", len(cols), len(cols_esp))
        elif cols_esp is not None:
            print(f"  [info] {nombre}: renglones generados {', '.join(filas)}; "
                  f"esperados {', '.join(filas_esp)} (otro orden de no componentes)")
            idx = [filas.index(sp) for sp in filas_esp]
            cols_re = [[c[i] for i in idx] for c in cols]
            for j, ce in enumerate(cols_esp):
                pos = [i for i, c in enumerate(cols_re) if c == ce]
                ok = len(pos) == 1
                print(f"  [{'PASS' if ok else 'FAIL'}] {nombre}: ν R{j + 1} esperada "
                      f"{ce} -> " + (f"generada como R{pos[0] + 1} (renglones "
                                     "reordenados como los esperados)"
                                     if ok else "NO está en ν"))
                if not ok:
                    FALLAS.append(f"{nombre} ν R{j + 1}")
            igual(f"{nombre}: número de reacciones", len(cols), len(cols_esp))
        # --- Aᵀν = 0 (Aᵀ con columnas en el orden de los renglones de ν)
        igual(f"{nombre}: Aᵀν = 0",
              comprobar_conservacion(red["inicial"], transpuesta(nu))[1], True)
        # --- [−I | −Z] = −1 · [I | Z] y bloque izquierdo = −I
        k = len(red["reducida"])
        menos = negar_renglones(red["reducida"])
        igual(f"{nombre}: [−I | −Z] = −1·[I | Z]", menos,
              [[-1 * x for x in f] for f in red["reducida"]])
        igual(f"{nombre}: bloque izquierdo = −I", [f[:k] for f in menos],
              [[-int(i == j) for j in range(k)] for i in range(k)])
        # --- signo -> lado, y coeficientes leídos sin pasar términos
        malos = []
        for j, c in enumerate(cols):
            txt = reaccion_desde_columna(c, filas, especies)
            izq, der = txt.split(" -> ")
            nombres = lambda lado: {t.split()[-1] for t in lado.split(" + ")}
            if (nombres(izq) != {sp for sp, v in zip(filas, c) if v < 0}
                    or nombres(der) != {sp for sp, v in zip(filas, c) if v > 0}
                    or matriz_reacciones([txt], filas)[0][0] != c):
                malos.append(f"R{j + 1}")
        ok = not malos
        print(f"  [{'PASS' if ok else 'FAIL'}] {nombre}: ν < 0 como reactivos y "
              f"ν > 0 como productos en R1–R{len(cols)}"
              + ("" if ok else f" (fallan {', '.join(malos)})"))
        if not ok:
            FALLAS.append(f"{nombre} signos")
        # --- mismas reacciones que antes (matriz_estequiometrica + reaccion_texto)
        antes = matriz_estequiometrica(red, especies)
        igual(f"{nombre}: ν por especie = matriz_estequiometrica",
              nu_por_especie(nu, filas, especies), antes)
        igual(f"{nombre}: reacciones escritas = las anteriores",
              [reaccion_desde_columna(c, filas, especies) for c in cols],
              [reaccion_texto(f, especies) for f in antes])
    # --- H2S: fracciones en −Z, cada columna escalada a enteros (×8)
    _, (red, nu_frac, filas, nu, fac), _, _ = casos["H2S"]
    igual("H2S: ν sin escalar tiene fracciones",
          any(x.denominator > 1 for f in nu_frac for x in f), True)
    igual("H2S: factor por columna", fac, [8, 8, 8])
    igual("H2S: ν escalada es entera",
          all(isinstance(x, int) for f in nu for x in f), True)
    from math import gcd
    igual("H2S: cada columna en su mínimo entero (mcd = 1)",
          [gcd(*c) for c in transpuesta(nu)], [1, 1, 1])
    igual("H2S: el escalado conserva los signos",
          [[(x > 0) - (x < 0) for x in f] for f in nu],
          [[(x > 0) - (x < 0) for x in f] for f in nu_frac])


# ---------------------------------------------------------------- caso t
def caso_t():
    head("t) Perclorato de amonio: estequiometría por matrices (s24_ej4)")
    from fractions import Fraction as Fr
    ej = cargar_ejercicio("s24_ej4_perclorato_matrices")
    igual("perclorato: rango", ej.rg, 4)
    igual("perclorato: r", ej.r, 5)
    igual("perclorato: componentes (del enunciado)", ej.COMPONENTES,
          ["NH4ClO4", "Cl2", "N2O", "NOCl"])
    ok = "Criterio 1" in ej.CRITERIO
    print(f"  [{'PASS' if ok else 'FAIL'}] perclorato: la salida muestra el "
          "criterio 1 (componentes del enunciado)")
    if not ok:
        FALLAS.append("perclorato criterio")
    igual("perclorato: columnas [I | Z]", ej.red["columnas"],
          ["NH4ClO4", "Cl2", "N2O", "NOCl", "HCl", "H2O", "N2", "O2", "ClO2"])
    igual("perclorato: matriz reducida [I | Z]", ej.red["reducida"],
          [[1, 0, 0, 0, Fr(1, 4), Fr(1, 2), 0, 0, 0],
           [0, 1, 0, 0, Fr(5, 4), Fr(1, 2), 1, -2, Fr(-3, 2)],
           [0, 0, 1, 0, Fr(3, 4), Fr(1, 2), 2, -2, -2],
           [0, 0, 0, 1, Fr(-7, 4), Fr(-3, 2), -2, 4, 4]])
    igual("perclorato: [−I | −Z] = −1·[I | Z]", ej.menos_IZ,
          [[-x for x in f] for f in ej.red["reducida"]])
    igual("perclorato: renglones de ν", ej.FILAS_NU,
          ["NH4ClO4", "Cl2", "N2O", "NOCl", "HCl", "H2O", "N2", "O2", "ClO2"])
    cols = [list(c) for c in zip(*ej.NU)]          # una lista por reacción
    for j, ce in enumerate([[-1, -5, -3, 7, 4, 0, 0, 0, 0],
                            [-1, -1, -1, 3, 0, 2, 0, 0, 0],
                            [0, -1, -2, 2, 0, 0, 1, 0, 0],
                            [0, 2, 2, -4, 0, 0, 0, 1, 0],
                            [0, 3, 4, -8, 0, 0, 0, 0, 2]]):
        igual(f"perclorato: ν R{j + 1}", cols[j], ce)
    igual("perclorato: número de reacciones", len(cols), 5)
    for j, txt in enumerate(["NH4ClO4 + 5 Cl2 + 3 N2O -> 7 NOCl + 4 HCl",
                             "NH4ClO4 + Cl2 + N2O -> 3 NOCl + 2 H2O",
                             "Cl2 + 2 N2O -> 2 NOCl + N2",
                             "4 NOCl -> 2 Cl2 + 2 N2O + O2",
                             "8 NOCl -> 3 Cl2 + 4 N2O + 2 ClO2"]):
        igual(f"perclorato: texto R{j + 1}", ej.REACCIONES[j], txt)
    igual("perclorato: Aᵀν = 0", ej.conserva, True)


# ---------------------------------------------------------------- caso u
def caso_u():
    """Minimización de Gibbs. Lo independiente: el libro, solve_extents, la
    solución analítica de 2A = B y el modelo de GAMS (dos etapas)."""
    from math import log, sqrt

    from core.gibbs import (cargar_dgf, dGf_janaf, dGf_sva, minimizar_gibbs,
                            reacciones_independientes, validar_gibbs)

    head("u) Minimización de Gibbs (S25, core/gibbs.py)")
    T = 1000.0
    ref = ["CH4", "H2O", "CO", "CO2", "H2"]
    alim = {"CH4": 2.0, "H2O": 3.0}

    # SVA Ej. 13.14: valores impresos en el libro
    libro = {"CH4": 19720.0, "H2O": -192420.0, "CO": -200240.0,
             "CO2": -395790.0, "H2": 0.0}
    r = minimizar_gibbs(ref, libro, alim, T)
    for sp, ye in zip(ref, (0.0196, 0.0980, 0.1743, 0.0371, 0.6710)):
        check(f"SVA 13.14: y_{sp}", r.y[sp], ye, atol=1.5e-4)
    for el, le in (("C", 0.7635), ("O", 25.068), ("H", 0.1994)):
        check(f"SVA 13.14: λ_{el}/RT", r.lam_RT[el], le, atol=1e-3)

    # Contraste con Método A, con las K que salen de los mismos ΔG°f
    def contra_metodo_a(especies, dGf, n0, T, P, etiqueta):
        r = minimizar_gibbs(especies, dGf, n0, T, P)
        sp, nu, _ = reacciones_independientes(r)
        K = np.exp(-nu @ np.array([dGf[s] for s in sp]) / (R_J * T))
        m = solve_extents(sp, n0, nu, K, P)
        igual(f"{etiqueta}: Método A convergió", m.converged, True)
        for s, ya in m.mole_fractions.items():
            check(f"{etiqueta}: y_{s} vs Método A", r.y[s], ya, rtol=1e-4)

    dGf, _ = cargar_dgf(ref, T)
    contra_metodo_a(ref, dGf, alim, T, 1.0, "reformado 1 bar")
    contra_metodo_a(ref, dGf, alim, T, 10.0, "reformado 10 bar")
    baja = {sp: dGf_sva(sp, 700.0) for sp in ref}
    contra_metodo_a(ref, baja, alim, 700.0, 1.0, "reformado 700 K")

    # Steam cracker: valores del modelo GAMS NLP + CNS (s25_ej2_gibbs.gms)
    sp2 = ["CH4", "C2H4", "C2H2", "CO2", "CO", "O2", "H2", "H2O", "C2H6"]
    dg2, fuentes = cargar_dgf(sp2, T)
    igual("cracker: fuente de C2H6", fuentes["C2H6"],
          "SVA ec. 13.18 (no está en JANAF)")
    igual("cracker: fuente de CH4", fuentes["CH4"], "JANAF C-067")
    n2 = {"C2H6": 1.0, "H2O": 4.0}
    r2 = minimizar_gibbs(sp2, dg2, n2, T)
    # Método A sin O2: n_O2 ~ 1e-21 queda bajo su límite de 1e-12 y no converge
    contra_metodo_a([sp for sp in sp2 if sp != "O2"], dg2, n2, T, 1.0,
                    "cracker sin O2")
    for sp, ye in (("CH4", 7.63432e-3), ("CO2", 6.28995e-2), ("CO", 0.155081),
                   ("H2", 0.604035), ("H2O", 0.170350), ("C2H2", 3.9725e-11),
                   ("C2H4", 1.04061e-8), ("O2", 6.02693e-22)):
        check(f"cracker: y_{sp} vs GAMS", r2.y[sp], ye, rtol=1e-4)

    # La traza se fija por K: 2 H2O = 2 H2 + O2
    k = exp(-2 * (0 - dg2["H2O"]) / (R_J * T))
    check("cracker: y_O2 vs K directa", r2.y["O2"],
          k * (r2.y["H2O"] / r2.y["H2"]) ** 2, rtol=1e-8)
    check("dGf_sva(C2H6, 1000 K)", dGf_sva("C2H6", T), 110873.0, atol=1.0)
    check("dGf_sva(CH4, 1000 K) vs JANAF", dGf_sva("CH4", T),
          dGf_janaf("CH4", T)[0], atol=150.0)

    # Validación completa y arranque
    igual("cracker: validar_gibbs sin FAIL", validar_gibbs(r2).ok, True)
    x0 = {sp: 1e-8 if sp != "C2H4" else 10.0 for sp in sp2}
    r3 = minimizar_gibbs(sp2, dg2, n2, T, arranque=x0)
    check("cracker: arranque extremo, y_O2", r3.y["O2"], r2.y["O2"], rtol=1e-8)

    # P y φ: φ = 2 en todas es lo mismo que duplicar P
    a = minimizar_gibbs(ref, dGf, alim, T, P=2.0)
    b = minimizar_gibbs(ref, dGf, alim, T, phi={sp: 2.0 for sp in ref})
    check("φ uniforme = P: y_H2", b.y["H2"], a.y["H2"], rtol=1e-9)
    check("más P, menos H2 (Le Chatelier)", float(a.y["H2"] < r.y["H2"]), 1.0,
          atol=0)

    # 2 A = B con K = 4: A sin rango completo y solución analítica
    d = {"C2H4": 0.0, "C4H8": -R_J * T * log(4.0)}
    r4 = minimizar_gibbs(["C2H4", "C4H8"], d, {"C2H4": 1.0}, T)
    check("2A = B: y_A analítica", r4.y["C2H4"], (sqrt(17) - 1) / 8, rtol=1e-9)
    igual("2A = B: validar_gibbs sin FAIL", validar_gibbs(r4).ok, True)

    # Especie imposible e inerte
    d5 = dict(dGf, NH3=1.0e5, N2=0.0)
    r5 = minimizar_gibbs(ref + ["NH3"], d5, alim, T)
    igual("NH3 sin N en la carga: excluida", r5.excluidas, ["NH3"])
    igual("NH3 sin N en la carga: n = 0", r5.n["NH3"], 0.0)
    igual("con especie excluida: validar sin FAIL", validar_gibbs(r5).ok, True)
    r6 = minimizar_gibbs(ref + ["N2"], d5, dict(alim, N2=5.0), T)
    check("N2 inerte conserva sus moles", r6.n["N2"], 5.0, rtol=1e-10)

    # Entradas inválidas
    def rechaza(etiqueta, error, f):
        try:
            f()
        except error:
            igual(etiqueta, True, True)
        else:
            igual(etiqueta, False, True)

    rechaza("rechaza n0 negativo", ValueError,
            lambda: minimizar_gibbs(ref, dGf, {"CH4": -1.0}, T))
    rechaza("rechaza alimentación nula", ValueError,
            lambda: minimizar_gibbs(ref, dGf, {}, T))
    rechaza("rechaza alimentación fuera de especies", ValueError,
            lambda: minimizar_gibbs(ref, dGf, {"Ar": 1.0}, T))
    rechaza("rechaza ΔG°f faltante", KeyError,
            lambda: minimizar_gibbs(ref, {"CH4": 0.0}, alim, T))
    rechaza("pide JANAF si falta la T", KeyError,
            lambda: cargar_dgf(["CH4"], 900.0))


if __name__ == "__main__":
    for f in (caso_a, caso_b, caso_c, caso_d, caso_e, caso_f, caso_g,
              caso_h, caso_i, caso_j, caso_k, caso_l, caso_m, caso_n, caso_o,
              caso_p, caso_q, caso_r, caso_s, caso_t, caso_u):
        f()
    print(f"\n{'=' * 72}")
    if FALLAS:
        print(f"RESULTADO: {len(FALLAS)} FALLAS -> " + "; ".join(FALLAS))
        sys.exit(1)
    print("RESULTADO: TODOS los tests PASS")
    sys.exit(0)

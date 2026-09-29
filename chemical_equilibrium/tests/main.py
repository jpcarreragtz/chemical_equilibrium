"""
main.py — Runs the built-in test cases for equilibrium.py and reports
PASS/FAIL against the values solved in the course notes.

Test 1  Method B: N2O4 <-> 2 NO2, rigid (constant-volume) gas batch,
        Kc = 0.1 mol/dm3, T = 340 K, P0 = 2 atm.  Expected X_e ~ 0.44.
        (The variable-volume flow answer, X_e ~ 0.51, is also shown.)
Test 2  Method B: SO2 + 1/2 O2 -> SO3, 28% SO2 / 72% air feed,
        P = 1485 kPa, T = 500 K. Symbolic + numeric table, X = 0..1.
Test 3  Method A cross-check of Test 1 (same composition required).
Test 4  Method A: propane steam reforming, 3 reactions, T = 700/900/1000 K.

Run:  python3 main.py
"""

from __future__ import annotations

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto

from typing import Dict, List, Tuple

import numpy as np

from core.equilibrium import (
    R_ATM, ExtentResult, StoichiometricTable, Kc_to_Ky, conversion,
    format_table, run_checks_method_a, run_checks_method_b, selectivity,
    solve_extents, to_atm, verify_reactions_balanced, yield_fraction,
)
from core.selfcheck import Problem, validate, validate_sweep

# Atomic compositions of every species used in the tests (for atom checks).
ATOMS = {
    "N2O4": {"N": 2, "O": 4},
    "NO2": {"N": 1, "O": 2},
    "SO2": {"S": 1, "O": 2},
    "O2": {"O": 2},
    "SO3": {"S": 1, "O": 3},
    "N2": {"N": 2},
    "C3H8": {"C": 3, "H": 8},
    "H2O": {"H": 2, "O": 1},
    "CO": {"C": 1, "O": 1},
    "H2": {"H": 2},
    "CO2": {"C": 1, "O": 2},
    "CH4": {"C": 1, "H": 4},
}


def header(title: str) -> None:
    print("\n" + "=" * 76)
    print(title)
    print("=" * 76)


def report(label: str, value: float, expected: float,
           rel_tol: float = 0.02) -> bool:
    ok = abs(value - expected) <= rel_tol * abs(expected)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: got {value:.5g}, "
          f"expected ~{expected:.5g} (tol {100 * rel_tol:g}%)")
    return ok


# ============================================================================
# Test 1 — Method B: N2O4 <-> 2 NO2 (rigid gas batch, Kc given)
# ============================================================================

def test_1() -> Tuple[bool, float, Dict[str, float]]:
    header("TEST 1 — Method B: N2O4 <-> 2 NO2, pure N2O4 feed, "
           "T = 340 K, P0 = 2 atm, Kc = 0.1 mol/dm3")

    Kc = 0.1  # mol/dm3  (delta = 1 -> Kc has units of concentration)
    batch = StoichiometricTable(
        nu={"N2O4": -1.0, "NO2": 2.0},
        theta={"N2O4": 1.0},
        limiting="N2O4",
        T0=340.0, P0=2.0, pressure_unit="atm",
        variable_volume=False,      # rigid batch vessel: V const, P rises
        atoms=ATOMS,
    )

    print(f"  Course parameters: delta = {batch.delta:g}, "
          f"y_A0 = {batch.y_A0:g}, epsilon = y_A0*delta = {batch.epsilon:g}")
    print(f"  C_A0 = y_A0*P0/(R*T0) = {batch.C_A0:.5g} mol/L, "
          f"C_T0 = {batch.C_T0:.5g} mol/L")

    Xe = batch.solve_Xe(Kc)
    print(f"\n  Equilibrium conversion (constant volume): X_e = {Xe:.4f}")

    print("\n  Stoichiometric table at X_e "
          "(Especie | Inicial | Cambio | Remanente | Concentracion mol/L):")
    print("  " + format_table(batch.numeric_rows(Xe),
          ["Especie", "Inicial (mol)", "Cambio (mol)", "Remanente (mol)",
           "C_i (mol/L)"]).replace("\n", "\n  "))

    ok = run_checks_method_b(batch, Xe, Kc, atoms=ATOMS)

    # Same reaction fed to a FLOW reactor at constant T, P (variable volume)
    flow = StoichiometricTable(
        nu={"N2O4": -1.0, "NO2": 2.0}, theta={"N2O4": 1.0},
        limiting="N2O4", T0=340.0, P0=2.0, pressure_unit="atm",
        variable_volume=True, atoms=ATOMS)
    Xe_flow = flow.solve_Xe(Kc)
    print(f"  For comparison, flow / variable volume "
          f"(v = v0(1+eps*X)): X_e = {Xe_flow:.4f}")

    print("\n  Verification vs course notes:")
    ok &= report("C_A0 [mol/L]", batch.C_A0, 0.0717)
    ok &= report("X_e (constant-volume batch)", Xe, 0.44)

    y_eq = {sp: n / sum(batch.moles(Xe).values())
            for sp, n in batch.moles(Xe).items()}
    return ok, Xe, y_eq


# ============================================================================
# Test 2 — Method B: SO2 + 1/2 O2 -> SO3, symbolic + numeric table
# ============================================================================

def test_2() -> bool:
    header("TEST 2 — Method B: SO2 + 1/2 O2 -> SO3, 28% SO2 / 72% air, "
           "P = 1485 kPa, T = 500 K")

    y_SO2, y_O2, y_N2 = 0.28, 0.72 * 0.21, 0.72 * 0.79
    tbl = StoichiometricTable(
        nu={"SO2": -1.0, "O2": -0.5, "SO3": 1.0},
        theta={"SO2": 1.0, "O2": y_O2 / y_SO2, "N2": y_N2 / y_SO2},
        limiting="SO2",
        T0=500.0, P0=1485.0, pressure_unit="kPa",
        variable_volume=True,       # flow at constant T, P
        atoms=ATOMS,
    )

    print(f"  Feed: y_SO2 = {y_SO2}, y_O2 = {y_O2:.4g}, y_N2 = {y_N2:.4g}"
          f"  (P0 = 1485 kPa = {to_atm(1485, 'kPa'):.4g} atm)")
    print(f"  Theta_O2 = {tbl.theta_of('O2'):.4f}, "
          f"Theta_N2 (inerte) = {tbl.theta_of('N2'):.4f}")
    print(f"  delta = {tbl.delta:g}, epsilon = y_A0*delta = {tbl.epsilon:g}")
    print(f"  C_A0 = {tbl.C_A0:.5g} mol/L, C_T = P/(RT) = {tbl.C_T0:.5g} "
          f"mol/L, X_max = {tbl.X_max():.4g}")
    for sp, exc in tbl.percent_excess().items():
        print(f"  % excess of {sp}: {exc:.2f} %")

    print("\n  Symbolic stoichiometric table (per mol of A = SO2):")
    print("  " + format_table(tbl.symbolic_rows(),
          ["Especie", "Inicial", "Cambio", "Remanente", "Concentracion"]
          ).replace("\n", "\n  "))

    print("\n  Numeric table, C_i(X) in mol/L for X = 0 .. 1:")
    rows: List[List[float]] = []
    for X in np.linspace(0.0, 1.0, 11):
        C = tbl.concentrations(X)
        rows.append([X, C["SO2"], C["O2"], C["SO3"], C["N2"], tbl.C_T(X)])
    print("  " + format_table(rows,
          ["X", "C_SO2", "C_O2", "C_SO3", "C_N2", "C_T"]
          ).replace("\n", "\n  "))

    ok = run_checks_method_b(tbl, 1.0, Kc=None, atoms=ATOMS)

    print("\n  Verification vs course notes:")
    ok &= report("Theta_O2", tbl.theta_of("O2"), 0.54)
    ok &= report("Theta_N2 (inerte)", tbl.theta_of("N2"), 2.04)
    ok &= report("delta", tbl.delta, -0.5)
    ok &= report("epsilon", tbl.epsilon, -0.14)
    ok &= report("C_A0 [mol/L]", tbl.C_A0, 0.1)
    ok &= report("C_T [mol/L]", tbl.C_T0, 0.357)
    return ok


# ============================================================================
# Test 3 — Method A sanity cross-check of Test 1
# ============================================================================

def test_3(Xe_ref: float, y_ref: Dict[str, float]) -> bool:
    header("TEST 3 — Method A cross-check: N2O4 <-> 2 NO2 by extent of "
           "reaction (must match Test 1)")

    # Convert Kc -> mole-fraction K:  K = Kc / C_T0^delta = Kc*(R*T)^delta
    # with P0 = 1 atm.  K = 0.1 * (0.08206*340) = 2.7900.
    T, delta = 340.0, 1.0
    K = Kc_to_Ky(0.1, delta, T0=T, P0=1.0, unit="atm")
    print(f"  K = Kc/C_T0^delta = Kc*(R*T/P0)^delta = {K:.5f}  (P0 = 1 atm)")

    # Test 1 is a RIGID vessel: as moles grow, P grows: P = P_init*n_tot/n0.
    # Method A takes P as an input, so iterate P to self-consistency.
    species, n0 = ["N2O4", "NO2"], {"N2O4": 1.0}
    nu = [[-1.0, 2.0]]
    P = 2.0
    res: ExtentResult = None  # type: ignore[assignment]
    for it in range(100):
        res = solve_extents(species, n0, nu, [K], P=P, P0=1.0, atoms=ATOMS)
        P_new = 2.0 * res.n_total / 1.0
        if abs(P_new - P) < 1e-12:
            break
        P = P_new
    print(f"  Self-consistent vessel pressure: P = {P:.5f} atm "
          f"(= 2 atm * n_total/n0, after {it + 1} iterations)")

    X = conversion(1.0, res.moles["N2O4"])
    print(f"  xi = {res.xi[0]:.5f} mol  ->  X = (n_A0-n_A)/n_A0 = {X:.4f}")
    print("  " + format_table(
        [[sp, res.moles[sp], res.mole_fractions[sp]] for sp in species],
        ["Especie", "n_i (mol)", "y_i"]).replace("\n", "\n  "))

    ok = run_checks_method_a(res, n0, atoms=ATOMS)

    problem = Problem(species=species, n0=n0, nu=nu, K=[K], P=P, P0=1.0,
                      atoms=ATOMS, T=T, pressure_unit="atm", label="Test 3")
    ok &= validate(problem, res).ok

    print("\n  Verification vs Test 1 (Method B):")
    ok &= report("X_e", X, Xe_ref, rel_tol=1e-4)
    for sp in species:
        ok &= report(f"y_{sp}", res.mole_fractions[sp], y_ref[sp],
                     rel_tol=1e-4)
    return ok


# ============================================================================
# Test 4 — Method A: propane steam reforming, 3 reactions, 3 temperatures
# ============================================================================

def test_4() -> bool:
    header("TEST 4 — Method A: propane steam reforming (3 reactions), "
           "feed 1 C3H8 + 4 H2O + 0.5 N2, P = 1 bar")

    species = ["C3H8", "H2O", "CO", "H2", "CO2", "CH4", "N2"]
    nu = [
        [-1.0, -3.0, 3.0, 7.0, 0.0, 0.0, 0.0],   # R1: C3H8+3H2O <-> 3CO+7H2
        [0.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0],   # R2: CO+H2O <-> CO2+H2
        [-1.0, 0.0, 0.0, -2.0, 0.0, 3.0, 0.0],   # R3: C3H8+2H2 <-> 3CH4
    ]
    n0 = {"C3H8": 1.0, "H2O": 4.0, "N2": 0.5}
    K_table = {  # T [K] -> (K1, K2, K3)
        700.0: (7.2285e-3, 7.5034, 2.4203e10),
        900.0: (1.3407e6, 1.5591, 2.4864e8),
        1000.0: (1.0501e9, 0.89956, 5.0080e7),
    }

    problems = verify_reactions_balanced(np.array(nu), species, ATOMS)
    if problems:
        print("  Atom checker found unbalanced reactions:")
        for p in problems:
            print("   ", p)
        print("  -> STOPPING this test (fix stoichiometry first).")
        return False
    print("  Atom checker: R1, R2, R3 are all balanced in C, H and O "
          "as given — nothing changed.")

    n_A0 = n0["C3H8"]
    # H2 ceiling: R1+R2 per mol C3H8 give 10 H2 but need 6 H2O; with only
    # 4 H2O the TRUE feed cap is the elemental-H limit,
    # (8*n_C3H8 + 2*n_H2O)/2 = 8 mol, reached by R1 to completion plus WGS
    # on the leftover water. It is NOT (10/6)*4 = 6.667: that splits the
    # water pro-rata, but reforming yields 7/3 H2 per H2O vs only 1 for WGS.
    h2_max_propane = 10.0 * n_A0
    h2_max_water = (8.0 * n_A0 + 2.0 * n0["H2O"]) / 2.0
    print(f"  Theoretical max H2: {h2_max_propane:g} mol (per propane, "
          f"R1+R2 complete) | {h2_max_water:.4g} mol (elemental-H cap of "
          f"this feed)")

    ok = True
    columns: Dict[float, List[float]] = {}
    xi_by_T: Dict[float, Tuple[float, ...]] = {}
    for T, K in K_table.items():
        print(f"\n  --- T = {T:g} K,  K = {K} ---")
        res = solve_extents(species, n0, nu, list(K), P=1.0, P0=1.0,
                            atoms=ATOMS)
        ok &= run_checks_method_a(res, n0, atoms=ATOMS)

        problem = Problem(species=species, n0=n0, nu=nu, K=list(K), P=1.0,
                          P0=1.0, atoms=ATOMS, inerts=["N2"], T=T,
                          pressure_unit="bar", label="Test 4")
        ok &= validate(problem, res).ok
        xi_by_T[T] = tuple(res.xi)

        n = res.moles
        X = conversion(n_A0, n["C3H8"])
        h2_yield_th = yield_fraction(n["H2"], n_A0, 10.0)
        h2_yield_water = n["H2"] / h2_max_water
        columns[T] = ([*res.xi, *res.y, X,
                       selectivity(n["CO"], n["CO2"]),
                       selectivity(n["CO"], n["CH4"]),
                       selectivity(n["CH4"], n["CO2"]),
                       h2_yield_th, h2_yield_water])

    labels = (["xi_1 (mol)", "xi_2 (mol)", "xi_3 (mol)"]
              + [f"y_{sp}" for sp in species]
              + ["X_C3H8", "S CO/CO2", "S CO/CH4", "S CH4/CO2",
                 "H2 yield vs 10/C3H8", "H2 yield vs water-max"])
    temps = sorted(columns)
    rows = [[lab] + [columns[T][k] for T in temps]
            for k, lab in enumerate(labels)]
    print("\n  Comparison across temperatures "
          "(selectivity S i/j = mol i / mol j among C-containing products):")
    print("  " + format_table(rows, ["Quantity"] +
          [f"T = {T:g} K" for T in temps]).replace("\n", "\n  "))

    ok &= validate_sweep(K_table, xi_by_T).ok

    print(f"\n  [{'PASS' if ok else 'FAIL'}] all equilibrium/atom/"
          f"non-negativity checks at the three temperatures")
    return ok


# ============================================================================

def main() -> int:
    print("Chemical equilibrium test suite "
          f"(R = {R_ATM} L·atm/mol/K; C_T = P/(RT) checks in mol/L)")
    results: Dict[str, bool] = {}

    ok1, Xe1, y1 = test_1()
    results["Test 1 (Method B, N2O4)"] = ok1
    results["Test 2 (Method B, SO2 oxidation table)"] = test_2()
    results["Test 3 (Method A cross-check of Test 1)"] = test_3(Xe1, y1)
    results["Test 4 (Method A, propane reforming)"] = test_4()

    header("SUMMARY")
    for name, ok in results.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    n_fail = sum(not ok for ok in results.values())
    print(f"\n  {len(results) - n_fail}/{len(results)} tests passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())

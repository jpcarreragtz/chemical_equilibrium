# Chemical Equilibrium (gas phase, ideal gas)

Two solution frameworks matching the two notations used in the course:

- **Method A — extent of reaction ξ** (Smith, Van Ness & Abbott), for 1..N
  simultaneous reactions:
  `K_j = Π y_i^ν_ij (P/P0)^δ_j`, solved with `scipy.optimize.least_squares`
  on log-form residuals (handles K from 1e-3 to 1e10), analytic Jacobian,
  multi-start over feasible extents, all `n_i ≥ 0` enforced.
- **Method B — stoichiometric table with conversion X** (Fogler), single
  reaction `A + (b/a)B -> (c/a)C + (d/a)D`, batch or flow, constant or
  variable volume. Builds Θ_i, δ, ε = y_A0·δ, the full table, %excess,
  and solves `Kc = Π C_i^ν_i` for X_e with `brentq` on `[0, X_max]`.

## Files

- `equilibrium.py` — core module: both solvers, unit handling
  (atm/bar/kPa, R = 0.08206 L·atm/mol/K), `Kc <-> K` conversion,
  element-balance verification, conversion/selectivity/yield, and the
  automatic "checks" report (atom balance in vs out, `n_i ≥ 0`,
  `|K_calc−K|/K < 1e-6`, `Σy_i = 1` / `C_T = P/(RT)` constant).
- `main.py` — runs the four built-in course test cases with PASS/FAIL:
  1. Method B: N2O4 ⇌ 2 NO2, rigid batch, Kc = 0.1 mol/dm³ → X_e ≈ 0.44
     (flow/variable-volume answer ≈ 0.51 also shown)
  2. Method B: SO2 + ½O2 → SO3, 28% SO2/72% air, 1485 kPa, 500 K —
     symbolic + numeric table for X = 0..1
  3. Method A cross-check of Test 1 (same composition, self-consistent
     vessel pressure P = P0·n_tot/n0)
  4. Method A: propane steam reforming (3 reactions) at 700/900/1000 K,
     with ξ_j, y_i, conversion, C-selectivities and H2 yields

## Run

```bash
python3 main.py        # needs numpy + scipy (tabulate optional, cosmetic)
```

Exit code 0 = all tests pass.

## Quick use

```python
from equilibrium import StoichiometricTable, solve_extents, Kc_to_Ky

# Method B (Fogler): reaction normalized per mol of limiting reactant A
tbl = StoichiometricTable(
    nu={"N2O4": -1, "NO2": 2},      # ν per mol A (A must be -1)
    theta={"N2O4": 1.0},            # Θ_i = n_i0/n_A0 for the feed
    limiting="N2O4",
    T0=340, P0=2.0, pressure_unit="atm",
    variable_volume=False,          # True = flow / constant-P (v = v0(1+εX))
)
Xe = tbl.solve_Xe(Kc=0.1)           # mol/L units consistent with δ

# Method A (Smith Van Ness): K_j = Π y_i^ν (P/P0)^δ
res = solve_extents(
    species=["N2O4", "NO2"], n0={"N2O4": 1.0},
    nu=[[-1, 2]], K=[2.79], P=2.883, P0=1.0,
    atoms={"N2O4": {"N": 2, "O": 4}, "NO2": {"N": 1, "O": 2}},
)
print(res.xi, res.mole_fractions)
```

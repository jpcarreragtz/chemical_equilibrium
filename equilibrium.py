"""
equilibrium.py — Gas-phase chemical equilibrium (ideal gas) for a Chemical
Equilibrium course.

Two solution frameworks are implemented, matching the two notations used in
the course:

METHOD A — Extent of reaction xi (Smith, Van Ness & Abbott notation),
for 1..N simultaneous reactions:

    n_i     = n_i0 + sum_j nu_ij * xi_j
    n_total = n0   + sum_j delta_j * xi_j ,   delta_j = sum_i nu_ij
    y_i     = n_i / n_total
    K_j     = prod_i y_i^nu_ij * (P/P0)^delta_j        (ideal gas)

    Solved with scipy.optimize.least_squares on LOG-form residuals

        r_j = sum_i nu_ij*ln(y_i) + delta_j*ln(P/P0) - ln(K_j)

    which behaves well for K spanning ~1e-3 to 1e10. Feasibility
    (all n_i >= 0) is enforced with penalty residuals plus multi-start.

METHOD B — Stoichiometric table with conversion X of the limiting reactant A
(Fogler notation), single reaction, batch or flow:

    A + (b/a) B -> (c/a) C + (d/a) D          (reversible allowed)

    Theta_i = n_i0/n_A0   (= F_i0/F_A0 = C_i0/C_A0 = y_i0/y_A0 for flow)
    delta   = d/a + c/a - b/a - 1
    epsilon = y_A0 * delta

    n_A = n_A0 (1 - X)          n_B = n_A0 (Theta_B - (b/a) X)
    n_C = n_A0 (Theta_C + (c/a) X)   n_D = n_A0 (Theta_D + (d/a) X)
    inerts: n_I = n_A0 * Theta_I

    Concentrations:
      constant volume (liquid, or gas in a rigid batch reactor):
          C_i = C_A0 (Theta_i + nu_i X)
      gas, variable volume (flow reactor or constant-P batch):
          v   = v0 (1 + epsilon X)(P0/P)(T/T0)
          C_i = C_A0 (Theta_i + nu_i X)/(1 + epsilon X) * (T0/T)(P/P0)

    Equilibrium conversion X_e from  Kc = prod_i C_i^nu_i  with a bounded
    scalar solver (scipy.optimize.brentq) on [0, X_max], where X_max is the
    largest X keeping every species non-negative.

Shared utilities: unit handling (atm/bar/kPa, R = 0.08206 L·atm/mol/K),
element balance verification, Kc <-> mole-fraction-K conversion, conversion /
selectivity / yield, table formatting, and an automatic "checks" report
(atom balance in vs out, non-negativity, residuals, sum(y)=1 or C_T constant).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq, least_squares

# ----------------------------------------------------------------------------
# Constants & unit handling
# ----------------------------------------------------------------------------

R_ATM: float = 0.08206  # L·atm / (mol·K)

_PRESSURE_TO_ATM = {
    "atm": 1.0,
    "bar": 1.0 / 1.01325,
    "kpa": 1.0 / 101.325,
    "kPa": 1.0 / 101.325,
}


def to_atm(P: float, unit: str = "atm") -> float:
    """Convert a pressure in atm / bar / kPa to atm."""
    key = unit if unit in _PRESSURE_TO_ATM else unit.lower()
    if key not in _PRESSURE_TO_ATM:
        raise ValueError(f"Unknown pressure unit '{unit}' (use atm, bar, kPa)")
    return P * _PRESSURE_TO_ATM[key]


def total_concentration(P: float, T: float, unit: str = "atm") -> float:
    """C_T = P/(R T) in mol/L for an ideal gas (P in atm/bar/kPa, T in K)."""
    return to_atm(P, unit) / (R_ATM * T)


def Kc_to_Ky(Kc: float, delta: float, T0: float, P0: float = 1.0,
             unit: str = "atm") -> float:
    """
    Convert a concentration-based Kc (mol/L units) to the mole-fraction /
    activity K used in Method A:

        K = Kc / C_T0^delta ,   C_T0 = P0 / (R T0)

    (equivalently K = Kc (R T0 / P0)^delta). Consistent with
    K = prod y_i^nu_i (P/P0)^delta evaluated at the reference P0.
    """
    C_T0 = total_concentration(P0, T0, unit)
    return Kc / C_T0 ** delta


def Ky_to_Kc(K: float, delta: float, T0: float, P0: float = 1.0,
             unit: str = "atm") -> float:
    """Inverse of :func:`Kc_to_Ky`:  Kc = K * C_T0^delta."""
    C_T0 = total_concentration(P0, T0, unit)
    return K * C_T0 ** delta


# ----------------------------------------------------------------------------
# Table formatting (tabulate if available, plain-text fallback)
# ----------------------------------------------------------------------------

try:  # pragma: no cover - cosmetic only
    from tabulate import tabulate as _tabulate

    def format_table(rows: Sequence[Sequence], headers: Sequence[str]) -> str:
        return _tabulate(rows, headers=headers, tablefmt="github",
                         floatfmt=".5g")
except ImportError:  # plain fallback, no extra dependency needed

    def format_table(rows: Sequence[Sequence], headers: Sequence[str]) -> str:
        def fmt(x) -> str:
            if isinstance(x, (int, np.integer)):
                return str(x)
            if isinstance(x, (float, np.floating)):
                return f"{x:.5g}"
            return str(x)

        str_rows = [[fmt(c) for c in row] for row in rows]
        widths = [max(len(h), *(len(r[k]) for r in str_rows)) if str_rows
                  else len(h) for k, h in enumerate(headers)]
        line = "| " + " | ".join(h.ljust(w) for h, w in zip(headers, widths)) + " |"
        sep = "|-" + "-|-".join("-" * w for w in widths) + "-|"
        body = ["| " + " | ".join(c.ljust(w) for c, w in zip(r, widths)) + " |"
                for r in str_rows]
        return "\n".join([line, sep, *body])


# ----------------------------------------------------------------------------
# Element (atom) balance utilities
# ----------------------------------------------------------------------------

AtomComposition = Dict[str, Dict[str, float]]  # e.g. {"C3H8": {"C": 3, "H": 8}}


def reaction_element_imbalance(nu_row: Sequence[float], species: Sequence[str],
                               atoms: AtomComposition) -> Dict[str, float]:
    """
    Net atoms produced by one reaction: sum_i nu_i * (atoms in species i),
    per element. A balanced reaction returns ~0 for every element.
    """
    net: Dict[str, float] = {}
    for nu_i, sp in zip(nu_row, species):
        if abs(nu_i) < 1e-14:
            continue
        if sp not in atoms:
            raise KeyError(f"No atomic composition given for species '{sp}'")
        for el, count in atoms[sp].items():
            net[el] = net.get(el, 0.0) + nu_i * count
    return net


def verify_reactions_balanced(nu: np.ndarray, species: Sequence[str],
                              atoms: AtomComposition,
                              tol: float = 1e-9) -> List[str]:
    """
    Element-wise balance check of every reaction (row of nu). Returns a list
    of human-readable problem descriptions; empty list = all balanced.
    """
    problems: List[str] = []
    for j, row in enumerate(np.atleast_2d(np.asarray(nu, dtype=float))):
        net = reaction_element_imbalance(row, species, atoms)
        bad = {el: v for el, v in net.items() if abs(v) > tol}
        if bad:
            problems.append(
                f"Reaction {j + 1} is NOT balanced: net atoms {bad}")
    return problems


def atom_totals(moles: Dict[str, float], atoms: AtomComposition) -> Dict[str, float]:
    """Total mol of each element contained in a {species: moles} mixture."""
    tot: Dict[str, float] = {}
    for sp, n in moles.items():
        for el, count in atoms.get(sp, {}).items():
            tot[el] = tot.get(el, 0.0) + n * count
    return tot


# ----------------------------------------------------------------------------
# Shared performance metrics
# ----------------------------------------------------------------------------

def conversion(n_A0: float, n_A: float) -> float:
    """X = (n_A0 - n_A)/n_A0 for the user-specified limiting reactant A."""
    return (n_A0 - n_A) / n_A0


def selectivity(n_desired: float, n_undesired: float) -> float:
    """S = mol desired product / mol undesired product (inf if none formed)."""
    return np.inf if n_undesired <= 0 else n_desired / n_undesired


def yield_fraction(n_desired: float, n_A0: float,
                   coef_desired_per_A: float) -> float:
    """
    Rendimiento: mol desired product / mol that WOULD form if ALL of the
    limiting reactant reacted through the desired reaction only
    (= n_A0 * coef_desired_per_A).
    """
    return n_desired / (n_A0 * coef_desired_per_A)


# ----------------------------------------------------------------------------
# METHOD A — extent of reaction, multi-reaction
# ----------------------------------------------------------------------------

@dataclass
class ExtentResult:
    """Solution of a Method-A (extent of reaction) equilibrium problem."""
    species: List[str]
    xi: np.ndarray            # extents xi_j, one per reaction [mol]
    n: np.ndarray             # equilibrium moles n_i [mol]
    y: np.ndarray             # equilibrium mole fractions
    n_total: float
    K_given: np.ndarray
    K_calc: np.ndarray        # K recomputed from the solution
    P: float
    P0: float
    converged: bool
    message: str = ""

    @property
    def moles(self) -> Dict[str, float]:
        return dict(zip(self.species, self.n))

    @property
    def mole_fractions(self) -> Dict[str, float]:
        return dict(zip(self.species, self.y))


def solve_extents(species: Sequence[str],
                  n0: Dict[str, float],
                  nu: Sequence[Sequence[float]],
                  K: Sequence[float],
                  P: float,
                  P0: float = 1.0,
                  atoms: Optional[AtomComposition] = None,
                  n_starts: int = 40,
                  seed: int = 0) -> ExtentResult:
    """
    Solve the multi-reaction equilibrium (Method A, Smith Van Ness).

    Parameters
    ----------
    species : names of all species (inerts included).
    n0      : initial moles per species (absent species -> 0; inerts have
              nu = 0 in every reaction).
    nu      : stoichiometric matrix, shape (n_reactions, n_species);
              nu[j][i] < 0 reactant, > 0 product, 0 inert/not involved.
    K       : equilibrium constants K_j (mole-fraction form with (P/P0)^delta,
              i.e. K_j = prod y_i^nu_ij (P/P0)^delta_j).
    P, P0   : system pressure and reference pressure, SAME units (default
              P0 = 1, e.g. 1 bar; only the ratio P/P0 enters).
    atoms   : optional atomic compositions; if given, every reaction is
              verified element-balanced BEFORE solving (raises otherwise).
    n_starts: multi-start attempts if the first solve fails.

    Returns
    -------
    ExtentResult with xi_j, n_i, y_i and recomputed K for validation.

    Notes
    -----
    Residuals are solved in log form,
        r_j = sum_i nu_ij ln y_i + delta_j ln(P/P0) - ln K_j,
    so K from 1e-3 to 1e10 is well conditioned. Moles inside the logs are
    clamped at 1e-12 and negative moles add large penalty residuals, which
    (with multi-start over random feasible extents) keeps n_i >= 0.
    """
    species = list(species)
    nu = np.atleast_2d(np.asarray(nu, dtype=float))
    n_rxn, n_sp = nu.shape
    if n_sp != len(species):
        raise ValueError("nu has a different number of columns than species")
    K = np.asarray(K, dtype=float)
    if K.shape != (n_rxn,):
        raise ValueError("need exactly one K per reaction")
    if np.any(K <= 0):
        raise ValueError("all K_j must be > 0")

    if atoms is not None:
        problems = verify_reactions_balanced(nu, species, atoms)
        if problems:
            raise ValueError("Unbalanced reaction(s):\n" + "\n".join(problems))

    n0_vec = np.array([float(n0.get(sp, 0.0)) for sp in species])
    if np.any(n0_vec < 0):
        raise ValueError("initial moles must be non-negative")
    delta = nu.sum(axis=1)
    lnK = np.log(K)
    ln_Pr = np.log(P / P0)
    EPS = 1e-12
    PEN = 1e6

    def residuals(xi: np.ndarray) -> np.ndarray:
        n = n0_vec + nu.T @ xi
        penalty = PEN * np.minimum(n, 0.0)          # pushes back to n_i >= 0
        n_pos = np.maximum(n, EPS)
        y = n_pos / n_pos.sum()
        r = nu @ np.log(y) + delta * ln_Pr - lnK
        return np.concatenate([r, penalty])

    def jacobian(xi: np.ndarray) -> np.ndarray:
        # r_j = sum_i nu_ij ln n_i - delta_j ln n_tot + const, hence
        #   dr_j/dxi_k = sum_i nu_ij nu_ik / n_i - delta_j delta_k / n_tot.
        # Analytic form is essential: near equilibrium a species can sit at
        # ~1e-10 mol, far below any finite-difference step size, which makes
        # numeric jacobians (and hence convergence) fail.
        n = n0_vec + nu.T @ xi
        n_pos = np.maximum(n, EPS)
        J_K = (nu / n_pos) @ nu.T - np.outer(delta, delta) / n_pos.sum()
        J_pen = PEN * np.where(n[:, None] < 0.0, nu.T, 0.0)
        return np.vstack([J_K, J_pen])

    rng = np.random.default_rng(seed)
    n_feed = max(float(n0_vec.sum()), 1.0)

    def sample_start(frac: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Build a feasible extent vector by advancing the reactions ONE AT A
        TIME (random order), each within the bounds allowed by the moles
        left after the previous ones. In coupled networks (a reactant of one
        reaction only produced by another) this always yields feasible
        samples, where independent box-sampling rejects almost everything.
        """
        xi = np.zeros(n_rxn)
        n = n0_vec.copy()
        order = rng.permutation(n_rxn) if frac is None else range(n_rxn)
        for j in order:
            cons = np.where(nu[j] < 0)[0]
            prod = np.where(nu[j] > 0)[0]
            hi = min([n[i] / -nu[j, i] for i in cons], default=n_feed)
            lo = -min([n[i] / nu[j, i] for i in prod], default=n_feed)
            f = rng.random() if frac is None else float(frac[j])
            xi[j] = lo + f * (hi - lo)
            n = n + nu[j] * xi[j]
        return xi

    starts: List[np.ndarray] = [np.zeros(n_rxn),
                                sample_start(np.full(n_rxn, 0.75)),
                                sample_start(np.full(n_rxn, 0.55))]
    while len(starts) < n_starts:
        starts.append(sample_start())
    scale = np.maximum(np.max(np.abs(np.array(starts)), axis=0), 1e-3)

    best: Optional[Tuple[float, np.ndarray]] = None
    for x0 in starts:
        try:
            sol = least_squares(residuals, x0, jac=jacobian, method="trf",
                                x_scale=scale, xtol=1e-15, ftol=1e-15,
                                gtol=1e-15, max_nfev=2000)
        except Exception:
            continue
        n = n0_vec + nu.T @ sol.x
        r_K = residuals(sol.x)[:n_rxn]
        err = float(np.max(np.abs(r_K)))
        # Acceptance at 5e-7 in log space (tighter than the 1e-6
        # relative-K check): a trace species computed as n0 - sum(nu*xi)
        # has an absolute cancellation floor of ~1e-16 mol, so for
        # n_i ~ 1e-10 mol the log-residual cannot go below ~1e-7.
        ok = err < 5e-7 and float(np.min(n)) > -1e-9
        if best is None or err < best[0]:
            best = (err, sol.x)
        if ok:
            break

    assert best is not None
    err, xi = best
    n = np.maximum(n0_vec + nu.T @ xi, 0.0)
    n_total = float(n.sum())
    y = n / n_total
    y_safe = np.maximum(y, EPS)
    K_calc = np.exp(nu @ np.log(y_safe) + delta * ln_Pr)
    converged = err < 5e-7
    msg = ("converged" if converged else
           f"NOT converged: max |log-residual| = {err:.3e} after "
           f"{len(starts)} starts")
    return ExtentResult(species=species, xi=xi, n=n, y=y, n_total=n_total,
                        K_given=K, K_calc=K_calc, P=P, P0=P0,
                        converged=converged, message=msg)


# ----------------------------------------------------------------------------
# METHOD B — stoichiometric table with conversion X (single reaction)
# ----------------------------------------------------------------------------

@dataclass
class StoichiometricTable:
    """
    Fogler-style stoichiometric table for a single reaction normalized to
    the limiting reactant A:

        A + (b/a) B -> (c/a) C + (d/a) D      (inerts I allowed)

    Parameters
    ----------
    nu : per-mol-of-A stoichiometric numbers {species: nu_i}; the limiting
         reactant MUST have nu = -1 (i.e. already divided through by a).
         Reactants < 0, products > 0. Species with nu = 0 (or absent from
         nu but present in theta) are inerts.
    theta : Theta_i = n_i0/n_A0 for every species present in the FEED
            (Theta_A = 1). Species that appear only as products may be
            omitted (Theta = 0).
    limiting : name of species A.
    T0, P0 : feed (entering) temperature [K] and pressure.
    T, P   : reaction conditions; default = T0, P0 (isothermal, isobaric).
    pressure_unit : 'atm' | 'bar' | 'kPa' for P0 and P.
    variable_volume : True  -> gas at constant T,P (flow / piston batch):
                               v = v0 (1+eps X)(P0/P)(T/T0)
                      False -> constant volume (liquid, or rigid gas batch
                               where P rises instead).
    n_A0 : basis, mol of A fed (mol/s for flow); default 1.

    Derived attributes: y_A0 = 1/sum(Theta_i), delta = sum nu_i,
    epsilon = y_A0*delta, C_A0 = y_A0*P0/(R T0), C_T0 = P0/(R T0).
    """
    nu: Dict[str, float]
    theta: Dict[str, float]
    limiting: str
    T0: float
    P0: float
    T: Optional[float] = None
    P: Optional[float] = None
    pressure_unit: str = "atm"
    variable_volume: bool = True
    n_A0: float = 1.0
    atoms: Optional[AtomComposition] = None

    species: List[str] = field(init=False)

    def __post_init__(self) -> None:
        if self.limiting not in self.nu or not np.isclose(self.nu[self.limiting], -1.0):
            raise ValueError("limiting reactant must have nu = -1 "
                             "(normalize the reaction per mol of A)")
        if self.theta.get(self.limiting) != 1.0:
            raise ValueError("Theta of the limiting reactant must be 1")
        self.T = self.T0 if self.T is None else self.T
        self.P = self.P0 if self.P is None else self.P
        # keep feed order: species in theta first, then product-only species
        self.species = list(self.theta) + [s for s in self.nu
                                           if s not in self.theta]
        if self.atoms is not None:
            row = [self.nu.get(sp, 0.0) for sp in self.species]
            problems = verify_reactions_balanced(
                np.array([row]), self.species, self.atoms)
            if problems:
                raise ValueError("Unbalanced reaction:\n" + "\n".join(problems))

    # -- course parameters ---------------------------------------------------
    def nu_of(self, sp: str) -> float:
        return self.nu.get(sp, 0.0)

    def theta_of(self, sp: str) -> float:
        return self.theta.get(sp, 0.0)

    @property
    def y_A0(self) -> float:
        """Feed mole fraction of A: 1/sum(Theta_i) since Theta_i = y_i0/y_A0."""
        return 1.0 / sum(self.theta.values())

    @property
    def delta(self) -> float:
        """delta = (d/a + c/a) - (b/a + 1) = sum_i nu_i (per mol A)."""
        return sum(self.nu.values())

    @property
    def epsilon(self) -> float:
        """epsilon = y_A0 * delta (fractional volume change at X = 1)."""
        return self.y_A0 * self.delta

    @property
    def C_T0(self) -> float:
        """Feed total concentration P0/(R T0) [mol/L]."""
        return total_concentration(self.P0, self.T0, self.pressure_unit)

    @property
    def C_A0(self) -> float:
        """C_A0 = y_A0 P0/(R T0) [mol/L]."""
        return self.y_A0 * self.C_T0

    # -- table entries -------------------------------------------------------
    def moles(self, X: float) -> Dict[str, float]:
        """n_i(X) = n_A0 (Theta_i + nu_i X); inerts unchanged."""
        return {sp: self.n_A0 * (self.theta_of(sp) + self.nu_of(sp) * X)
                for sp in self.species}

    def concentrations(self, X: float) -> Dict[str, float]:
        """
        C_i(X) [mol/L].

        variable volume (gas, const T,P):
            C_i = C_A0 (Theta_i + nu_i X)/(1+eps X) * (T0/T)(P/P0)
        constant volume:
            C_i = C_A0 (Theta_i + nu_i X)
        """
        if self.variable_volume:
            Pr = (to_atm(self.P, self.pressure_unit)
                  / to_atm(self.P0, self.pressure_unit))
            f = (self.T0 / self.T) * Pr / (1.0 + self.epsilon * X)
        else:
            f = 1.0
        return {sp: self.C_A0 * (self.theta_of(sp) + self.nu_of(sp) * X) * f
                for sp in self.species}

    def C_T(self, X: float) -> float:
        """Total concentration at conversion X (should equal P/(RT) for the
        variable-volume gas case at fixed T, P)."""
        return sum(self.concentrations(X).values())

    def X_max(self) -> float:
        """
        Largest X keeping every species non-negative:
        X_max = min(1, Theta_i/|nu_i|) over reactants (detects a reactant in
        deficit relative to A, e.g. Theta_B < b/a).
        """
        limits = [1.0]
        for sp in self.species:
            nu_i = self.nu_of(sp)
            if nu_i < 0 and sp != self.limiting:
                limits.append(self.theta_of(sp) / -nu_i)
        return min(limits)

    def percent_excess(self) -> Dict[str, float]:
        """% excess of each non-limiting reactant:
        (Theta_B - b/a)/(b/a) * 100."""
        out = {}
        for sp in self.species:
            nu_i = self.nu_of(sp)
            if nu_i < 0 and sp != self.limiting:
                coef = -nu_i
                out[sp] = (self.theta_of(sp) - coef) / coef * 100.0
        return out

    # -- equilibrium ---------------------------------------------------------
    def Kc_calc(self, X: float) -> float:
        """Kc(X) = prod_i C_i^nu_i (units of (mol/L)^delta)."""
        C = self.concentrations(X)
        return float(np.prod([C[sp] ** self.nu_of(sp)
                              for sp in self.species if self.nu_of(sp) != 0]))

    def solve_Xe(self, Kc: float) -> float:
        """
        Equilibrium conversion X_e from Kc = prod C_i^nu_i, using brentq on
        the log-form residual g(X) = sum nu_i ln C_i - ln Kc over
        [~0, ~X_max]. g -> -inf as products vanish (X->0 with no products in
        the feed) and +inf as a reactant is exhausted (X->X_max), so a root
        is bracketed whenever equilibrium lies inside (0, X_max).
        """
        if Kc <= 0:
            raise ValueError("Kc must be > 0")

        def g(X: float) -> float:
            C = self.concentrations(X)
            tot = 0.0
            for sp in self.species:
                nu_i = self.nu_of(sp)
                if nu_i != 0:
                    tot += nu_i * np.log(max(C[sp], 1e-300))
            return tot - np.log(Kc)

        lo, hi = 1e-12, self.X_max() - 1e-12
        g_lo, g_hi = g(lo), g(hi)
        if g_lo * g_hi > 0:
            raise RuntimeError(
                f"No equilibrium in (0, X_max={self.X_max():.4g}): "
                f"g(0)={g_lo:.3g}, g(X_max)={g_hi:.3g}. "
                "Check Kc units/feed (products already past equilibrium?).")
        return float(brentq(g, lo, hi, xtol=1e-14, rtol=1e-14))

    # -- pretty printing -----------------------------------------------------
    def symbolic_rows(self) -> List[List[str]]:
        """Fogler table with symbolic entries (Especie|Inicial|Cambio|
        Remanente|Concentracion), coefficients shown numerically."""
        rows = []
        for sp in self.species:
            th, nu_i = self.theta_of(sp), self.nu_of(sp)
            ini = "n_A0" if sp == self.limiting else f"{th:g}·n_A0"
            if nu_i == 0:
                cam, rem = "—", ini
            else:
                sign = "-" if nu_i < 0 else "+"
                cam = f"{sign}{abs(nu_i):g}·n_A0·X"
                rem = (f"n_A0(1 - X)" if sp == self.limiting
                       else f"n_A0({th:g} {sign} {abs(nu_i):g}X)")
            if self.variable_volume:
                conc = f"C_A0({th:g}{nu_i:+g}X)/(1{self.epsilon:+g}X)·(T0/T)(P/P0)"
            else:
                conc = f"C_A0({th:g}{nu_i:+g}X)"
            rows.append([sp, ini, cam, rem, conc])
        return rows

    def numeric_rows(self, X: float) -> List[List]:
        """Table evaluated at a given conversion X."""
        n, C = self.moles(X), self.concentrations(X)
        rows = []
        for sp in self.species:
            n0_i = self.n_A0 * self.theta_of(sp)
            rows.append([sp, n0_i, n[sp] - n0_i, n[sp], C[sp]])
        rows.append(["TOTAL", self.n_A0 * sum(self.theta.values()),
                     self.n_A0 * self.delta * X,
                     sum(n.values()), self.C_T(X)])
        return rows


# ----------------------------------------------------------------------------
# Automatic validation ("checks" section)
# ----------------------------------------------------------------------------

def run_checks_method_a(res: ExtentResult, n0: Dict[str, float],
                        atoms: Optional[AtomComposition] = None,
                        rel_tol: float = 1e-6) -> bool:
    """
    Print the standard checks for a Method-A result:
      1. atom balance in vs out (needs `atoms`)
      2. all n_i >= 0
      3. |K_calc - K_given|/K_given < rel_tol per reaction
      4. sum(y_i) = 1
    Returns True iff every check passes.
    """
    ok = True
    print("  Checks:")

    if atoms is not None:
        in_tot = atom_totals({sp: n0.get(sp, 0.0) for sp in res.species}, atoms)
        out_tot = atom_totals(res.moles, atoms)
        worst = max((abs(in_tot.get(el, 0) - out_tot.get(el, 0))
                     for el in set(in_tot) | set(out_tot)), default=0.0)
        good = worst < 1e-8
        ok &= good
        print(f"    [{'OK' if good else 'FAIL'}] atom balance in vs out "
              f"(max |Δ| = {worst:.2e} mol): "
              + ", ".join(f"{el}: {in_tot.get(el, 0):.6g} -> "
                          f"{out_tot.get(el, 0):.6g}"
                          for el in sorted(set(in_tot) | set(out_tot))))
    else:
        print("    [--] atom balance skipped (no atomic compositions given)")

    good = bool(np.all(res.n >= -1e-9))
    ok &= good
    print(f"    [{'OK' if good else 'FAIL'}] all n_i >= 0 "
          f"(min n_i = {res.n.min():.3e} mol)")

    rel = np.abs(res.K_calc - res.K_given) / res.K_given
    good = bool(np.all(rel < rel_tol))
    ok &= good
    print(f"    [{'OK' if good else 'FAIL'}] residuals "
          f"|K_calc-K|/K < {rel_tol:g} (max = {rel.max():.2e})")

    sy = float(res.y.sum())
    good = abs(sy - 1.0) < 1e-9
    ok &= good
    print(f"    [{'OK' if good else 'FAIL'}] sum(y_i) = {sy:.12f}")

    if not res.converged:
        ok = False
        print(f"    [FAIL] solver: {res.message}")
    if not ok:
        print("    *** WARNING: one or more checks FAILED ***")
    return ok


def run_checks_method_b(tbl: StoichiometricTable, X: float, Kc: Optional[float],
                        atoms: Optional[AtomComposition] = None,
                        rel_tol: float = 1e-6) -> bool:
    """
    Print the standard checks for a Method-B result at conversion X:
      1. atom balance feed vs outlet (needs `atoms`)
      2. all C_i (and n_i) >= 0
      3. |Kc_calc - Kc|/Kc < rel_tol (if Kc given)
      4. C_T constant = P/(RT) (gas, variable volume) — or reported P rise
         for the rigid-vessel case.
    Returns True iff every check passes.
    """
    ok = True
    print("  Checks:")

    n_out = tbl.moles(X)
    if atoms is not None:
        n_in = {sp: tbl.n_A0 * tbl.theta_of(sp) for sp in tbl.species}
        in_tot, out_tot = atom_totals(n_in, atoms), atom_totals(n_out, atoms)
        worst = max((abs(in_tot.get(el, 0) - out_tot.get(el, 0))
                     for el in set(in_tot) | set(out_tot)), default=0.0)
        good = worst < 1e-8
        ok &= good
        print(f"    [{'OK' if good else 'FAIL'}] atom balance in vs out "
              f"(max |Δ| = {worst:.2e} mol): "
              + ", ".join(f"{el}: {in_tot.get(el, 0):.6g} -> "
                          f"{out_tot.get(el, 0):.6g}"
                          for el in sorted(set(in_tot) | set(out_tot))))
    else:
        print("    [--] atom balance skipped (no atomic compositions given)")

    C = tbl.concentrations(X)
    min_c = min(C.values())
    good = min_c >= -1e-12 and min(n_out.values()) >= -1e-12
    ok &= good
    print(f"    [{'OK' if good else 'FAIL'}] all n_i, C_i >= 0 "
          f"(min C_i = {min_c:.3e} mol/L)")

    if Kc is not None:
        rel = abs(tbl.Kc_calc(X) - Kc) / Kc
        good = rel < rel_tol
        ok &= good
        print(f"    [{'OK' if good else 'FAIL'}] residual "
              f"|Kc_calc-Kc|/Kc = {rel:.2e} < {rel_tol:g}")

    if tbl.variable_volume:
        C_T_expected = total_concentration(tbl.P, tbl.T, tbl.pressure_unit)
        drift = max(abs(tbl.C_T(x) - C_T_expected)
                    for x in np.linspace(0.0, X, 11) if x <= tbl.X_max())
        good = drift < 1e-9
        ok &= good
        print(f"    [{'OK' if good else 'FAIL'}] C_T = P/(RT) = "
              f"{C_T_expected:.5g} mol/L constant over [0, X] "
              f"(max drift {drift:.2e})")
    else:
        P_eq = to_atm(tbl.P0, tbl.pressure_unit) * (1 + tbl.epsilon * X) \
            * tbl.T / tbl.T0
        print(f"    [OK] constant volume: C_T varies with X by design; "
              f"P rises to {P_eq:.4g} atm at X = {X:.4g}")

    if not ok:
        print("    *** WARNING: one or more checks FAILED ***")
    return ok

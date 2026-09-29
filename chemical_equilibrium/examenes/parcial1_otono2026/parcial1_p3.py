"""Parcial 1 — Problema 3 (4.5 pts)

Enunciado: en un reactor de flujo continuo ocurren simultáneamente
    R1: CH4 + H2O = CO + 3 H2      (balanceada; el examen la da esquelética)
    R2: CO + H2O = CO2 + H2        (WGS)
Alimentación 1 mol CH4 + 1 mol H2O; salida en equilibrio a 700 K y 1000 K.
Determinar avances y composiciones; discutir carácter endo/exotérmico.

Supuestos: gas ideal; P = 1 bar (el examen no la da; = P° de las K);
H2O como VAPOR (Hf298 = -241 818 J/mol, gas ideal, Tabla C.4); K(T) por
la ec. 13.18 de SVA con datos de data/sva_tables.json (Tablas C.1/C.4).
"""

import csv
import os
import pathlib
import sys

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))  # raíz del proyecto

import numpy as np

from core.equilibrium import atom_totals, format_table, solve_extents
from core.reactions import build_system, parse_reaction
from core.selfcheck import Problem, validate, validate_sweep
from core.thermo import cargar_datos, reaction_thermo

# ══════════════════════════ INPUTS ════════════════════════════════════════
reacciones = ["CH4 + H2O = CO + 3 H2", "CO + H2O = CO2 + H2"]
alim = {"CH4": 1.0, "H2O": 1.0}
P, P0, unidad = 1.0, 1.0, "bar"
temperaturas = [700.0, 1000.0]
especies_termo = ["CH4", "H2O", "CO", "H2", "CO2"]
XI0_NEWTON = (0.1, 0.01)               # valores iniciales que da el examen
# ══════════════════════════ FIN DE INPUTS ═════════════════════════════════

especies, nu, atomos, _ = build_system(reacciones, alim)
nu_dicts = [parse_reaction(rx) for rx in reacciones]
datos = cargar_datos(especies_termo)
nu_arr = np.array(nu, float)
delta = nu_arr.sum(axis=1)
n0_vec = np.array([alim.get(sp, 0.0) for sp in especies])

print("=" * 70)
print("  PARCIAL 1 — P3: reformado de metano, 1 CH4 + 1 H2O, 700/1000 K")
print("=" * 70)

print("\n1) ENUNCIADO Y SUPUESTOS")
print("   Gas ideal; P = 1 bar (supuesta, = P° de las K); H2O = VAPOR")
print("   (Tabla C.4, gas ideal). K(T) con la ec. 13.18 de SVA.")
for j, rx in enumerate(reacciones):
    print(f"   R{j + 1}: {rx}   (delta_{j + 1} = {delta[j]:+g})")

print("\n2) TABLA DE AVANCES (n_i = n_i0 + sum nu_ij·xi_j)")
expr = {"CH4": "1 - xi1", "H2O": "1 - xi1 - xi2", "CO": "xi1 - xi2",
        "H2": "3·xi1 + xi2", "CO2": "xi2"}
print("  " + format_table([[sp, alim.get(sp, 0.0), expr[sp]]
                           for sp in especies],
      ["Especie", "n_i0 (mol)", "n_i(xi1, xi2)"]).replace("\n", "\n  "))
print("   n_T = 2 + 2·xi1")

print("\n3) PARÁMETROS (Tablas C.4/C.1 vía ec. 13.18)")
t298 = [reaction_thermo(nud, datos, 298.15) for nud in nu_dicts]
for j, t in enumerate(t298):
    print(f"   R{j + 1}: dH°298 = {t['dH298']:+.0f} J/mol   "
          f"dG°298 = {t['dG298']:+.0f} J/mol   K298 = {t['K']:.3e}")

# resolver ambas temperaturas
resultados = {}
for T in temperaturas:
    termo = [reaction_thermo(nud, datos, T) for nud in nu_dicts]
    K = [t["K"] for t in termo]
    res = solve_extents(especies, alim, nu, K, P=P, P0=P0, atoms=atomos)
    resultados[T] = {"termo": termo, "K": K, "res": res}

print("\n4) ECUACIONES A RESOLVER (números sustituidos; P = P° -> factor 1)")
for T in temperaturas:
    K1, K2 = resultados[T]["K"]
    print(f"   T = {T:g} K:")
    print(f"     y_CO·y_H2³/(y_CH4·y_H2O)·(P/P°)² = K1 = {K1:.4e}")
    print(f"     y_H2·y_CO2/(y_H2O·y_CO)          = K2 = {K2:.4f}")

print("\n5) RESULTADOS")
for T in temperaturas:
    res = resultados[T]["res"]
    X = (alim["CH4"] - res.moles["CH4"]) / alim["CH4"]
    resultados[T]["X"] = X
    print(f"   T = {T:g} K:  xi1 = {res.xi[0]:.5f}  xi2 = {res.xi[1]:.5f}  "
          f"n_T = {res.n_total:.4f} mol   X_CH4 = {X:.4f}")
    print("  " + format_table(
        [[sp, res.moles[sp], res.mole_fractions[sp]] for sp in especies],
        ["Especie", "n_i (mol)", "y_i"]).replace("\n", "\n  "))

print("\n6) COMPROBACIÓN")
for T in temperaturas:
    res = resultados[T]["res"]
    rel = np.abs(res.K_calc - res.K_given) / res.K_given
    ent = atom_totals(alim, atomos)
    sal = atom_totals(res.moles, atomos)
    peor = max(abs(ent.get(e, 0.0) - sal.get(e, 0.0))
               for e in set(ent) | set(sal))
    print(f"   T = {T:g} K: |K_calc-K|/K máx = {rel.max():.1e};  átomos "
          f"máx |Δ| = {peor:.1e} mol;  sum y_i = {res.y.sum():.6f}")

print("\n7) SELFCHECK (20 arranques aleatorios por T + Le Chatelier)")
for T in temperaturas:
    prob = Problem(species=especies, n0=alim, nu=nu, K=resultados[T]["K"],
                   P=P, P0=P0, atoms=atomos, T=T, pressure_unit=unidad,
                   label="P3 parcial")
    validate(prob, resultados[T]["res"], n_starts=20, k_tol=1e-4)
validate_sweep({T: resultados[T]["K"] for T in temperaturas},
               {T: tuple(resultados[T]["res"].xi) for T in temperaturas})

print("\n8) TABLA TERMODINÁMICA POR REACCIÓN")
filas = []
for T in temperaturas:
    for j, t in enumerate(resultados[T]["termo"]):
        filas.append([f"{T:g}", f"R{j + 1}", t["dCp_R"], t["dH"] / 1000.0,
                      t["dG"] / 1000.0, t["K"]])
print("  " + format_table(filas, ["T (K)", "Rxn", "dCp/R", "dH° (kJ/mol)",
                                  "dG° (kJ/mol)", "K"]).replace("\n", "\n  "))

print("\n9) TRAYECTORIA DE NEWTON desde (xi1, xi2) = (0.1, 0.01)")
print("   (forma logarítmica de las MISMAS ecuaciones: r_j = ln(Pi y^nu_j)"
      " - ln K_j; Jacobiano por diferencias finitas; paso amortiguado para"
      " mantener n_i > 0)")


def _residuo(xi, lnK):
    n = n0_vec + nu_arr.T @ xi
    if np.min(n) <= 0:
        return None
    y = n / n.sum()
    return nu_arr @ np.log(y) - lnK


def newton(xi0, lnK, tol=1e-10, itmax=25):
    xi = np.array(xi0, float)
    r = _residuo(xi, lnK)
    historia = [(0, xi.copy(), float(np.max(np.abs(r))), 1.0)]
    for it in range(1, itmax + 1):
        J = np.zeros((2, 2))
        for k in range(2):                    # Jacobiano numérico
            h = 1e-8 * max(1.0, abs(xi[k]))
            xp = xi.copy()
            xp[k] += h
            J[:, k] = (_residuo(xp, lnK) - r) / h
        paso = np.linalg.solve(J, -r)
        lam = 1.0
        while lam > 1e-6:                     # amortiguamiento: n_i > 0
            r_new = _residuo(xi + lam * paso, lnK)
            if r_new is not None and np.max(np.abs(r_new)) < np.max(np.abs(r)):
                break
            lam *= 0.5
        xi = xi + lam * paso
        r = _residuo(xi, lnK)
        historia.append((it, xi.copy(), float(np.max(np.abs(r))), lam))
        if np.max(np.abs(r)) < tol:
            break
    return xi, historia


for T in temperaturas:
    lnK = np.log(np.array(resultados[T]["K"]))
    xi_n, hist = newton(XI0_NEWTON, lnK)
    print(f"\n   T = {T:g} K:")
    print(f"     {'it':>3s} {'xi1':>10s} {'xi2':>10s} {'max|r|':>10s} "
          f"{'paso lam':>9s}")
    for it, xi, err, lam in hist:
        print(f"     {it:3d} {xi[0]:10.6f} {xi[1]:10.6f} {err:10.2e} "
              f"{lam:9.4f}")
    dif = np.max(np.abs(xi_n - resultados[T]["res"].xi))
    print(f"     -> coincide con el solver multi-arranque "
          f"(dif máx = {dif:.1e})")

# ---------------- CSV y gráfica ----------------
ruta_csv = os.path.join(_DIR, "resultados_P3.csv")
with open(ruta_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["T_K", "xi1", "xi2", "X_CH4", "n_T"]
               + [f"n_{sp}" for sp in especies]
               + [f"y_{sp}" for sp in especies]
               + ["K1", "K2", "dH_R1_J_mol", "dH_R2_J_mol",
                  "dG_R1_J_mol", "dG_R2_J_mol"])
    for T in temperaturas:
        res, tm = resultados[T]["res"], resultados[T]["termo"]
        w.writerow([T, f"{res.xi[0]:.5f}", f"{res.xi[1]:.5f}",
                    f"{resultados[T]['X']:.5f}", f"{res.n_total:.4f}"]
                   + [f"{res.moles[sp]:.5f}" for sp in especies]
                   + [f"{res.mole_fractions[sp]:.5f}" for sp in especies]
                   + [f"{tm[0]['K']:.5g}", f"{tm[1]['K']:.5g}",
                      f"{tm[0]['dH']:.0f}", f"{tm[1]['dH']:.0f}",
                      f"{tm[0]['dG']:.0f}", f"{tm[1]['dG']:.0f}"])
print(f"\n   CSV: {ruta_csv}")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.arange(len(especies))
    ancho, TINTA = 0.38, "#1a1a19"
    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    for k, (T, color) in enumerate(zip(temperaturas,
                                       ("#2a78d6", "#eb6834"))):
        ys = [resultados[T]["res"].mole_fractions[sp] for sp in especies]
        barras = ax.bar(x + (k - 0.5) * ancho, ys, ancho * 0.94,
                        color=color, label=f"{T:g} K")
        for b, yv in zip(barras, ys):
            ax.annotate(f"{yv:.3f}", (b.get_x() + b.get_width() / 2, yv),
                        ha="center", va="bottom", fontsize=7.5, color=TINTA)
    ax.set_xticks(x, especies)
    ax.set_ylabel("y_i (fracción mol)", color=TINTA)
    ax.set_title("P3: composición de equilibrio a 700 y 1000 K (1 bar)",
                 color=TINTA, fontsize=10.5)
    ax.grid(axis="y", color="#dddddd", lw=0.7)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color("#bbbbbb")
    ax.tick_params(colors="#555555", labelsize=8.5)
    ax.legend(frameon=False, fontsize=9, labelcolor=TINTA)
    fig.tight_layout()
    ruta_png = os.path.join(_DIR, "grafica_P3.png")
    fig.savefig(ruta_png, facecolor="white")
    plt.close(fig)
    print(f"   Gráfica: {ruta_png}")
except ImportError:
    print("   (sin gráfica: instala matplotlib)")

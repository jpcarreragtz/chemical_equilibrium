"""Ejercicios adicionales 2do parcial — PROBLEMA 1 (síntesis de amoniaco)

    1/2 N2(g) + 3/2 H2(g) -> NH3(g)        alimentación estequiométrica

a) y_NH3 de equilibrio a 1 bar y 300 K                     (gas ideal)
b) T a la que y_NH3 = 0.50 a 1 bar                          (gas ideal)
c) T a la que y_NH3 = 0.50 a 100 bar, mezcla = gas ideal
d) T a la que y_NH3 = 0.50 a 100 bar, mezcla = SOLUCIÓN IDEAL de gases
   (regla de Lewis-Randall: phi^_i = phi_i de especie pura a T y P;
    phi por la virial generalizada de Pitzer, Sesión 13)

Planteamiento (a mano)
----------------------
Base: n_N2,0 = 1/2, n_H2,0 = 3/2 (proporción estequiométrica 1:3).
ν = {NH3: +1, N2: −1/2, H2: −3/2}  →  δ = Σν = −1
    n_NH3 = ξ     n_N2 = ½(1 − ξ)     n_H2 = 3/2 (1 − ξ)     n_T = 2 − ξ
    y_NH3 = ξ/(2 − ξ)    y_N2 = ½(1 − ξ)/(2 − ξ)    y_H2 = 3/2(1 − ξ)/(2 − ξ)

    K = K_phi · K_y · (P/P°)^δ ,   K_y = y_NH3 / (y_N2^½ · y_H2^3/2)
    gas ideal:  K_phi = 1                →  K_y = K·(P/P°)
    sol. ideal: K_phi = φ_NH3/(φ_N2^½ φ_H2^3/2)  →  K_y = K·(P/P°)/K_phi

    y_NH3 = 0.5  →  ξ = 2/3, y_N2 = 0.125, y_H2 = 0.375  →  K_y = 6.1584
    (la composición NO depende de T ni de P: sólo cambia qué K se necesita)

K(T) por la ec. 13.18 de SVA (Tablas C.1/C.4, data/sva_tables.json), como
en el 1er parcial. b), c), d) se resuelven con brentq sobre T:
    ln K(T) − ln[K_y · (P/P°)^(−δ)... ]  = 0   (ver funciones abajo)
"""

import sys, pathlib, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz

from math import exp, log
from scipy.optimize import brentq

from core.equilibrium import format_table
from core.nonideal import K_phi_lewis, cargar_criticas
from core.thermo import R_J, T_REF, cargar_datos, delta_property, reaction_thermo

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
NU = {"NH3": 1.0, "N2": -0.5, "H2": -1.5}
N0 = {"N2": 0.5, "H2": 1.5, "NH3": 0.0}     # alimentación estequiométrica
P0 = 1.0                                    # bar
Y_OBJ = 0.50
CLAVE = {"a": 0.9349, "b": 399.5, "c": 577.6, "d": 568.6}   # maestra
COMPANERO = {"c": 551.0, "d": 557.0}
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

DELTA = sum(NU.values())
datos = cargar_datos(list(NU))
crit = cargar_criticas(list(NU))


def comp(xi):
    """n_i, y_i en función de ξ (tabla del planteamiento)."""
    n = {sp: N0[sp] + NU[sp] * xi for sp in NU}
    nT = sum(n.values())
    return n, {sp: v / nT for sp, v in n.items()}, nT


def Ky(xi):
    _, y, _ = comp(xi)
    return y["NH3"] / (y["N2"] ** 0.5 * y["H2"] ** 1.5)


def lnK(T):
    return log(reaction_thermo(NU, datos, T)["K"])


def resolver_xi(K, P, Kphi=1.0):
    """ξ tal que K_y(ξ) = K·(P/P°)^(−δ)/K_phi."""
    obj = log(K) - DELTA * log(P / P0) - log(Kphi)
    return brentq(lambda x: log(Ky(x)) - obj, 1e-12, 1 - 1e-12, xtol=1e-15)


def T_para_y(P, modelo="gi", T_lo=300.0, T_hi=1500.0, invertir_phi=False,
             criticas=None):
    """T tal que y_NH3 = Y_OBJ. modelo 'gi' (K_phi = 1) o 'lewis'."""
    xi = 2 * Y_OBJ / (1 + Y_OBJ)                       # y = ξ/(2−ξ)
    ky = Ky(xi)

    def f(T):
        Kphi = 1.0
        if modelo == "lewis":
            Kphi = K_phi_lewis(NU, T, P, criticas or crit, avisar=False)[0]
            if invertir_phi:
                Kphi = 1.0 / Kphi
        # K requerida = K_phi · K_y · (P/P°)^δ
        return lnK(T) - (log(Kphi) + log(ky) + DELTA * log(P / P0))

    return brentq(f, T_lo, T_hi, xtol=1e-10), xi, ky


print("=" * 72)
print("  PROBLEMA 1 — 1/2 N2 + 3/2 H2 -> NH3  (alimentación estequiométrica)")
print("=" * 72)

t298 = reaction_thermo(NU, datos, T_REF)
dA, dB, dC, dD = (delta_property(NU, datos, k) for k in "ABCD")
print("\n0) DATOS TERMO (Tablas C.1/C.4 SVA)")
print(f"   ΔH°298 = {t298['dH298']:.0f} J/mol   ΔG°298 = {t298['dG298']:.0f} "
      f"J/mol   K298 = {t298['K']:.2f}")
print(f"   ΔA = {dA:.4f}  ΔB = {dB:.4e}  ΔC = {dC:.1e}  ΔD = {dD:.4e}")
print(f"   δ = Σν = {DELTA:+.0f}  → K = K_phi·K_y·(P/P°)^δ")

# ---------------------------------------------------------------- a)
Ta, Pa = 300.0, 1.0
Ka = exp(lnK(Ta))
xia = resolver_xi(Ka, Pa)
na, ya, nTa = comp(xia)
print("\na) T = 300 K, P = 1 bar, gas ideal")
print(f"   K(300) = {Ka:.3f}   (ln K = {log(Ka):.4f})")
print(f"   ξ/(2−ξ) ... resuelve K_y(ξ) = K → ξ = {xia:.6f},  n_T = {nTa:.6f}")
print("  " + format_table([[sp, na[sp], ya[sp]] for sp in NU],
                          ["especie", "n_i", "y_i"]).replace("\n", "\n  "))
print(f"   → y_NH3 = {ya['NH3']:.4f}     (maestra {CLAVE['a']})")

# ---------------------------------------------------------------- b) c)
Tb, xi_obj, ky = T_para_y(1.0)
Tc, _, _ = T_para_y(100.0)
print(f"\n   Para y_NH3 = 0.50: ξ = {xi_obj:.4f}, y_N2 = 0.125, y_H2 = 0.375, "
      f"K_y = {ky:.4f}")
print(f"\nb) P = 1 bar, gas ideal: K requerida = K_y·(1)^(-1) = {ky:.4f}")
print(f"   → T = {Tb:.2f} K     (maestra {CLAVE['b']} K)")
print(f"\nc) P = 100 bar, gas ideal: K requerida = K_y/100 = {ky / 100:.5f}")
print(f"   → T = {Tc:.2f} K     (maestra {CLAVE['c']} K; compañero "
      f"{COMPANERO['c']} K)")

# ---------------------------------------------------------------- d)
Td, _, _ = T_para_y(100.0, "lewis")
Kphi_d, det = K_phi_lewis(NU, Td, 100.0, crit, avisar=False)
print("\nd) P = 100 bar, SOLUCIÓN IDEAL (Lewis-Randall, φ_i puro, Pitzer)")
filas = [[sp, *crit[sp], det[sp]["Tr"], det[sp]["Pr"], det[sp]["B0"],
          det[sp]["B1"], det[sp]["phi"], "sí" if det[sp]["valido"] else "NO*"]
         for sp in NU]
print("  " + format_table(filas, ["esp", "Tc", "Pc", "ω", "Tr", "Pr", "B0",
                                  "B1", "φ", "Vr≥2"]).replace("\n", "\n  "))
print(f"   K_phi = φ_NH3/(φ_N2^½ φ_H2^3/2) = {Kphi_d:.5f}")
print(f"   K requerida = K_phi·K_y/100 = {Kphi_d * ky / 100:.5f}"
      f"  = K({Td:.2f} K) = {exp(lnK(Td)):.5f}")
print(f"   → T = {Td:.2f} K   (maestra {CLAVE['d']} K; compañero "
      f"{COMPANERO['d']} K)")

Td_inv, _, _ = T_para_y(100.0, "lewis", invertir_phi=True)
crit_sin_w = {sp: (Tc_, Pc_, 0.0) for sp, (Tc_, Pc_, _w) in crit.items()}
Td_clave, _, _ = T_para_y(100.0, "lewis", invertir_phi=True,
                          criticas=crit_sin_w)
# Lee-Kesler (Tablas E.13/E.14) para el NH3 como contraste del φ de Pitzer
TrN, PrN = Td / crit["NH3"][0], 100.0 / crit["NH3"][1]
def _bil(v00, v01, v10, v11, tr, pr):          # Tr 1.40/1.50, Pr 0.8/1.0
    a = (tr - 1.40) / 0.10; b = (pr - 0.8) / 0.2
    return (v00 * (1 - a) * (1 - b) + v01 * (1 - a) * b + v10 * a * (1 - b)
            + v11 * a * b)
phi0 = _bil(0.9120, 0.8892, 0.9290, 0.9141, TrN, PrN)
phi1 = _bil(1.0642, 1.0815, 1.0666, 1.0865, TrN, PrN)
phiLK = phi0 * phi1 ** crit["NH3"][2]

print("\n   NOTA sobre d)  (DISCREPANCIA CON LA CLAVE)")
print(f"   * φ_NH3 < 1 y φ_N2, φ_H2 > 1  → K_phi < 1 → se necesita una K MENOR")
print(f"     que en c) → T_d > T_c. Físicamente: la no idealidad a 100 bar")
print(f"     FAVORECE al NH3 (por eso T_d = {Td:.1f} K > T_c = {Tc:.1f} K).")
print(f"   * Si K_phi se aplica INVERTIDO (K_y = K·P·K_phi en vez de ÷K_phi)")
print(f"     sale T = {Td_inv:.1f} K ≈ clave de la maestra ({CLAVE['d']} K).")
print(f"   * Invertido Y sin el término ω·B1 (ln φ = Pr/Tr·B0) sale "
      f"T = {Td_clave:.1f} K = clave exacta.")
print(f"     → la clave combina 2 cosas: K_phi del lado equivocado (~14 K) "
      f"y omitir ω·B1 (~2 K).")
print(f"   * φ_NH3 por Lee-Kesler (Tablas E.13/E.14) a Tr = {TrN:.3f}, "
      f"Pr = {PrN:.3f}: {phiLK:.4f}")
print(f"     vs Pitzer {det['NH3']['phi']:.4f} → el método de φ NO explica la "
      f"diferencia.")
print("   * Vr < 2 para NH3 y N2: la virial está en el borde de validez")
print("     (warning de SVA Fig. 3.14); es lo que se ve en clase (Sesión 13).")

# ---------------------------------------------------------------- checks
print("\nAUTOVALIDACIÓN")
ok = True
for nombre, T, P, Kphi in (("b", Tb, 1.0, 1.0), ("c", Tc, 100.0, 1.0),
                           ("d", Td, 100.0, Kphi_d)):
    xi = resolver_xi(exp(lnK(T)), P, Kphi)
    y = comp(xi)[1]["NH3"]
    bien = abs(y - Y_OBJ) < 1e-8
    ok &= bien
    print(f"  [{'PASS' if bien else 'FAIL'}] {nombre}) re-resolviendo a "
          f"T = {T:.2f} K: y_NH3 = {y:.8f}")
n_at = {"N": 2 * na["N2"] + na["NH3"], "H": 2 * na["H2"] + 3 * na["NH3"]}
bien = abs(n_at["N"] - 1.0) < 1e-12 and abs(n_at["H"] - 3.0) < 1e-12
ok &= bien
print(f"  [{'PASS' if bien else 'FAIL'}] balance atómico a): N = "
      f"{n_at['N']:.12f} (1), H = {n_at['H']:.12f} (3)")
for k, v, cl in (("a", ya["NH3"], CLAVE["a"]), ("b", Tb, CLAVE["b"]),
                 ("c", Tc, CLAVE["c"])):
    bien = abs(v - cl) / cl < 5e-4
    ok &= bien
    print(f"  [{'PASS' if bien else 'FAIL'}] {k}) vs clave: {v:.4f} vs {cl}")
bien = Td > Tc
ok &= bien
print(f"  [{'PASS' if bien else 'FAIL'}] tendencia física d) > c) "
      f"({Td:.1f} > {Tc:.1f} K)")
print(f"  [info] d) clave {CLAVE['d']} K ≠ {Td:.1f} K: ver NOTA (K_phi invertido)")
print(f"\nVEREDICTO: {'OK' if ok else 'REVISAR'}")

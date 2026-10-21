"""SIMULACRO a) — Presión + fugacidad: síntesis de metanol (Sesión 13)

    CO(g) + 2 H2(g) -> CH3OH(g)      alimentación estequiométrica (1 : 2)
    Incisos (elige con el argumento:  python3 a_metanol_presion.py 1|2|3):
      1) y_CH3OH a 300 K y 1 bar (gas ideal)
      2) T a la que y_CH3OH = 0.5 a 1 bar (gas ideal)
      3) T a la que y_CH3OH = 0.5 a 100 bar, Lewis-Randall (φ puro, Pitzer),
         convención de K_phi por default (core/convenciones.py → "maestra")
    En los JSON el metanol se llama CH4O (data/sva_tables.json y
    data/criticas_b1.json: Tc = 512.6 K, Pc = 80.97 bar, ω = 0.564).

PLANTEAMIENTO A MANO
--------------------
ν = {CO: −1, H2: −2, CH4O: +1}   →   δ = Σν = −2
Base: n_CO,0 = 1, n_H2,0 = 2
    n_CO = 1 − ξ      n_H2 = 2 − 2ξ      n_CH4O = ξ      n_T = 3 − 2ξ
    y_CO = (1 − ξ)/(3 − 2ξ)   y_H2 = 2(1 − ξ)/(3 − 2ξ)   y_CH4O = ξ/(3 − 2ξ)
    K_y = y_CH4O/(y_CO·y_H2²) = ξ(3 − 2ξ)² / [4(1 − ξ)³]
K(T) por la ec. 13.18 (ΔH°298 = −90 135 J/mol, ΔG°298 = −24 791 J/mol).
    gas ideal:  K = K_y·(P/P°)^δ = K_y/P²      →  K_y = K·P²
    Lewis-Randall, convención "maestra":  K = K_y·(P/P°)^δ / K_phi
                                          →  K_y = K·P²·K_phi
    K_phi = φ_CH4O/(φ_CO·φ_H2²), φ_i de especie pura a (T, P), Pitzer:
            ln φ = (Pr/Tr)(B⁰ + ωB¹), B⁰ = 0.083 − 0.422/Tr^1.6, B¹ = 0.139 − 0.172/Tr^4.2
1) T = 300 K, P = 1 bar: resolver ξ de K_y(ξ) = K(300).
2) y = 0.5 → ξ/(3 − 2ξ) = 0.5 → ξ = 0.75 → y_CO = 1/6, y_H2 = 1/3
   → K_y = 0.5/[(1/6)(1/3)²] = 27 (fijo).   Resolver K(T) = 27/P² con P = 1.
3) Resolver K(T) = 27/(P²·K_phi(T)) con P = 100 (K_phi depende de T: iterar).
   δ < 0 → la presión favorece al metanol → T(100 bar) > T(1 bar).
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz (simulacro/ está 3 niveles abajo)

import contextlib, io, os
from math import exp, log
from scipy.optimize import brentq

from core.convenciones import (CONVENCION_KPHI, describir, factor_Kphi,
                               K_efectiva, otra_convencion, validar_convencion)
from core.equilibrium import format_table
from core.exercise import run_exercise
from core.nonideal import K_phi_lewis, cargar_criticas
from core.reactions import parse_formula, parse_reaction
from core.thermo import T_REF, cargar_datos, delta_property, reaction_thermo

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
INCISO = sys.argv[1] if len(sys.argv) > 1 else "3"
CASOS = {
    "1": dict(titulo="Simulacro a1) metanol: y_CH3OH a 300 K y 1 bar (gas ideal)",
              P=1.0, modo="y_a_T", modelo="gi"),
    "2": dict(titulo="Simulacro a2) metanol: T para y_CH3OH = 0.5 a 1 bar (gas ideal)",
              P=1.0, modo="T_para_y", modelo="gi"),
    "3": dict(titulo="Simulacro a3) metanol: T para y_CH3OH = 0.5 a 100 bar (Lewis-Randall)",
              P=100.0, modo="T_para_y", modelo="lewis"),
}
_c = CASOS[INCISO]
titulo = _c["titulo"]
reaccion = "CO + 2 H2 = CH4O"              # CH4O = CH3OH (nombre en los JSON)
alimentacion = {"CO": 1.0, "H2": 2.0}      # estequiométrica
P, P0 = _c["P"], 1.0                       # bar

modo = _c["modo"]                          # "y_a_T" | "T_para_y"
modelo = _c["modelo"]                      # "gi" | "lewis"
convencion_Kphi = None                     # None = core.convenciones (default "maestra")

T = 300.0                                  # modo "y_a_T": temperatura [K]
especie_obj, y_obj = "CH4O", 0.50          # modo "T_para_y": objetivo
T_lo, T_hi = 300.0, 1500.0                 # intervalo de búsqueda de T [K]

especies_termo = None     # None → las de la reacción (data/sva_tables.json)
criticas = None           # None → data/criticas_b1.json
clave = None              # sin clave de la maestra: se verifica abajo a mano
tol_clave = 0.0
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

# ─── A partir de aquí: COPIA LITERAL de plantilla/plantilla_presion_fugacidad.py (no editar) ───

# ---------------------------------------------------------------- sistema
nu = parse_reaction(reaccion)
especies = list(nu) + [sp for sp in alimentacion if sp not in nu]
inertes = [sp for sp in especies if sp not in nu]
n0 = {sp: float(alimentacion.get(sp, 0.0)) for sp in especies}
delta = sum(nu.values())
datos = cargar_datos(especies_termo or list(nu))
con_phi = modelo == "lewis"
if modelo not in ("gi", "lewis"):
    raise ValueError("modelo debe ser 'gi' o 'lewis'")
if con_phi:
    conv = validar_convencion(convencion_Kphi)
    crit = criticas or cargar_criticas(list(nu))
else:
    conv, crit = None, None

# límites físicos de ξ: todas las n_i ≥ 0
xi_hi = min(n0[sp] / -v for sp, v in nu.items() if v < 0)
xi_lo = max([-n0[sp] / v for sp, v in nu.items() if v > 0] + [-xi_hi * 1e3]) + 0.0
EPS = 1e-12 * (xi_hi - xi_lo)


def comp(xi):
    n = {sp: n0[sp] + nu.get(sp, 0.0) * xi for sp in especies}
    nT = sum(n.values())
    return n, {sp: v / nT for sp, v in n.items()}, nT


def lnKy(xi):
    _, y, _ = comp(xi)
    return sum(v * log(y[sp]) for sp, v in nu.items())


def lnK(T_):
    return log(reaction_thermo(nu, datos, T_)["K"])


def Kphi(T_):
    return K_phi_lewis(nu, T_, P, crit, avisar=False) if con_phi else (1.0, {})


def F_de(T_):
    return factor_Kphi(Kphi(T_)[0], conv) if con_phi else 1.0


def resolver_xi(T_):
    """ξ tal que ln K_y(ξ) = ln K(T) − δ·ln(P/P°) − ln F(T)."""
    obj = lnK(T_) - delta * log(P / P0) - log(F_de(T_))
    return brentq(lambda x: lnKy(x) - obj, xi_lo + EPS, xi_hi - EPS, xtol=1e-15)


def xi_para_y(sp, y):
    """ξ tal que y_sp(ξ) = y."""
    return brentq(lambda x: comp(x)[1][sp] - y, xi_lo + EPS, xi_hi - EPS,
                  xtol=1e-15)


def fmt_n(sp):
    v = nu.get(sp, 0.0)
    if v == 0:
        return f"{n0[sp]:g}"
    return f"{n0[sp]:g} {'+' if v > 0 else '−'} {abs(v):g}·ξ"


# ---------------------------------------------------------------- impresión
print("=" * 72 + f"\n  {titulo}\n" + "=" * 72)
print(f"\n1) REACCIÓN Y ν:   {reaccion}")
print("   ν = {" + ", ".join(f"{sp}: {v:+g}" for sp, v in nu.items())
      + f"}}   →   δ = Σν = {delta:+g}")
print(f"   Inertes: {', '.join(inertes) if inertes else '(ninguno)'}")
print("\n2) TABLA EN ξ (n_i = n_i0 + ν_i·ξ):")
nT_expr = f"{sum(n0.values()):g}" + (f" {'+' if delta > 0 else '−'} {abs(delta):g}·ξ"
                                   if abs(delta) > 1e-12 else "")
print("  " + format_table([[sp, n0[sp], nu.get(sp, 0.0), fmt_n(sp),
                            f"({fmt_n(sp)})/({nT_expr})"] for sp in especies],
                          ["especie", "n_i0", "ν_i", "n_i(ξ)", "y_i(ξ)"]
                          ).replace("\n", "\n  "))
print(f"   n_T = {nT_expr};   ξ ∈ ({xi_lo:.4g}, {xi_hi:.4g})")
print("   K_y = Π y_i^ν_i = " + " · ".join(f"y_{sp}^({v:g})" for sp, v in nu.items()))
print(f"   K = F·K_y·(P/P°)^δ,  P/P° = {P / P0:g}")
if con_phi:
    print(f"   Modelo: solución ideal (Lewis-Randall, φ_i puro por Pitzer)")
    print(f"   Convención '{conv}':  {describir(conv)}")
else:
    print("   Modelo: gas ideal (K_phi = 1, F = 1)")

t298 = reaction_thermo(nu, datos, T_REF)
dA, dB, dC, dD = (delta_property(nu, datos, k) for k in "ABCD")
print("\n3) K(T) POR LA EC. 13.18 (Tablas C.1/C.4 SVA)")
print(f"   ΔH°298 = {t298['dH298']:.0f} J/mol   ΔG°298 = {t298['dG298']:.0f} "
      f"J/mol   K298 = {t298['K']:.4g}")
print(f"   ΔA = {dA:.4f}  ΔB = {dB:.4e}  ΔC = {dC:.4e}  ΔD = {dD:.4e}")

# ---------------------------------------------------------------- resolver
if modo == "y_a_T":
    T_sol = float(T)
    xi_sol = resolver_xi(T_sol)
elif modo == "T_para_y":
    xi_sol = xi_para_y(especie_obj, y_obj)
    ky_obj = lnKy(xi_sol)
    T_sol = brentq(lambda T_: lnK(T_) - (log(F_de(T_)) + ky_obj
                                        + delta * log(P / P0)),
                   T_lo, T_hi, xtol=1e-10)
else:
    raise ValueError("modo debe ser 'y_a_T' o 'T_para_y'")

n_sol, y_sol, nT_sol = comp(xi_sol)
K_T = exp(lnK(T_sol))
Kphi_sol, det = Kphi(T_sol)
F_sol = F_de(T_sol)
K_ef = K_efectiva(K_T, Kphi_sol, conv) if con_phi else K_T
termo = reaction_thermo(nu, datos, T_sol)

print(f"\n4) SOLUCIÓN  (modo '{modo}', modelo '{modelo}')")
print(f"   T = {T_sol:.2f} K   P = {P:g} bar   K(T) = {K_T:.5g}   "
      f"ΔH°(T) = {termo['dH'] / 1000:.2f} kJ/mol   ΔG°(T) = {termo['dG'] / 1000:.2f} kJ/mol")
if con_phi:
    filas = [[sp, *crit[sp], det[sp]["Tr"], det[sp]["Pr"], det[sp]["B0"],
              det[sp]["B1"], det[sp]["phi"], "sí" if det[sp]["valido"] else "NO*"]
             for sp in nu]
    print("  " + format_table(filas, ["esp", "Tc", "Pc", "ω", "Tr", "Pr", "B0",
                                      "B1", "φ", "Vr≥2"]).replace("\n", "\n  "))
    if any(not det[sp]["valido"] for sp in nu):
        print("   * Vr < 2: φ extrapolado fuera de la zona virial (Fig. 3.14 SVA)")
    print(f"   K_phi = Π φ_i^ν_i = {Kphi_sol:.5f}   →   F = {F_sol:.5f}   "
          f"(convención '{conv}')")
    print(f"   K_ef = K/F = {K_ef:.5g}   (lo que resuelve el solver ideal: "
          f"K_y·(P/P°)^δ = K_ef)")
print(f"   K_y requerida = K·(P/P°)^(−δ)/F = {exp(lnKy(xi_sol)):.5g}   →   ξ = {xi_sol:.6f}")
print("  " + format_table([[sp, n_sol[sp], y_sol[sp]] for sp in especies],
                          ["especie", "n_i", "y_i"]).replace("\n", "\n  "))
print(f"   n_T = {nT_sol:.6f}")
if modo == "T_para_y":
    print(f"\n   → T = {T_sol:.2f} K  para  y_{especie_obj} = {y_obj}"
          + (f"     (clave {clave})" if clave is not None else ""))
    res_principal = T_sol
else:
    print(f"\n   → y_i a {T_sol:.2f} K: " + ", ".join(f"{sp} {y_sol[sp]:.4f}"
                                                     for sp in especies)
          + (f"     (clave {clave})" if clave is not None else ""))
    res_principal = max(y_sol.values()) if clave is None else y_sol.get(
        especie_obj, max(y_sol.values()))

# comparación de modelos / convenciones (informativa)
if con_phi:
    otra = otra_convencion(conv)
    if modo == "T_para_y":
        T_gi = brentq(lambda T_: lnK(T_) - (ky_obj + delta * log(P / P0)),
                      T_lo, T_hi, xtol=1e-10)
        T_otra = brentq(lambda T_: lnK(T_) - (
            log(factor_Kphi(Kphi(T_)[0], otra)) + ky_obj + delta * log(P / P0)),
            T_lo, T_hi, xtol=1e-10)
        print(f"   Comparación: gas ideal T = {T_gi:.2f} K;  '{otra}' T = {T_otra:.2f} K")
    else:
        obj_gi = lnK(T_sol) - delta * log(P / P0)
        y_gi = comp(brentq(lambda x: lnKy(x) - obj_gi, xi_lo + EPS, xi_hi - EPS))[1]
        obj_ot = obj_gi - log(factor_Kphi(Kphi_sol, otra))
        y_ot = comp(brentq(lambda x: lnKy(x) - obj_ot, xi_lo + EPS, xi_hi - EPS))[1]
        sp = especie_obj if especie_obj in y_sol else especies[0]
        print(f"   Comparación y_{sp}: gas ideal {y_gi[sp]:.4f};  '{otra}' {y_ot[sp]:.4f}")

# ---------------------------------------------------------------- checks
print("\nAUTOVALIDACIÓN")
ok = True
def chk(cond, txt):
    global ok
    ok &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {txt}")

res_ln = lnKy(xi_sol) - (lnK(T_sol) - delta * log(P / P0) - log(F_sol))
chk(abs(res_ln) < 1e-8, f"residual ln K_y − [ln K − δ ln(P/P°) − ln F] = {res_ln:.1e}")
chk(abs(sum(y_sol.values()) - 1) < 1e-12, f"Σ y_i = {sum(y_sol.values()):.12f}")
chk(all(v >= -1e-12 for v in n_sol.values()), "todas las n_i ≥ 0")
try:
    at_in = {}; at_out = {}
    for sp in especies:
        for el, k in parse_formula(sp).items():
            at_in[el] = at_in.get(el, 0.0) + n0[sp] * k
            at_out[el] = at_out.get(el, 0.0) + n_sol[sp] * k
    err_at = max(abs(at_out[el] - at_in[el]) for el in at_in)
    chk(err_at < 1e-10, f"balance atómico cierra (máx |Δ| = {err_at:.1e})")
except Exception as e:                      # especies abstractas (A, B...)
    print(f"  [info] balance atómico no evaluado ({e})")
if modo == "T_para_y":
    chk(abs(y_sol[especie_obj] - y_obj) < 1e-10,
        f"y_{especie_obj}(T = {T_sol:.2f} K) = {y_sol[especie_obj]:.10f} = objetivo")
if con_phi:
    chk(abs(K_ef * F_sol - K_T) / K_T < 1e-12,
        f"K_ef·F = K  (convención '{conv}' aplicada una sola vez)")
    tend = ("K_phi < 1 → " + ("F > 1" if F_sol > 1 else "F < 1")) if Kphi_sol < 1 \
        else ("K_phi > 1 → " + ("F > 1" if F_sol > 1 else "F < 1"))
    print(f"  [info] K_phi = {Kphi_sol:.4f}, F = {F_sol:.4f} ({tend}); "
          f"convención '{conv}'")
# contraste con el motor general (run_exercise) a la T de la solución
with contextlib.redirect_stdout(io.StringIO()):
    rx = run_exercise(titulo="contraste plantilla_presion_fugacidad",
                      reacciones=[reaccion], alimentacion=alimentacion,
                      fase="gas", P=P, P0=P0, unidad="bar",
                      temperaturas=[T_sol], modo_K="termo",
                      datos_termo=datos, csv_out=os.devnull,
                      fugacidad="lewis" if con_phi else None,
                      criticas=crit, convencion_Kphi=conv)
y_rx = rx["resultados"][float(T_sol)]["res"].mole_fractions
dif = max(abs(y_rx[sp] - y_sol[sp]) for sp in especies)
chk(dif < 1e-6 and rx["veredicto"] != "FAIL",
    f"contraste con run_exercise (misma T, fugacidad={'lewis' if con_phi else None}):"
    f" máx |Δy| = {dif:.1e}, selfcheck {rx['veredicto']}")
if clave is not None:
    d = abs(res_principal - clave)
    chk(d <= tol_clave, f"vs clave: {res_principal:.4g} vs {clave} "
                        f"(|Δ| = {d:.3g} ≤ {tol_clave})")
print(f"\nVEREDICTO: {'OK' if ok else 'REVISAR'}")

# ══════════════ AUTOVALIDACIÓN DEL SIMULACRO (cálculo independiente) ══════════
from core.convenciones import CONVENCION_KPHI as _CONV_DEF
_ok = True
def _chk(cond, txt):
    global _ok
    _ok &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {txt}")

def _K(T_):                      # ec. 13.18 (mismos datos C.1/C.4)
    return reaction_thermo({"CO": -1.0, "H2": -2.0, "CH4O": 1.0}, datos, T_)["K"]

def _phi(T_, P_, sp):            # Pitzer escrito a mano (sin core.nonideal)
    Tc_, Pc_, w_ = cargar_criticas([sp])[sp]
    Tr, Pr = T_ / Tc_, P_ / Pc_
    return exp(Pr / Tr * ((0.083 - 0.422 / Tr ** 1.6) + w_ * (0.139 - 0.172 / Tr ** 4.2)))

def _Kphi(T_, P_):
    return _phi(T_, P_, "CH4O") / (_phi(T_, P_, "CO") * _phi(T_, P_, "H2") ** 2)

_Ky = lambda x: x * (3 - 2 * x) ** 2 / (4 * (1 - x) ** 3)
print()
print("CHEQUEO A MANO (fórmulas del planteamiento, sin la plantilla):")
xi1 = brentq(lambda x: log(_Ky(x)) - log(_K(300.0)), 1e-12, 1 - 1e-12, xtol=1e-15)
y1 = xi1 / (3 - 2 * xi1)
T2 = brentq(lambda T_: log(_K(T_)) - log(27.0), 300, 1500, xtol=1e-10)
T3gi = brentq(lambda T_: log(_K(T_)) - log(27.0 / 100.0 ** 2), 300, 1500, xtol=1e-10)
assert _CONV_DEF == "maestra", "este chequeo a mano está escrito para la convención maestra"
T3 = brentq(lambda T_: log(_K(T_)) - log(27.0 / (100.0 ** 2 * _Kphi(T_, 100.0))),
            300, 1500, xtol=1e-10)
print(f"   1) K(300) = {_K(300.0):.5g} → ξ = {xi1:.6f} → y_CH3OH = {y1:.4f}")
print(f"   2) K_y(ξ = 0.75) = {_Ky(0.75):.4f} (= 27) → K(T) = 27 → T = {T2:.2f} K")
print(f"   3) gas ideal a 100 bar: K(T) = 27/100² → T = {T3gi:.2f} K")
print(f"      Lewis-Randall 'maestra': K(T) = 27/(100²·K_phi) → T = {T3:.2f} K"
      f"  (K_phi = {_Kphi(T3, 100.0):.4f})")
_chk(abs(_Ky(0.75) - 27.0) < 1e-12, "K_y para y = 0.5 vale 27 (planteamiento)")
if INCISO == "1":
    _chk(abs(y_sol["CH4O"] - y1) < 1e-8, f"plantilla y_CH3OH = {y_sol['CH4O']:.6f} = a mano {y1:.6f}")
    print(f"\nRESULTADO a1): y_CH3OH(300 K, 1 bar) = {y_sol['CH4O']:.4f}")
elif INCISO == "2":
    _chk(abs(T_sol - T2) < 1e-6, f"plantilla T = {T_sol:.4f} K = a mano {T2:.4f} K")
    print(f"\nRESULTADO a2): T(y = 0.5, 1 bar) = {T_sol:.1f} K")
else:
    _chk(abs(T_sol - T3) < 1e-6, f"plantilla T = {T_sol:.4f} K = a mano {T3:.4f} K")
    _chk(T3gi > T2, f"δ < 0: T(100 bar) = {T3gi:.1f} K > T(1 bar) = {T2:.1f} K (Le Chatelier)")
    _chk(Kphi_sol < 1, f"K_phi = {Kphi_sol:.4f} < 1 (φ_CH3OH < 1, φ_CO y φ_H2 > 1)")
    fuera = [sp for sp in nu if not det[sp]["valido"]]
    print(f"  [info] Vr < 2 (virial extrapolada, Fig. 3.14): {', '.join(fuera) or 'ninguna'}")
    print(f"\nRESULTADO a3): T(y = 0.5, 100 bar, Lewis-Randall '{conv}') = {T_sol:.1f} K"
          f"   (gas ideal {T3gi:.1f} K)")
print(f"VEREDICTO SIMULACRO: {'OK' if _ok else 'REVISAR'}")

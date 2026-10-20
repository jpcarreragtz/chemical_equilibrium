"""PLANTILLA — UNA reacción en fase gas con efecto de la presión y fugacidad
(2do parcial, Sesión 13 / Problema 1). Copia, edita SOLO el bloque DATOS:

    cp plantilla/plantilla_presion_fugacidad.py ejercicios/ejN_tema.py
    cd ejercicios && python3 ejN_tema.py

Modos (elige uno en `modo`):
  "y_a_T"     → composición de equilibrio a T y P dadas.
  "T_para_y"  → T a la que `especie_obj` llega a la fracción mol `y_obj`
                a la P dada (K_y queda fija por la composición objetivo).
Modelos (`modelo`):
  "gi"     → gas ideal, K_phi = 1.
  "lewis"  → solución ideal de gases (Lewis-Randall): phi^_i = phi_i de
             especie PURA a (T, P) por la virial de Pitzer
             (Tc, Pc, ω de data/criticas_b1.json). K_phi = Π phi_i^ν_i.

De qué lado entra K_phi lo decide ÚNICAMENTE core/convenciones.py:
    K = F · K_y · (P/P°)^δ ,  F = factor_Kphi(K_phi, convencion)
    "maestra"  (default): F = 1/K_phi  →  K_y = K·(P/P°)^(−δ)·K_phi
    "rigurosa" (SVA 13.28): F = K_phi  →  K_y = K·(P/P°)^(−δ)/K_phi
Aquí NO hay lógica de convención: solo se llama a core.convenciones.

Planteamiento que se imprime: ν y δ, n_i(ξ) y y_i(ξ), K_y(ξ), K(T) por la
ec. 13.18 (ΔH°, ΔG°, ΔA..ΔD), tabla de φ (Tr, Pr, B⁰, B¹, φ, Vr≥2), K_phi,
convención usada y K efectiva; luego AUTOVALIDACIÓN con PASS/FAIL (incluye
un contraste con run_exercise, el motor de plantilla_ejercicio.py).

Precargado: Problema 1 d) — ½N2 + 3/2 H2 = NH3, 100 bar, y_NH3 = 0.5,
Lewis-Randall → T = 570.4 K con la convención "maestra" (584.5 K rigurosa).
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # raíz

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
titulo = "Problema 1 d): síntesis de NH3 a 100 bar, y_NH3 = 0.5 (Lewis-Randall)"
reaccion = "1/2 N2 + 3/2 H2 = NH3"         # texto con coeficientes
alimentacion = {"N2": 0.5, "H2": 1.5}      # mol; los inertes también van aquí
P, P0 = 100.0, 1.0                         # bar (phi de Pitzer necesita bar)

modo = "T_para_y"                          # "y_a_T" | "T_para_y"
modelo = "lewis"                           # "gi" | "lewis"
convencion_Kphi = None                     # None = core.convenciones (default "maestra")

T = 300.0                                  # modo "y_a_T": temperatura [K]
especie_obj, y_obj = "NH3", 0.50           # modo "T_para_y": objetivo
T_lo, T_hi = 300.0, 1500.0                 # intervalo de búsqueda de T [K]

especies_termo = None     # None → las de la reacción (data/sva_tables.json)
criticas = None           # None → data/criticas_b1.json; o {sp: (Tc, Pc, ω)}
clave = 568.6             # valor esperado (K o fracción mol) o None
tol_clave = 2.5           # tolerancia absoluta para reportar PASS vs clave
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

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

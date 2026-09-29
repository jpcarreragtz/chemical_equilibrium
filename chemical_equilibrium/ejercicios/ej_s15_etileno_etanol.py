"""Sesión 15 — Equilibrio en sistemas heterogéneos (SVA Ej. 13.13 7a ed.)

Enunciado: estime las composiciones de las fases líquida y vapor cuando el
etileno reacciona con agua para formar etanol a 200 °C y 34.5 bar,
condiciones que aseguran la presencia de ambas fases. El recipiente se
mantiene a 34.5 bar por conexión a una fuente de etileno. No hay otras
reacciones:   C2H4(g) + H2O(g) -> C2H5OH(g)      (nu = -1, P° = 1 bar)

Modelo (Sesión 15):
* Regla de fases: F = 2 - pi + N - r = 2 - 2 + 3 - 1 = 2. Con T y P fijas
  las composiciones quedan DETERMINADAS: no hay tabla estequiométrica ni
  avance xi (la fuente de etileno solo repone lo que reacciona y mantiene
  P constante).
* Vapor: C2H4 + H2O + EtOH; líquido: SOLO H2O + EtOH (se desprecia la
  solubilidad del etileno en el líquido, como hace SVA).
* f^v_i = y_i·phi_i·P con phi_i de especie PURA a (T, P) — regla de
  Lewis-Randall (LEWIS_RANDALL = True), correlación virial de Pitzer.
* f^l_i = x_i·gamma_i·phi_i_sat·Psat_i, Poynting ~ 1 (se reporta cuánto
  vale pero no entra al cálculo base). gamma por Wilson (1=EtOH, 2=H2O).
* K = [x_E g_E phis_E Psat_E · P°] / [(y_C phi_C P)(x_W g_W phis_W Psat_W)]
* ELV (ec. 14): y_i = x_i g_i phis_i Psat_i/(phi_i P);  x_W = 1 - x_E;
  y_C = 1 - y_E - y_W.  Incógnita única: x_E; residual ln K_calc - ln K.
"""

import csv
import os
import pathlib
import sys
from math import exp, log

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto

import numpy as np
from scipy.optimize import brentq

from core.equilibrium import format_table
from core.nonideal import phi_pitzer, wilson_gamma, wilson_gamma_inf, wilson_lambdas
from core.selfcheck import Report
from core.thermo import R_J, T_REF, cargar_datos, icph, icps, K_vant_hoff, \
    delta_property, reaction_thermo

LEWIS_RANDALL = True     # phi^_i ~ phi_i(T,P) de especie pura (curso)

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
DATOS = {
    "T": 473.15,          # K  (200 °C)
    "P": 34.5,            # bar
    "P0": 1.0,            # bar (estado estándar gas ideal)
    # críticas (Tabla B.1 SVA):        Tc [K]   Pc [bar]  omega
    "criticas": {"C2H4": (282.3, 50.40, 0.087),
                 "H2O": (647.1, 220.55, 0.345),
                 "EtOH": (513.9, 61.48, 0.645)},
    # Psat a 200 °C [bar]: H2O de tablas de vapor; EtOH valor usado por SVA
    # (la Antoine B.2 del EtOH, rango 3-96 °C, extrapolada da ~32.6 bar ->
    #  se corre como sensibilidad)
    "Psat_H2O": 15.549,
    "Psat_EtOH": 30.22,
    "Psat_EtOH_alt": 32.6,
    # Wilson etanol(1)/agua(2), DECHEMA: a_ij [cal/mol], V_i [cm3/mol]
    "a12": 382.30, "a21": 955.45, "V1": 58.69, "V2": 18.07,
    # volúmenes líquidos para el factor de Poynting [cm3/mol]
    "Vl_H2O": 18.0, "Vl_EtOH": 58.7,
    # K de referencia de SVA/curso (ver nota en sección 2): sensibilidad
    "K_ref_SVA": 0.0323,
}
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

T, P, P0 = DATOS["T"], DATOS["P"], DATOS["P0"]

print("=" * 70)
print("  SESIÓN 15 — C2H4(g) + H2O(g) <-> EtOH(g) con fases L + V")
print("  T = 200 °C = 473.15 K,  P = 34.5 bar (fuente de etileno)")
print("=" * 70)

print("\n1) SUPUESTOS Y REGLA DE FASES")
print("   F = 2 - pi + N - r = 2 - 2 + 3 - 1 = 2  ->  con T y P fijas las")
print("   composiciones quedan determinadas: NO aplica tabla")
print("   estequiométrica ni avance de reacción xi.")
print("   Vapor: C2H4/H2O/EtOH. Líquido: solo H2O/EtOH (solubilidad de")
print("   C2H4 despreciada, como SVA). phi^_i ~ phi_i puro (Lewis-Randall,")
print(f"   LEWIS_RANDALL = {LEWIS_RANDALL}). Poynting ~ 1 (se reporta).")

# ---------------- 2) K por ec. 13.18 (thermo.py + tablas verificadas) -----
nu = {"C2H4": -1.0, "H2O": -1.0, "C2H6O": 1.0}   # etanol = C2H6O en el JSON
datos_t = cargar_datos(["C2H4", "H2O", "C2H6O"])
dA, dB, dC, dD = (delta_property(nu, datos_t, k) for k in "ABCD")
IDCPH = icph(T_REF, T, dA, dB, dC, dD)
IDCPS = icps(T_REF, T, dA, dB, dC, dD)
t473 = reaction_thermo(nu, datos_t, T)
t298 = reaction_thermo(nu, datos_t, T_REF)
K = t473["K"]
K_vh = K_vant_hoff(T, t298["K"], t298["dH298"])

print("\n2) CONSTANTE DE EQUILIBRIO (ec. 13.18, Tablas C.1/C.4 verificadas)")
print(f"   dH°298 = {t298['dH298']:.0f} J/mol    dG°298 = {t298['dG298']:.0f} "
      f"J/mol    K298 = {t298['K']:.4g}")
print(f"   dA = {dA:.3f}  dB = {dB:.4g}  dC = {dC:.4g}  dD = {dD:.4g}")
print(f"   IDCPH(298.15->473.15) = {IDCPH:.3f} K      IDCPS = {IDCPS:.6f}")
print(f"   ln K(473.15) = {log(K):.5f}   ->   K = {K:.5f}")
print(f"   Van't Hoff (dH cte): K = {K_vh:.5f}")
print(f"   [NOTA] La referencia SVA/curso usa K = {DATOS['K_ref_SVA']} "
      f"(ln K = {log(DATOS['K_ref_SVA']):.4f}); ese valor equivale a la")
print("   ec. 13.18 con el SIGNO de la corrección de Cp invertido y además")
print("   contradice Van't Hoff (reacción exotérmica con |dH| creciente =>")
print("   K(473) < K_vantHoff = 0.0317). Base = 0.03104; la fila de")
print("   sensibilidad 'K = 0.0323' reproduce los números del libro.")

# ---------------- 3) phi por la virial de Pitzer ---------------------------
esp_v = ["C2H4", "H2O", "EtOH"]
phi_TP = {sp: phi_pitzer(T, P, *DATOS["criticas"][sp]) for sp in esp_v}
phi_sat = {"H2O": phi_pitzer(T, DATOS["Psat_H2O"], *DATOS["criticas"]["H2O"]),
           "EtOH": phi_pitzer(T, DATOS["Psat_EtOH"],
                              *DATOS["criticas"]["EtOH"])}

print("\n3) COEFICIENTES DE FUGACIDAD (virial generalizada de Pitzer)")
filas = []
for sp in esp_v:
    r = phi_TP[sp]
    filas.append([f"{sp} @(T,P)", r["Tr"], r["Pr"], r["B0"], r["B1"],
                  r["phi"], "sí" if r["valido"] else "NO (extrapolada)"])
for sp in ("H2O", "EtOH"):
    r = phi_sat[sp]
    filas.append([f"{sp} @Psat", r["Tr"], r["Pr"], r["B0"], r["B1"],
                  r["phi"], "sí" if r["valido"] else "NO (extrapolada)"])
print("  " + format_table(filas, ["caso", "Tr", "Pr", "B0", "B1", "phi",
                                  "zona virial"]).replace("\n", "\n  "))

# ---------------- 4) Poynting (solo para reportar) --------------------------
poy = {sp: exp(V * 1e-6 * (P - Ps) * 1e5 / (R_J * T))
       for sp, V, Ps in (("H2O", DATOS["Vl_H2O"], DATOS["Psat_H2O"]),
                         ("EtOH", DATOS["Vl_EtOH"], DATOS["Psat_EtOH"]))}
print("\n4) FACTOR DE POYNTING (ec. 8, con V^l constante)")
print(f"   H2O: exp[18.0·(34.5-15.549)·1e-1/(8.314·473.15)] = "
      f"{poy['H2O']:.5f}")
print(f"   EtOH: exp[58.7·(34.5-30.22)·1e-1/(8.314·473.15)] = "
      f"{poy['EtOH']:.5f}")
print("   Ambos ~ 1.007-1.009 -> se justifica tomarlos = 1 en el cálculo")
print("   base (la fila 'con Poynting' de la sensibilidad los incluye).")

# ---------------- 5) Wilson -------------------------------------------------
L12, L21 = wilson_lambdas(T, DATOS["V1"], DATOS["V2"], DATOS["a12"],
                          DATOS["a21"])
g1_inf, g2_inf = wilson_gamma_inf(T, DATOS["V1"], DATOS["V2"], DATOS["a12"],
                                  DATOS["a21"])
print("\n5) WILSON etanol(1)/agua(2), DECHEMA, R = 1.987 cal/(mol·K)")
print(f"   L12 = (V2/V1)·exp(-a12/RT) = {L12:.5f}      "
      f"L21 = (V1/V2)·exp(-a21/RT) = {L21:.5f}")
print(f"   (dilución infinita: g_EtOH_inf = {g1_inf:.3f}, "
      f"g_H2O_inf = {g2_inf:.3f})")


# ---------------- modelo y solución -----------------------------------------
def estado(xE, prm):
    """Todo el estado para un x_EtOH dado, con los parámetros prm."""
    xW = 1.0 - xE
    if prm["gammas"]:
        gE, gW = wilson_gamma(xE, T, DATOS["V1"], DATOS["V2"], DATOS["a12"],
                              DATOS["a21"])
    else:
        gE = gW = 1.0
    fE = xE * gE * prm["phis_E"] * prm["psat_E"] * prm["poy_E"]  # f^l [bar]
    fW = xW * gW * prm["phis_W"] * prm["psat_W"] * prm["poy_W"]
    yE = fE / (prm["phi_E"] * P)                                  # ec. 14
    yW = fW / (prm["phi_W"] * P)
    yC = 1.0 - yE - yW                                            # ec. 20
    return {"xE": xE, "xW": xW, "gE": gE, "gW": gW, "fE": fE, "fW": fW,
            "yE": yE, "yW": yW, "yC": yC}


def lnK_calc(e, prm):
    return log(e["fE"] * P0 / ((e["yC"] * prm["phi_C"] * P) * e["fW"]))


def resolver(prm, x_lo=1e-6, frac_hi=0.999):
    """brentq sobre F(xE) = lnK_calc - lnK en la región con y_C > 0."""
    malla = np.linspace(1e-6, 0.999, 4000)
    con_gas = [x for x in malla if estado(x, prm)["yC"] > 0.0]
    x_max = max(con_gas)
    lo, hi = x_lo, frac_hi * x_max
    F = lambda x: lnK_calc(estado(x, prm), prm) - log(prm["K"])
    Flo, Fhi = F(lo), F(hi)
    if Flo * Fhi > 0:
        raise RuntimeError(f"Sin cambio de signo en [{lo:.3g}, {hi:.3g}]: "
                           f"F = {Flo:.3g}, {Fhi:.3g}")
    xE = brentq(F, lo, hi, xtol=1e-14, rtol=8.9e-16)
    e = estado(xE, prm)
    e["residual"] = lnK_calc(e, prm) - log(prm["K"])
    e["x_max"] = x_max
    return e


BASE = {"K": K, "gammas": True,
        "phi_C": phi_TP["C2H4"]["phi"], "phi_E": phi_TP["EtOH"]["phi"],
        "phi_W": phi_TP["H2O"]["phi"], "phis_E": phi_sat["EtOH"]["phi"],
        "phis_W": phi_sat["H2O"]["phi"], "psat_E": DATOS["Psat_EtOH"],
        "psat_W": DATOS["Psat_H2O"], "poy_E": 1.0, "poy_W": 1.0}

sol = resolver(BASE)

print("\n6) SOLUCIÓN (incógnita única x_EtOH; brentq, residual logarítmico)")
print(f"   x_EtOH = {sol['xE']:.5f}    x_H2O = {sol['xW']:.5f}")
print(f"   gamma_EtOH = {sol['gE']:.4f}    gamma_H2O = {sol['gW']:.4f}")
print(f"   y_C2H4 = {sol['yC']:.5f}    y_H2O = {sol['yW']:.5f}    "
      f"y_EtOH = {sol['yE']:.5f}")
print(f"   residual |ln K_calc - ln K| = {abs(sol['residual']):.2e}")
print("   Fugacidades [bar] (isofugacidad L = V para H2O y EtOH):")
fil = [["EtOH", sol["fE"], sol["yE"] * phi_TP["EtOH"]["phi"] * P],
       ["H2O", sol["fW"], sol["yW"] * phi_TP["H2O"]["phi"] * P],
       ["C2H4", "-", sol["yC"] * phi_TP["C2H4"]["phi"] * P]]
print("  " + format_table(fil, ["Especie", "f^ liq (bar)", "f^ vap (bar)"]
                          ).replace("\n", "\n  "))

# ---------------- 7) sensibilidad --------------------------------------------
def variante(nombre, **cambios):
    prm = dict(BASE)
    prm.update(cambios)
    e = resolver(prm)
    return [nombre, e["xE"], e["xW"], e["yC"], e["yW"], e["yE"],
            e["gE"], e["gW"]]


sens = [variante("base (K = 0.03104)"),
        variante("gamma = 1 (líq. ideal)", gammas=False),
        variante("phi = 1 en todo (gas ideal)", phi_C=1.0, phi_E=1.0,
                 phi_W=1.0, phis_E=1.0, phis_W=1.0),
        variante("Psat_EtOH = 32.6 bar (Antoine extrapolada)",
                 psat_E=DATOS["Psat_EtOH_alt"],
                 phis_E=phi_pitzer(T, DATOS["Psat_EtOH_alt"],
                                   *DATOS["criticas"]["EtOH"])["phi"]),
        variante("con Poynting", poy_E=poy["EtOH"], poy_W=poy["H2O"]),
        variante(f"K = {DATOS['K_ref_SVA']} (valor SVA/curso)",
                 K=DATOS["K_ref_SVA"])]
CAB = ["caso", "x_EtOH", "x_H2O", "y_C2H4", "y_H2O", "y_EtOH",
       "g_EtOH", "g_H2O"]
print("\n7) SENSIBILIDAD")
print("  " + format_table(sens, CAB).replace("\n", "\n  "))

# ---------------- 8) checks automáticos --------------------------------------
print("\n8) CHECKS AUTOMÁTICOS")
rep = Report()
sy = sol["yC"] + sol["yW"] + sol["yE"]
sx = sol["xE"] + sol["xW"]
fracs = [sol[k] for k in ("xE", "xW", "yC", "yW", "yE")]
rep.add("pass" if abs(sy - 1) < 1e-10 and abs(sx - 1) < 1e-10
        and all(0 < f < 1 for f in fracs) else "fail", 1,
        f"sum y = {sy:.12f}, sum x = {sx:.12f}, fracciones en (0,1)")
rep.add("pass" if abs(sol["residual"]) < 1e-8 else "fail", 2,
        f"residual |ln K_calc - ln K| = {abs(sol['residual']):.2e} < 1e-8")
iso = []
for sp, fl, y, phv in (("EtOH", sol["fE"], sol["yE"], phi_TP["EtOH"]["phi"]),
                       ("H2O", sol["fW"], sol["yW"], phi_TP["H2O"]["phi"])):
    fv = y * phv * P
    iso.append(abs(fl - fv) / fv)
rep.add("pass" if max(iso) < 1e-8 else "fail", 3,
        f"isofugacidad: |f^l - f^v|/f^v = EtOH {iso[0]:.1e}, "
        f"H2O {iso[1]:.1e}")
K_fug = (sol["fE"] * P0) / ((sol["yC"] * phi_TP["C2H4"]["phi"] * P)
                            * sol["fW"])
rep.add("pass" if abs(K_fug - K) / K < 1e-8 else "fail", 4,
        f"K desde fugacidades = {K_fug:.6f} vs K = {K:.6f} "
        f"(dif rel {abs(K_fug - K) / K:.1e})")
g1_puro, _ = wilson_gamma(1.0, T, DATOS["V1"], DATOS["V2"], DATOS["a12"],
                          DATOS["a21"])
_, g2_puro = wilson_gamma(0.0, T, DATOS["V1"], DATOS["V2"], DATOS["a12"],
                          DATOS["a21"])
rep.add("pass" if abs(g1_puro - 1) < 1e-12 and abs(g2_puro - 1) < 1e-12
        else "fail", 5, f"Wilson: g1(x1=1) - 1 = {g1_puro - 1:.1e}, "
        f"g2(x2=1) - 1 = {g2_puro - 1:.1e}")
raices = [resolver(BASE, x_lo=lo, frac_hi=fh)["xE"]
          for lo, fh in ((1e-6, 0.999), (1e-5, 0.9), (1e-4, 0.8),
                         (0.01, 0.6), (0.02, 0.5))]
spread = max(raices) - min(raices)
rep.add("pass" if spread < 1e-10 else "fail", 6,
        f"robustez: 5 intervalos distintos -> misma raíz "
        f"(dispersión = {spread:.1e})")
rep.add("pass", 7, "regla de fases: F = 2 con (T, P) dadas -> composiciones "
                   "fijas; xi NO aplica")
veredicto = rep.verdict()

# ---------------- salidas ----------------------------------------------------
out = os.path.join(_DIR, "entrega_s15")
os.makedirs(out, exist_ok=True)
ruta_csv = os.path.join(out, "sensibilidad_s15.csv")
with open(ruta_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(CAB)
    for fila in sens:
        w.writerow([fila[0]] + [f"{v:.5f}" for v in fila[1:]])

ruta_md = os.path.join(out, "RESUMEN_s15.md")
with open(ruta_md, "w", encoding="utf-8") as f:
    f.write(f"""# Sesión 15 — Hidratación de etileno con fases L+V (SVA Ej. 13.13)

**Condiciones:** T = 200 °C (473.15 K), P = 34.5 bar (fuente de etileno).
**Regla de fases:** F = 2−2+3−1 = 2 → con T y P fijas las composiciones
quedan determinadas; **no aplica ξ** ni tabla estequiométrica.
**Supuestos:** líquido sin C2H4 (SVA); Lewis-Randall (φ̂ᵢ ≈ φᵢ puro,
virial de Pitzer); Poynting ≈ 1 (H2O {poy['H2O']:.4f}, EtOH {poy['EtOH']:.4f});
Wilson etanol(1)/agua(2) DECHEMA.

## Constante de equilibrio (ec. 13.18, tablas verificadas)
ΔH°298 = {t298['dH298']:.0f} J/mol, ΔG°298 = {t298['dG298']:.0f} J/mol,
K298 = {t298['K']:.4g}; IDCPH = {IDCPH:.3f} K, IDCPS = {IDCPS:.6f};
**K(473.15) = {K:.5f}** (ln K = {log(K):.4f}); Van't Hoff ΔH cte: {K_vh:.4f}.
> Nota: el valor SVA/curso 0.0323 equivale a invertir el signo de la
> corrección de Cp y contradice Van't Hoff (K debe ser < 0.0317 por ser
> exotérmica con |ΔH| creciente); se incluye como sensibilidad y ahí se
> reproducen los números del libro.

## Coeficientes de fugacidad (Pitzer)
φ(T,P): C2H4 {phi_TP['C2H4']['phi']:.4f}, H2O {phi_TP['H2O']['phi']:.4f},
EtOH {phi_TP['EtOH']['phi']:.4f} (EtOH fuera de la zona virial → extrapolado);
φsat: H2O {phi_sat['H2O']['phi']:.4f}, EtOH {phi_sat['EtOH']['phi']:.4f}.
Wilson: Λ12 = {L12:.4f}, Λ21 = {L21:.4f}.

## Solución (base, K = {K:.5f})
| | EtOH | H2O | C2H4 |
|---|---|---|---|
| x (líquido) | {sol['xE']:.4f} | {sol['xW']:.4f} | — |
| γ | {sol['gE']:.3f} | {sol['gW']:.3f} | — |
| y (vapor) | {sol['yE']:.4f} | {sol['yW']:.4f} | {sol['yC']:.4f} |
| f̂ [bar] | {sol['fE']:.3f} | {sol['fW']:.3f} | {sol['yC'] * phi_TP['C2H4']['phi'] * P:.3f} |

Residual |ln K_calc − ln K| = {abs(sol['residual']):.1e}. Checks: {veredicto}.

## Sensibilidad (CSV: sensibilidad_s15.csv)
""")
    f.write("| " + " | ".join(CAB) + " |\n")
    f.write("|" + "---|" * len(CAB) + "\n")
    for fila in sens:
        f.write("| " + fila[0] + " | "
                + " | ".join(f"{v:.4f}" for v in fila[1:]) + " |\n")

print(f"\n   Resumen: {ruta_md}")
print(f"   CSV sensibilidad: {ruta_csv}")

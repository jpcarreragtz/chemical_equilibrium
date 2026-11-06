"""2do parcial práctico (Grupo B, 22/10/2026) — PROBLEMA 1
Deshidrogenación de etilbenceno a estireno

    C8H10(g) -> C8H8(g) + H2(g)          alimentación: 2 mol/s de C8H10 puro

a) X y fracciones mol a 800 K y 1 bar                         (gas ideal)
b) T a la que X = 0.80 a 15 bar                  (SOLUCIÓN IDEAL de gases)
c) P a la que X = 0.80 a 1100 K                  (SOLUCIÓN IDEAL de gases)

Planteamiento (a mano)
----------------------
ν = {C8H10: −1, C8H8: +1, H2: +1}   →   δ = Σν = +1
Tabla en X (ξ = 2X mol/s):
    n_EB = 2(1 − X)    n_ST = 2X    n_H2 = 2X    n_T = 2(1 + X)
    y_EB = (1 − X)/(1 + X)      y_ST = y_H2 = X/(1 + X)
    K_y = y_ST·y_H2 / y_EB = X²/[(1 + X)(1 − X)] = X²/(1 − X²)
    ⇒ X = sqrt(K_y / (1 + K_y))   (solución analítica)

    K = F · K_y · (P/P°)^δ ,   K(T) por la ec. 13.18 de SVA (Tablas C.1/C.4)
    gas ideal:   F = 1
    sol. ideal:  K_phi = φ_ST·φ_H2/φ_EB  (Lewis-Randall: φ_i de especie PURA
                 a T y P, virial truncada de Pitzer
                 ln φ = (Pr/Tr)(B0 + ω·B1), B0 = 0.083 − 0.422/Tr^1.6,
                 B1 = 0.139 − 0.172/Tr^4.2);  F según core/convenciones.py:
        "maestra"  (default): F = 1/K_phi  →  K_y = K·(P/P°)^(−δ)·K_phi
        "rigurosa" (SVA 13.28): F = K_phi  →  K_y = K·(P/P°)^(−δ)/K_phi

    X = 0.80  →  K_y = 0.64/0.36 = 1.7778 (no depende de T ni de P)
    b) brentq sobre T:  ln K(T) − [ln F(T, 15) + ln K_y + ln 15] = 0
    c) brentq sobre P:  ln K(1100) − [ln F(1100, P) + ln K_y + ln P] = 0
       (gas ideal: P = K(1100)/K_y directo)
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz

from math import exp, log, sqrt
from scipy.optimize import brentq

from core.convenciones import (CONVENCION_KPHI, CONVENCIONES, describir,
                               factor_Kphi, otra_convencion)
from core.equilibrium import format_table
from core.nonideal import K_phi_lewis, cargar_criticas
from core.thermo import cargar_datos, delta_property, reaction_thermo

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
NU = {"C8H10": -1.0, "C8H8": 1.0, "H2": 1.0}
F_EB0 = 2.0                     # mol/s de etilbenceno alimentado
P0 = 1.0                        # bar (estado estándar)
TA, PA = 800.0, 1.0             # inciso a)
X_OBJ = 0.80                    # incisos b) y c)
PB = 15.0                       # inciso b): P [bar]
TC_ = 1100.0                    # inciso c): T [K]
ESPERADO = {"dH298": 117440.0, "dG298": 83010.0,
            "dA": 4.175, "dB": -4.766e-3, "dC": 1.814e-6, "dD": 8300.0,
            "K800": 0.05701, "Xa": 0.2322, "y_EB": 0.6231, "y_ST": 0.1885,
            "Tb": 1174.0, "Tb_gi": 1173.97, "K1100": 10.879, "Pc": 6.12}
CONVENCION = None       # None → core.convenciones.CONVENCION_KPHI ("maestra")
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

CONV = CONVENCION or CONVENCION_KPHI
OTRA = otra_convencion(CONV)
DELTA = sum(NU.values())
datos = cargar_datos(list(NU))
crit = cargar_criticas(list(NU))
NOMBRE = {"C8H10": "EB", "C8H8": "ST", "H2": "H2"}


def comp(X):
    """Flujos [mol/s] y fracciones mol en función de X (tabla en X)."""
    n = {"C8H10": F_EB0 * (1 - X), "C8H8": F_EB0 * X, "H2": F_EB0 * X}
    nT = sum(n.values())
    return n, {sp: v / nT for sp, v in n.items()}, nT


def Ky_de_X(X):
    return X * X / (1 - X * X)


def X_de_Ky(ky):
    return sqrt(ky / (1 + ky))


def K_T(T):
    return reaction_thermo(NU, datos, T)["K"]


def Kphi(T, P):
    return K_phi_lewis(NU, T, P, crit, avisar=False)[0]


def F_de(T, P, modelo, convencion):
    return factor_Kphi(Kphi(T, P), convencion) if modelo == "lewis" else 1.0


def X_equilibrio(T, P, modelo="gi", convencion=None):
    """X de equilibrio: K_y = K/(F·(P/P°)^δ) y X = sqrt(K_y/(1+K_y))."""
    ky = K_T(T) / (F_de(T, P, modelo, convencion) * (P / P0) ** DELTA)
    return X_de_Ky(ky), ky


def T_para_X(X, P, modelo="gi", convencion=None, T_lo=600.0, T_hi=1500.0):
    ky = Ky_de_X(X)
    f = lambda T: (log(K_T(T)) - log(F_de(T, P, modelo, convencion))
                   - log(ky) - DELTA * log(P / P0))
    return brentq(f, T_lo, T_hi, xtol=1e-10)


def P_para_X(X, T, modelo="gi", convencion=None, P_lo=1e-3, P_hi=200.0):
    """Modo 'P_para_X' (no está en la plantilla): brentq sobre P."""
    ky = Ky_de_X(X)
    f = lambda P: (log(K_T(T)) - log(F_de(T, P, modelo, convencion))
                   - log(ky) - DELTA * log(P / P0))
    return brentq(f, P_lo, P_hi, xtol=1e-12)


def tabla_phi(T, P):
    Kp, det = K_phi_lewis(NU, T, P, crit, avisar=False)
    filas = [[f"{sp} ({NOMBRE[sp]})", *crit[sp], det[sp]["Tr"], det[sp]["Pr"],
              det[sp]["B0"], det[sp]["B1"], det[sp]["phi"], det[sp]["Vr"],
              "sí" if det[sp]["valido"] else "NO*"] for sp in NU]
    print("  " + format_table(filas, ["especie", "Tc", "Pc", "ω", "Tr", "Pr",
                                      "B0", "B1", "φ", "Vr", "Vr≥2"]
                              ).replace("\n", "\n  "))
    print(f"   K_phi = φ_ST·φ_H2/φ_EB = {Kp:.6f}")
    return Kp, det


print("=" * 72)
print("  PROBLEMA 1 — C8H10(g) -> C8H8(g) + H2(g)   (2 mol/s de etilbenceno)")
print("=" * 72)

# ---------------------------------------------------------------- 0) termo
t298 = reaction_thermo(NU, datos, 298.15)
dA, dB, dC, dD = (delta_property(NU, datos, k) for k in "ABCD")
print("\n0) DATOS TERMO (Tablas C.1/C.4 SVA, data/sva_tables.json)")
filas = [[sp, datos[sp]["Hf298"], datos[sp]["Gf298"], datos[sp]["A"],
          datos[sp]["B"], datos[sp]["C"], datos[sp]["D"]] for sp in NU]
print("  " + format_table(filas, ["especie", "ΔH°f J/mol", "ΔG°f J/mol", "A",
                                  "B", "C", "D"]).replace("\n", "\n  "))
print(f"   ΔH°298 = 147 360 + 0 − 29 920 = {t298['dH298']:.0f} J/mol")
print(f"   ΔG°298 = 213 900 + 0 − 130 890 = {t298['dG298']:.0f} J/mol")
print(f"   ΔA = {dA:.4f}  ΔB = {dB:.4e}  ΔC = {dC:.4e}  ΔD = {dD:.0f}")
print(f"   δ = Σν = {DELTA:+.0f}  → K = F·K_y·(P/P°)^δ   (F = 1 en gas ideal)")
print(f"   K_y = X²/(1 − X²)  ⇔  X = sqrt(K_y/(1 + K_y))")

# ---------------------------------------------------------------- a)
Ka = K_T(TA)
Xa, kya = X_equilibrio(TA, PA)
na, ya, nTa = comp(Xa)
Kphi_a = Kphi(TA, PA)
print(f"\na) T = {TA:g} K, P = {PA:g} bar, GAS IDEAL")
print(f"   ΔG°({TA:g}) = {reaction_thermo(NU, datos, TA)['dG']:.1f} J/mol"
      f"  →  K({TA:g}) = {Ka:.5f}   (ln K = {log(Ka):.5f})")
print(f"   K_y = K·(P/P°)^(−1) = {kya:.5f}")
print(f"   X = sqrt({kya:.5f}/(1 + {kya:.5f})) = {Xa:.5f}")
filas = [[f"{sp} ({NOMBRE[sp]})", F_EB0 if sp == "C8H10" else 0.0, na[sp],
          ya[sp]] for sp in NU]
filas.append(["TOTAL", F_EB0, nTa, sum(ya.values())])
print("  " + format_table(filas, ["especie", "entra mol/s", "sale mol/s",
                                  "y_i"]).replace("\n", "\n  "))
print(f"   → X = {Xa:.4f};  y_EB = {ya['C8H10']:.4f},  y_ST = y_H2 = "
      f"{ya['C8H8']:.4f}")
print(f"   K_phi(800 K, 1 bar) = {Kphi_a:.5f} ≈ 1  → a 1 bar el gas ideal es "
      f"adecuado")

# ---------------------------------------------------------------- b)
ky_obj = Ky_de_X(X_OBJ)
_, y_obj, _ = comp(X_OBJ)
print(f"\n   Para X = {X_OBJ}: y_EB = {y_obj['C8H10']:.5f}, y_ST = y_H2 = "
      f"{y_obj['C8H8']:.5f}, K_y = 0.64/0.36 = {ky_obj:.5f}")

Tb_gi = T_para_X(X_OBJ, PB)
Tb_por = {c: T_para_X(X_OBJ, PB, "lewis", c) for c in CONVENCIONES}
Tb = Tb_por[CONV]
print(f"\nb) X = {X_OBJ}, P = {PB:g} bar, SOLUCIÓN IDEAL (Lewis-Randall)")
print(f"   Convención '{CONV}' (core/convenciones.py):")
print(f"   {describir(CONV)}")
print(f"   gas ideal: K requerida = K_y·(P/P°) = {ky_obj * PB:.4f} → "
      f"T = {Tb_gi:.2f} K")
print(f"   Tabla de φ a T = {Tb:.2f} K, P = {PB:g} bar:")
Kphi_b, det_b = tabla_phi(Tb, PB)
F_b = factor_Kphi(Kphi_b, CONV)
print(f"   F = {F_b:.6f}  →  K requerida = F·K_y·(P/P°) = "
      f"{F_b * ky_obj * PB:.4f} = K({Tb:.2f} K) = {K_T(Tb):.4f}")
print(f"   → T = {Tb:.2f} K  ({Tb - 273.15:.2f} °C)   ['{OTRA}': "
      f"{Tb_por[OTRA]:.2f} K]")

# ---------------------------------------------------------------- c)
Kc = K_T(TC_)
Pc_gi = Kc / ky_obj
Pc_por = {c: P_para_X(X_OBJ, TC_, "lewis", c) for c in CONVENCIONES}
Pc = Pc_por[CONV]
print(f"\nc) X = {X_OBJ}, T = {TC_:g} K, SOLUCIÓN IDEAL (Lewis-Randall)")
print(f"   K({TC_:g}) = {Kc:.4f}   (ec. 13.18)")
print(f"   gas ideal: P/P° = K/K_y = {Kc:.4f}/{ky_obj:.5f} → P = {Pc_gi:.4f} bar")
print(f"   Tabla de φ a T = {TC_:g} K, P = {Pc:.4f} bar:")
Kphi_c, det_c = tabla_phi(TC_, Pc)
F_c = factor_Kphi(Kphi_c, CONV)
print(f"   F = {F_c:.6f}  →  P/P° = K/(F·K_y) = {Kc / (F_c * ky_obj):.4f}")
print(f"   → P = {Pc:.4f} bar   ['{OTRA}': {Pc_por[OTRA]:.4f} bar]")
print("   (Kφ ≈ 1 porque las especies están a Tr alta y Pr baja: el efecto de "
      "la\n    no idealidad es despreciable en b) y c))")

# ---------------------------------------------------------------- checks
print("\nAUTOVALIDACIÓN")
ok = True


def chk(cond, txt):
    global ok
    ok &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {txt}")


def rel(v, ref):
    return abs(v - ref) / abs(ref)


chk(abs(t298["dH298"] - ESPERADO["dH298"]) < 0.5
    and abs(t298["dG298"] - ESPERADO["dG298"]) < 0.5,
    f"ΔH°298 = {t298['dH298']:.0f}, ΔG°298 = {t298['dG298']:.0f} J "
    f"(esperado {ESPERADO['dH298']:.0f}, {ESPERADO['dG298']:.0f})")
chk(all(rel(v, ESPERADO[k]) < 1e-3 for k, v in
        (("dA", dA), ("dB", dB), ("dC", dC), ("dD", dD))),
    f"ΔA, ΔB, ΔC, ΔD = {dA:.3f}, {dB:.4e}, {dC:.4e}, {dD:.0f}")
chk(rel(Ka, ESPERADO["K800"]) < 1e-3,
    f"a) K(800) = {Ka:.5f} vs {ESPERADO['K800']}")
chk(rel(Xa, ESPERADO["Xa"]) < 1e-3 and rel(ya["C8H10"], ESPERADO["y_EB"]) < 1e-3
    and rel(ya["C8H8"], ESPERADO["y_ST"]) < 1e-3,
    f"a) X = {Xa:.4f}, y_EB = {ya['C8H10']:.4f}, y_ST = {ya['C8H8']:.4f} "
    f"(esperado {ESPERADO['Xa']}, {ESPERADO['y_EB']}, {ESPERADO['y_ST']})")
chk(abs(sum(ya.values()) - 1) < 1e-12
    and abs(ya["C8H8"] * ya["H2"] / ya["C8H10"] * PA / P0 - Ka) < 1e-12,
    f"a) Σy = 1 y sustituyendo: y_ST·y_H2/y_EB·(P/P°) = K")
at_in = {"C": 8 * F_EB0, "H": 10 * F_EB0}
at_out = {"C": 8 * na["C8H10"] + 8 * na["C8H8"],
          "H": 10 * na["C8H10"] + 8 * na["C8H8"] + 2 * na["H2"]}
chk(all(abs(at_in[e] - at_out[e]) < 1e-12 for e in at_in),
    f"a) balance atómico: C {at_out['C']:.6f} ({at_in['C']:g}), "
    f"H {at_out['H']:.6f} ({at_in['H']:g}) átomo-mol/s")
chk(abs(Kphi_a - 1) < 0.01, f"a) K_phi(800 K, 1 bar) = {Kphi_a:.5f} ≈ 1 "
    f"(|K_phi − 1| < 0.01)")
for c in CONVENCIONES:
    Xb = X_equilibrio(Tb_por[c], PB, "lewis", c)[0]
    Xc = X_equilibrio(TC_, Pc_por[c], "lewis", c)[0]
    chk(abs(Xb - X_OBJ) < 1e-9 and abs(Xc - X_OBJ) < 1e-9,
        f"'{c}': re-resolviendo X en b) ({Tb_por[c]:.2f} K) = {Xb:.10f} y "
        f"en c) ({Pc_por[c]:.4f} bar) = {Xc:.10f}")
chk(abs(Tb - ESPERADO["Tb"]) < 0.1 and abs(Tb_gi - ESPERADO["Tb_gi"]) < 0.01,
    f"b) T = {Tb:.2f} K (esperado {ESPERADO['Tb']}); gas ideal {Tb_gi:.2f} K "
    f"(esperado {ESPERADO['Tb_gi']})")
chk(abs(det_b["C8H10"]["phi"] - 0.9936) < 5e-5
    and abs(det_b["C8H8"]["phi"] - 0.9920) < 5e-5
    and abs(det_b["H2"]["phi"] - 1.0017) < 5e-5,
    f"b) φ_EB = {det_b['C8H10']['phi']:.4f}, φ_ST = {det_b['C8H8']['phi']:.4f}"
    f", φ_H2 = {det_b['H2']['phi']:.4f} (esperado 0.9936, 0.9920, 1.0017)")
chk(rel(Kc, ESPERADO["K1100"]) < 1e-3 and abs(Pc - ESPERADO["Pc"]) < 0.005,
    f"c) K(1100) = {Kc:.4f}, P = {Pc:.4f} bar (esperado {ESPERADO['K1100']},"
    f" {ESPERADO['Pc']})")
chk(abs(Kphi_b - 1) < 1e-3 and abs(Kphi_c - 1) < 1e-3,
    f"K_phi b) = {Kphi_b:.6f}, c) = {Kphi_c:.6f} (≈ 1)")
chk(all(d["valido"] for d in (*det_b.values(), *det_c.values())),
    "Vr ≥ 2 para todas las especies en b) y c) → virial truncada válida")
dT = abs(Tb_por["maestra"] - Tb_por["rigurosa"])
dP = abs(Pc_por["maestra"] - Pc_por["rigurosa"])
chk(dT < 0.01 and dP < 0.001,
    f"convenciones coinciden: |ΔT| b) = {dT:.4f} K (< 0.01), |ΔP| c) = "
    f"{dP:.5f} bar (< 0.001)")

print(f"\n>>> RESULTADOS (convención '{CONV}')")
print(f"    a) X = {Xa:.4f};  y_EB = {ya['C8H10']:.4f}, y_ST = y_H2 = "
      f"{ya['C8H8']:.4f};  salida: EB {na['C8H10']:.4f}, ST {na['C8H8']:.4f},"
      f" H2 {na['H2']:.4f} mol/s")
print(f"    b) T = {Tb:.1f} K ({Tb - 273.15:.1f} °C)")
print(f"    c) P = {Pc:.2f} bar")
print(f"\nVEREDICTO: {'OK' if ok else 'REVISAR'}")

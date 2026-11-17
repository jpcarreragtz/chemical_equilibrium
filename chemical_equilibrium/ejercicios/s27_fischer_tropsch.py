"""Sesión 27 — Fischer–Tropsch por minimización de Gibbs a 500 K y 25 bar.

Especies: parafinas C1..C15 (CkH2k+2), CO2, H2O, H2, CO. ΔG a 500 K y φ de la
maestra (data/ft_s27_maestra.json, "Sesión 27, diapositiva Datos").
Datos de la maestra (Sesión 27), confirmados por ella el 17-nov-2026; se usan tal cual.
  Análisis 1: CH4, CO2, H2O y suma C5–C15, con y sin ASF (α = 0.80).
  Análisis 2: lo mismo, con y sin φ en la función objetivo.
  Análisis 3: b de biomasa (C=1, H=2, O=1; H2/CO = 1) y de gas natural
              (C=1, H=5, O=1; H2/CO = 2.5).
  Opcional:   % en peso C1 / C2–C4 / C5–C11 / C12–C15 con PM = 12k + (2k+2).
ASF (diapositiva, E6_ASF): n(Ck) = n(C1)·α^(k−1), k = 2..15, como 14
restricciones lineales extra de A n = b (core.gibbs, `restricciones`).

SUPUESTO: los Análisis 1 y 2 no dicen qué b usar; se toma H2/CO = 2
(C=1, H=4, O=1) en B_BASE. Cámbialo ahí.
"""

import contextlib
import io
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto
import core.entorno  # noqa: F401

import numpy as np

from core import DATA
from core.gibbs import minimizar_gibbs, prueba_de_fase, validar_gibbs

SI_NO = {True: "sí", False: "no"}   # etiquetas; compatible con Python < 3.12

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
MAESTRA = json.load(open(DATA / "ft_s27_maestra.json", encoding="utf-8"))
T = float(MAESTRA["_meta"]["T_K"])        # 500 K
P = float(MAESTRA["_meta"]["P_bar"])      # 25 bar (P0 = 1 bar)
ALFA = 0.80                                # ASF
B_BASE = {"C": 1.0, "H": 4.0, "O": 1.0}    # SUPUESTO (H2/CO = 2) de los Análisis 1 y 2
B_BIOMASA = {"C": 1.0, "H": 2.0, "O": 1.0}     # H2/CO = 1
B_GAS_NATURAL = {"C": 1.0, "H": 5.0, "O": 1.0}  # H2/CO = 2.5

ESPECIES = list(MAESTRA["especies"])
DGF = {sp: float(d["dG_500K"]) for sp, d in MAESTRA["especies"].items()}
PHI = {sp: float(d["phi"]) for sp, d in MAESTRA["especies"].items()}
ETIQ = {sp: d["etiqueta"] for sp, d in MAESTRA["especies"].items()}
PARAFINAS = [sp for sp in ESPECIES if ETIQ[sp].startswith("C") and ETIQ[sp][1:].isdigit()]
NCARB = {sp: int(ETIQ[sp][1:]) for sp in PARAFINAS}
ALIMENTACIONES = {"base H2/CO=2 (supuesto)": B_BASE, "biomasa H2/CO=1": B_BIOMASA,
                  "gas natural H2/CO=2.5": B_GAS_NATURAL}


def titulo(t: str) -> None:
    print(f"\n{'=' * 84}\n{t}\n{'=' * 84}")


def alimentacion_de_b(b):
    """b (átomos C, H, O) -> CO + H2. Exige C = O, que es lo que da CO + H2."""
    if abs(b["C"] - b["O"]) > 1e-12:
        raise ValueError("Con CO + H2 la alimentación tiene C = O")
    return {"CO": b["C"], "H2": b["H"] / 2.0}


def restricciones_asf(alfa):
    """n(Ck) − α^(k−1)·n(C1) = 0 para k = 2..15."""
    c1 = PARAFINAS[0]
    return [({sp: 1.0, c1: -alfa ** (NCARB[sp] - 1)}, 0.0) for sp in PARAFINAS[1:]]


def pct_peso(n):
    """% en peso de C1 / C2–C4 / C5–C11 / C12–C15 entre las parafinas."""
    masa = {sp: n[sp] * (12 * NCARB[sp] + (2 * NCARB[sp] + 2)) for sp in PARAFINAS}
    tot = sum(masa.values())
    grupos = {"C1": (1, 1), "C2–C4": (2, 4), "C5–C11": (5, 11), "C12–C15": (12, 15)}
    return {g: 100 * sum(m for sp, m in masa.items() if a <= NCARB[sp] <= z) / tot
            for g, (a, z) in grupos.items()}


def pct_peso_asf_analitico(alfa):
    """Con ASF, n_k ∝ α^(k−1): el % en peso no depende de la alimentación ni de φ."""
    return pct_peso({sp: alfa ** (NCARB[sp] - 1) for sp in PARAFINAS})


# %% 1. Las 12 corridas: 3 alimentaciones × ASF sí/no × φ sí/no
RES, VALID = {}, {}
for nombre, b in ALIMENTACIONES.items():
    for asf in (False, True):
        for con_phi in (False, True):
            r = minimizar_gibbs(ESPECIES, DGF, alimentacion_de_b(b), T, P,
                                phi=PHI if con_phi else None,
                                restricciones=restricciones_asf(ALFA) if asf else None)
            with contextlib.redirect_stdout(io.StringIO()):
                VALID[(nombre, asf, con_phi)] = validar_gibbs(r, n_arranques=5)
            RES[(nombre, asf, con_phi)] = r


def fila(clave):
    r = RES[clave]
    c5 = sum(r.n[sp] for sp in PARAFINAS if NCARB[sp] >= 5)
    vals = (r.n["CH4"], r.n["CO2"], r.n["H2O"], c5, r.y["CH4"], r.y["CO2"],
            r.y["H2O"], c5 / r.n_total, r.n["H2"], r.n["CO"])
    return "".join(f" {v:>10.5g}" for v in vals)


ENC = f"{'':<40}" + "".join(f" {c:>10}" for c in (
    "n CH4", "n CO2", "n H2O", "n C5–C15", "y CH4", "y CO2", "y H2O", "y C5–C15",
    "n H2", "n CO"))
BASE = "base H2/CO=2 (supuesto)"

titulo(f"Análisis 1: base (SUPUESTO b = C 1, H 4, O 1), sin φ — {T:g} K, {P:g} bar")
print(ENC)
for asf in (False, True):
    print(f"{'ASF ' + ('sí' if asf else 'no'):<40}" + fila((BASE, asf, False)))

titulo("Análisis 2: base, con y sin φ en la función objetivo")
print(ENC)
for asf in (False, True):
    for con_phi in (False, True):
        print(f"{f'ASF {SI_NO[bool(asf)]}, φ {SI_NO[bool(con_phi)]}':<40}"
              + fila((BASE, asf, con_phi)))

titulo("Análisis 3: biomasa (H2/CO = 1) y gas natural (H2/CO = 2.5)")
print(ENC)
for nombre in list(ALIMENTACIONES)[1:]:
    for asf in (False, True):
        for con_phi in (False, True):
            print(f"{f'{nombre}, ASF {SI_NO[bool(asf)]}, φ {SI_NO[bool(con_phi)]}':<40}"
                  + fila((nombre, asf, con_phi)))

# %% 2. % en peso (opcional de la diapositiva)
titulo("Opcional: % en peso entre parafinas, PM = 12k + (2k+2)")
print(f"{'':<46}{'C1':>8}{'C2–C4':>8}{'C5–C11':>9}{'C12–C15':>9}")
for clave, r in RES.items():
    nombre, asf, con_phi = clave
    w = pct_peso(r.n)
    print(f"{f'{nombre}, ASF {SI_NO[bool(asf)]}, φ {SI_NO[bool(con_phi)]}':<46}"
          + "".join(f"{v:>8.2f}" if g != "C5–C11" else f"{v:>9.2f}" for g, v in w.items()))
w_an = pct_peso_asf_analitico(ALFA)
print(f"{'ASF analítico (α = 0.80, cualquier b y φ)':<46}"
      + "".join(f"{v:>8.2f}" if g != "C5–C11" else f"{v:>9.2f}" for g, v in w_an.items()))

# %% 3. Pregunta A
r0, r1 = RES[(BASE, False, False)], RES[(BASE, True, False)]
c5_0 = sum(r0.n[sp] for sp in PARAFINAS if NCARB[sp] >= 5)
titulo("Pregunta A")
print(f"""Sin ASF el equilibrio es casi puro metano: n_C1 = {r0.n['CH4']:.4f} mol (y = {r0.y['CH4']:.3f},
{pct_peso(r0.n)['C1']:.2f} % en peso de las parafinas) y la suma C5–C15 es {c5_0:.2e} mol.
Cada CH2 extra cuesta +16.3 kJ/mol de ΔG a 500 K, así que la termodinámica empuja todo el
carbono a C1 (y a CO2 + H2O): las cadenas largas son muy inestables frente al metano.
Con ASF la distribución la fija la cinética de crecimiento de cadena (α), no el equilibrio:
C5–C15 = {sum(r1.n[sp] for sp in PARAFINAS if NCARB[sp] >= 5):.4f} mol y el % en peso es el de la ASF analítica.
Por eso un reactor de equilibrio puro no sirve para simular combustibles líquidos en FT:
predeciría metano, cuando el producto real es la distribución cinética (ASF) que el
catalizador produce lejos del equilibrio entre parafinas. Hay que imponer la ASF (u otro
modelo cinético) como restricción.""")

# %% 4. Validación y carbono
titulo("Validación (7 puntos + restricciones) y prueba de fase del grafito")
print("Con ASF el punto 4 (K) es n/a: el equilibrio está restringido.")
for clave, v in VALID.items():
    nombre, asf, con_phi = clave
    d_c = prueba_de_fase(RES[clave], "C(s)", 0.0)
    print(f"  {f'{nombre}, ASF {SI_NO[bool(asf)]}, φ {SI_NO[bool(con_phi)]}':<46}"
          f"{v.n_pass} PASS {v.n_warn} WARN {v.n_fail} FAIL;  Δ_C = {d_c:+.3f} (a_C = {np.exp(-d_c):.3g})")
print("El WARN es el punto 7: a_C > 1. A 500 K y 25 bar el grafito (2 CO -> C + CO2) bajaría G;")
print("la maestra no lo incluye, así que estos equilibrios son solo-gas, metaestables frente a C(s).")

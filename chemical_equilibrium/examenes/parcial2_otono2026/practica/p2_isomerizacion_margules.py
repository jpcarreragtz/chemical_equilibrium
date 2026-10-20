"""Ejercicios adicionales 2do parcial — PROBLEMA 2 (isomerización líquida)

    A(l) -> B(l),   A y B líquidos miscibles,   G^E/RT = 0.1·x_A·x_B
    ΔG°298 = −1000 J.  ¿Composición de equilibrio a 25 °C? ¿Qué error se
    comete si se supone solución ideal?

Planteamiento (a mano)
----------------------
Fase líquida, estado estándar = líquido puro a T y 1 bar  →  a_i = x_i·γ_i
(la presión no entra: f_i/f_i° ≈ 1 para líquidos)

    K = a_B/a_A = (x_B γ_B)/(x_A γ_A),     ln K = −ΔG°/RT = 1000/(8.314·298.15)
    Margules 1 parámetro (EQ. Físico):  G^E/RT = A·x_A·x_B,  A = 0.1
        ln γ_A = A·x_B²      ln γ_B = A·x_A²
    Ecuación a resolver (x_B = 1 − x_A):
        ln[(1 − x_A)/x_A] + A·(x_A² − x_B²) = ln K
        (x_A² − x_B² = x_A − x_B = 2x_A − 1)
        ⇒  ln[(1 − x_A)/x_A] + A(2x_A − 1) = ln K
    Solución ideal (γ = 1):  x_B/x_A = K  ⇒  x_A = 1/(1 + K)
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))   # raíz

from math import exp, log
from scipy.optimize import brentq

from core.nonideal import margules_GE_RT, margules_gamma
from core.thermo import R_J

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
T = 298.15            # K
DG0 = -1000.0         # J/mol
A12, A21 = 0.1, None  # Margules: 1 parámetro (A21 = None) o 2 parámetros
                      # (convención: 1 = A, 2 = B; ln γA∞ = A12, ln γB∞ = A21)
CLAVE = {"a": 0.3955, "b": 0.4005}
# ══════════════════════════ FIN DE DATOS ══════════════════════════════════

K = exp(-DG0 / (R_J * T))


def F(xA):
    gA, gB = margules_gamma(xA, A12, A21)
    return log((1 - xA) * gB / (xA * gA)) - log(K)


xA = brentq(F, 1e-12, 1 - 1e-12, xtol=1e-15)
gA, gB = margules_gamma(xA, A12, A21)
xA_id = 1.0 / (1.0 + K)
err = (xA_id - xA) / xA * 100

print("=" * 72)
print("  PROBLEMA 2 — A(l) -> B(l), Margules G^E/RT = 0.1 x_A x_B")
print("=" * 72)
print(f"\n1) ln K = −ΔG°/RT = {-DG0:.0f}/({R_J}·{T}) = {log(K):.5f}  → "
      f"K = {K:.5f}")
print(f"2) Margules: A12 = {A12}, A21 = {A21 if A21 is not None else A12}"
      f"  ({'1 parámetro' if A21 is None else '2 parámetros'})")
print("   ln γA = A·xB²,  ln γB = A·xA²")
print("\n3) NO IDEAL (actividades):")
print(f"   x_A = {xA:.4f}   x_B = {1 - xA:.4f}")
print(f"   γ_A = {gA:.5f}   γ_B = {gB:.5f}   G^E/RT = "
      f"{margules_GE_RT(xA, A12, A21):.5f}")
print(f"   comprobación: x_Bγ_B/(x_Aγ_A) = {(1 - xA) * gB / (xA * gA):.5f} = K")
print(f"   (maestra a) x_A = {CLAVE['a']})")
print("\n4) SOLUCIÓN IDEAL (γ = 1):  x_A = 1/(1 + K)")
print(f"   x_A = {xA_id:.4f}   x_B = {1 - xA_id:.4f}   (maestra b) "
      f"x_A = {CLAVE['b']})")
print(f"\n5) Error al suponer solución ideal: (x_A,id − x_A)/x_A = "
      f"{err:+.2f} %  (en x_B: {((1 - xA_id) - (1 - xA)) / (1 - xA) * 100:+.2f} %)")

print("\nAUTOVALIDACIÓN")
ok = True
for nom, v, cl in (("a", xA, CLAVE["a"]), ("b", xA_id, CLAVE["b"])):
    b = abs(v - cl) < 6e-5
    ok &= b
    print(f"  [{'PASS' if b else 'FAIL'}] {nom}) {v:.5f} vs clave {cl}")
b = abs(F(xA)) < 1e-12
ok &= b
print(f"  [{'PASS' if b else 'FAIL'}] residual ln K_calc − ln K = {F(xA):.1e}")
g1, g2 = margules_gamma(1.0, A12, A21)
b = abs(g1 - 1) < 1e-15 and abs(margules_gamma(0.0, A12, A21)[1] - 1) < 1e-15
ok &= b
print(f"  [{'PASS' if b else 'FAIL'}] límites γ_i → 1 cuando x_i → 1")
print(f"\nVEREDICTO: {'OK' if ok else 'REVISAR'}")

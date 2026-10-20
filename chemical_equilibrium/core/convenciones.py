"""
convenciones.py — ÚNICO lugar donde se decide cómo entra el término de
coeficientes de fugacidad K_phi en la constante de equilibrio en fase gas.

K_phi = Π phi_i^nu_i, con phi_i de especie PURA a (T, P) (regla de
Lewis-Randall, virial de Pitzer: ln phi = (Pr/Tr)(B0 + ω·B1)); se calcula
con core.nonideal.K_phi_lewis y los datos de data/criticas_b1.json.

Dos convenciones (CONVENCION_KPHI elige la que se usa por omisión):

  "maestra"  (la que pide la profesora del curso — DEFAULT)
      K = K_y · (P/P°)^δ / K_phi        ⇔   K_y = K · (P/P°)^(−δ) · K_phi
      Problema 1 d) del 2do parcial: T = 570.4 K

  "rigurosa" (de f^_i = phi_i·y_i·P, SVA ec. 13.28)
      K = K_phi · K_y · (P/P°)^δ        ⇔   K_y = K · (P/P°)^(−δ) / K_phi
      Problema 1 d) del 2do parcial: T = 584.5 K

En ambos casos se escribe  K = F · K_y · (P/P°)^δ  con el FACTOR
      F = factor_Kphi(K_phi, convencion)   ("rigurosa" → K_phi; "maestra" → 1/K_phi)
y el solver de gas ideal (que resuelve K_y·(P/P°)^δ = K) recibe la
"K efectiva"  K_ef = K_efectiva(K, K_phi, convencion) = K / F.

Para cambiar el default de TODO el proyecto edita solo CONVENCION_KPHI.
Cualquier script puede forzar la otra pasando convencion="..." .
"""

from __future__ import annotations

from typing import Optional

CONVENCIONES = ("maestra", "rigurosa")
CONVENCION_KPHI = "maestra"          # ← default del proyecto (examen práctico)


def validar_convencion(convencion: Optional[str] = None) -> str:
    """Devuelve la convención a usar (None → CONVENCION_KPHI); ValueError
    si no es 'maestra' ni 'rigurosa'."""
    c = CONVENCION_KPHI if convencion is None else convencion
    if c not in CONVENCIONES:
        raise ValueError(f"convención de K_phi desconocida: {c!r}; opciones: "
                         f"{CONVENCIONES} (ver core/convenciones.py)")
    return c


def otra_convencion(convencion: Optional[str] = None) -> str:
    """La convención contraria a la dada (o a la default)."""
    c = validar_convencion(convencion)
    return CONVENCIONES[1] if c == CONVENCIONES[0] else CONVENCIONES[0]


def factor_Kphi(Kphi: float, convencion: Optional[str] = None) -> float:
    """F tal que  K = F · K_y · (P/P°)^δ .
    'rigurosa' → F = K_phi ;  'maestra' → F = 1/K_phi ."""
    c = validar_convencion(convencion)
    if Kphi <= 0.0:
        raise ValueError(f"K_phi debe ser > 0 (recibido {Kphi})")
    return Kphi if c == "rigurosa" else 1.0 / Kphi


def K_efectiva(K: float, Kphi: float, convencion: Optional[str] = None) -> float:
    """K 'tipo gas ideal' que hay que pasarle al solver ideal:
        K_ef = K / F   →   el solver resuelve  K_y·(P/P°)^δ = K_ef
    y la solución ya incluye el efecto de la fugacidad."""
    return K / factor_Kphi(Kphi, convencion)


def describir(convencion: Optional[str] = None) -> str:
    """Texto corto de la fórmula de la convención (para imprimir)."""
    c = validar_convencion(convencion)
    if c == "maestra":
        return ("K = K_y·(P/P°)^δ / K_phi   ⇔   K_y = K·(P/P°)^(−δ)·K_phi"
                "   [convención de la maestra; F = 1/K_phi]")
    return ("K = K_phi·K_y·(P/P°)^δ   ⇔   K_y = K·(P/P°)^(−δ)/K_phi"
            "   [SVA ec. 13.28; F = K_phi]")

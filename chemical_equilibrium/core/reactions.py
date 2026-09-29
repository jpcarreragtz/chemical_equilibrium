"""
reactions.py — Parser de reacciones escritas como texto.

    "C3H8 + 3 H2O = 3 CO + 7 H2"   ->   {"C3H8": -1, "H2O": -3, "CO": 3, "H2": 7}

* Separadores aceptados: "=", "<=>", "<->", "⇌", "->", "→".
* Coeficientes: enteros ("3 H2O"), decimales ("0.5 O2") o fracciones
  ("1/2 O2"); sin coeficiente = 1. El espacio es opcional ("3H2O").
* El NOMBRE de la especie es su fórmula ("C3H8", "H2O", "Ca(OH)2"), de la
  que se deduce automáticamente la composición atómica — por eso NO uses
  nombres como "propano" o "n-C4H10".
* Cada reacción se verifica balanceada elemento por elemento al parsear;
  si no balancea se lanza ValueError diciendo reacción y elemento
  (escribe TÚ los coeficientes: el parser no balancea por ti).

build_system() arma todo lo que pide equilibrium.solve_extents:
especies (orden reproducible), matriz nu, dict de átomos e inertes.
"""

from __future__ import annotations

import re
from typing import Dict, List, Sequence, Tuple

_SEPARADORES = ("<=>", "<->", "⇌", "→", "->", "=")
_RE_COEF = re.compile(r"^\s*(\d+(?:\.\d+)?(?:\s*/\s*\d+(?:\.\d+)?)?)?\s*(.+?)\s*$")
_RE_ELEMENTO = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> Dict[str, int]:
    """Composición atómica de una fórmula química: 'C3H8' -> {'C':3,'H':8}.
    Soporta paréntesis: 'Ca(OH)2' -> {'Ca':1,'O':2,'H':2}."""
    def _rec(s: str, i: int) -> Tuple[Dict[str, int], int]:
        counts: Dict[str, int] = {}
        while i < len(s):
            c = s[i]
            if c == "(":
                sub, i = _rec(s, i + 1)
                m = re.match(r"\d+", s[i:])
                mult = int(m.group()) if m else 1
                i += m.end() if m else 0
                for el, n in sub.items():
                    counts[el] = counts.get(el, 0) + n * mult
            elif c == ")":
                return counts, i + 1
            else:
                m = _RE_ELEMENTO.match(s, i)
                if not m or m.start() != i or not m.group(1):
                    raise ValueError(f"Fórmula química inválida: '{formula}' "
                                     f"(problema en '...{s[i:]}')")
                el, num = m.group(1), m.group(2)
                counts[el] = counts.get(el, 0) + (int(num) if num else 1)
                i = m.end()
        return counts, i

    formula = formula.strip()
    if not formula:
        raise ValueError("Fórmula vacía")
    counts, _ = _rec(formula, 0)
    if not counts:
        raise ValueError(f"Fórmula química inválida: '{formula}'")
    return counts


def _parse_termino(term: str) -> Tuple[float, str]:
    """'3 H2O' -> (3.0, 'H2O');  '1/2 O2' -> (0.5, 'O2');  'CO' -> (1.0, 'CO')."""
    m = _RE_COEF.match(term)
    if not m or not m.group(2):
        raise ValueError(f"Término de reacción inválido: '{term}'")
    coef_txt, especie = m.group(1), m.group(2).strip()
    if coef_txt is None:
        coef = 1.0
    elif "/" in coef_txt:
        num, den = coef_txt.split("/")
        coef = float(num) / float(den)
    else:
        coef = float(coef_txt)
    if not re.match(r"^[A-Z(]", especie):
        raise ValueError(f"Término de reacción inválido: '{term}' "
                         f"(la especie debe empezar con mayúscula)")
    return coef, especie


def parse_reaction(reaccion: str) -> Dict[str, float]:
    """'C3H8 + 3 H2O = 3 CO + 7 H2' -> {especie: nu} (reactivos < 0)."""
    for sep in _SEPARADORES:
        if sep in reaccion:
            izq, der = reaccion.split(sep, 1)
            break
    else:
        raise ValueError(f"La reacción '{reaccion}' no tiene separador "
                         f"(usa '=', '<=>' o '->')")
    nu: Dict[str, float] = {}
    for lado, signo in ((izq, -1.0), (der, +1.0)):
        for term in lado.split("+"):
            if not term.strip():
                raise ValueError(f"Término vacío en '{reaccion}' "
                                 f"(¿un '+' de más?)")
            coef, sp = _parse_termino(term)
            nu[sp] = nu.get(sp, 0.0) + signo * coef
    nu = {sp: c for sp, c in nu.items() if abs(c) > 1e-12}
    if not any(c < 0 for c in nu.values()) or not any(c > 0 for c in nu.values()):
        raise ValueError(f"La reacción '{reaccion}' necesita al menos un "
                         f"reactivo y un producto")
    return nu


def verificar_balance(reaccion: str, nu: Dict[str, float],
                      atomos: Dict[str, Dict[str, float]] = None) -> None:
    """Lanza ValueError nombrando el elemento si la reacción no balancea.
    Con `atomos` explícito se usa esa composición (especies abstractas);
    si no, se deduce de la fórmula del nombre."""
    neto: Dict[str, float] = {}
    for sp, coef in nu.items():
        comp = atomos[sp] if atomos is not None else parse_formula(sp)
        for el, cnt in comp.items():
            neto[el] = neto.get(el, 0.0) + coef * cnt
    malos = {el: v for el, v in neto.items() if abs(v) > 1e-9}
    if malos:
        detalle = ", ".join(f"{el}: {v:+g}" for el, v in sorted(malos.items()))
        raise ValueError(f"Reacción NO balanceada: '{reaccion}' "
                         f"(átomos netos: {detalle}). Escribe los "
                         f"coeficientes estequiométricos correctos.")


def build_system(reacciones: Sequence[str], alimentacion: Dict[str, float],
                 atomos: Dict[str, Dict[str, float]] = None
                 ) -> Tuple[List[str], List[List[float]], Dict[str, Dict[str, int]],
                            List[str]]:
    """
    A partir de reacciones en texto y la alimentación construye:
        especies : orden reproducible (aparición en las reacciones, luego
                   las especies solo-alimentadas, p. ej. inertes)
        nu       : matriz (n_reacciones x n_especies) para solve_extents
        atomos   : {especie: {elemento: cuenta}} deducido de las fórmulas,
                   o el dict `atomos` dado (para especies ABSTRACTAS tipo
                   A, R, S: declara tú los fragmentos conservados, p. ej.
                   A = {"Rm": 1, "Sm": 2}, R = {"Rm": 1}, S = {"Sm": 1})
        inertes  : especies alimentadas que no aparecen en ninguna reacción
    """
    nu_dicts = [parse_reaction(rx) for rx in reacciones]

    especies: List[str] = []
    for nud in nu_dicts:
        for sp in nud:
            if sp not in especies:
                especies.append(sp)
    inertes = [sp for sp in alimentacion if sp not in especies]
    especies += inertes

    if atomos is not None:
        sin_comp = [sp for sp in especies if sp not in atomos]
        if sin_comp:
            raise KeyError("Diste `atomos` explícito pero faltan especies: "
                           + ", ".join(sin_comp))
    for rx, nud in zip(reacciones, nu_dicts):
        verificar_balance(rx, nud, atomos)

    nu = [[float(nud.get(sp, 0.0)) for sp in especies] for nud in nu_dicts]
    atomos_out = (dict(atomos) if atomos is not None
                  else {sp: parse_formula(sp) for sp in especies})
    return especies, nu, atomos_out, inertes

"""
independencia.py — Reacciones independientes por el método matricial
(Sesiones 23 y 24). Aritmética EXACTA con fractions.Fraction: nunca flotantes.

    A    : matriz atómica (especies × elementos), A[i][k] = átomos de k en i
    Aᵀ   : (elementos × especies); rango(Aᵀ) = rango(A)
    r    = N − rango                       reacciones independientes
    Aᵀ con columnas [componentes | no componentes]  --renglón-->  [I | Z]
    ν    : una reacción por especie NO componente j:
           ν_j = +1 (se forma j),  ν_componente = −Z[:, j],  resto 0;
           luego se escala a enteros mínimos.
    Aᵀν  = 0                               (conservación de cada elemento)
    F    = 2 − π + N − r − s               (regla de las fases)

Las matrices son listas de renglones (List[List[Fraction]]). ν se guarda
con UN RENGLÓN POR REACCIÓN y columnas en el orden de `especies`; la
comprobación Aᵀν = 0 usa ν con una columna por reacción (su transpuesta).

La composición atómica sale de core.reactions.parse_formula, y las
reacciones en texto se validan con core.reactions.parse_reaction y
verificar_balance (mismos separadores y mensajes de error que el solver).
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, lcm
from typing import Dict, List, Sequence, Tuple, Union

from core.reactions import (_RE_COEF, _SEPARADORES, parse_formula,
                            parse_reaction, verificar_balance)

Matriz = List[List[Fraction]]
Reaccion = Union[str, Dict[str, object]]


# ------------------------------------------------------------ utilidades
def transpuesta(M: Sequence[Sequence[Fraction]]) -> Matriz:
    return [list(col) for col in zip(*M)] if M else []


def _copia(M: Sequence[Sequence[object]]) -> Matriz:
    return [[Fraction(x) for x in fila] for fila in M]


def _fmt(x: Fraction) -> str:
    return str(x)       # Fraction: '3', '-4', '1/2'


def formatear_matriz(M: Sequence[Sequence[Fraction]],
                     filas: Sequence[str] = None,
                     columnas: Sequence[str] = None,
                     corte: int = None) -> str:
    """Matriz como texto alineado. `corte` = índice de columna donde se
    dibuja la barra vertical (p. ej. [I | Z])."""
    n_col = len(M[0]) if M else (len(columnas) if columnas else 0)
    celdas = [[_fmt(x) for x in fila] for fila in M]
    ancho = max([len(c) for fila in celdas for c in fila]
                + [len(c) for c in (columnas or [])] + [1])
    ancho_f = max([len(f) for f in (filas or [])] + [0])

    def linea(valores: Sequence[str], izq: str) -> str:
        partes = []
        for j, v in enumerate(valores):
            if corte is not None and j == corte:
                partes.append("|")
            partes.append(v.rjust(ancho))
        return f"{izq.ljust(ancho_f)}  " + " ".join(partes)

    out = []
    if columnas:
        out.append(linea(list(columnas), ""))
    for i, fila in enumerate(celdas):
        etiqueta = filas[i] if filas else ""
        out.append(linea(fila, etiqueta) if n_col else etiqueta)
    return "\n".join(out)


# ------------------------------------------------------- matriz atómica
def matriz_atomica(especies: Sequence[str],
                   elementos: Sequence[str] = None) -> Tuple[Matriz, List[str]]:
    """A (especies × elementos) a partir de las fórmulas. Sin `elementos`
    se usan en orden de primera aparición. Devuelve (A, elementos)."""
    comps = [parse_formula(sp) for sp in especies]
    if elementos is None:
        elementos = []
        for c in comps:
            for el in c:
                if el not in elementos:
                    elementos.append(el)
    else:
        elementos = list(elementos)
        faltan = sorted({el for c in comps for el in c} - set(elementos))
        if faltan:
            raise ValueError("Las especies contienen elementos que no "
                             "están en `elementos`: " + ", ".join(faltan))
    A = [[Fraction(c.get(el, 0)) for el in elementos] for c in comps]
    return A, elementos


def rango(M: Sequence[Sequence[object]]) -> int:
    """Rango por eliminación gaussiana exacta."""
    M = _copia(M)
    if not M:
        return 0
    n_f, n_c = len(M), len(M[0])
    r = 0
    for j in range(n_c):
        p = next((i for i in range(r, n_f) if M[i][j] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(r + 1, n_f):
            if M[i][j] != 0:
                f = M[i][j] / M[r][j]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
        if r == n_f:
            break
    return r


def num_reacciones_independientes(A: Sequence[Sequence[object]]) -> int:
    """r = N − rango(A), con N = número de especies (renglones de A)."""
    return len(A) - rango(A)


# ------------------------------------------- elección de componentes
def elegir_componentes(especies: Sequence[str],
                       A: Sequence[Sequence[object]],
                       alimentacion: Sequence[str],
                       complemento: Sequence[str] = None,
                       componentes_enunciado: Sequence[str] = None
                       ) -> Tuple[List[str], str]:
    """
    Elige las especies componentes (pivotes de [I | Z]) con el criterio del
    proyecto (CLAUDE.md):
      1. Si el enunciado las da, se usan esas.
      2. Si no, primero las especies de alimentación, en el orden dado.
      3. Una de alimentación linealmente dependiente de las anteriores no
         se usa como pivote, y se avisa.
      4. Si faltan pivotes para llegar al rango, se completa con
         `complemento`; sin esa lista, con las demás especies en el orden
         de `especies`, mostrando cuáles se agregaron.
    `A` es la matriz atómica (especies × elementos) alineada con `especies`;
    `alimentacion` puede ser una lista o un dict (se usan sus claves en
    orden). Independencia lineal con aritmética exacta (Fraction).

    Devuelve (componentes, texto que explica cómo se eligieron). Lanza
    ValueError si una especie no pertenece al sistema o si no se alcanza
    el rango.
    """
    especies = list(especies)
    if len(A) != len(especies):
        raise ValueError("A debe tener un renglón por especie")
    fila = {sp: [Fraction(x) for x in A[i]] for i, sp in enumerate(especies)}
    listas = {"alimentacion": list(alimentacion or []),
              "complemento": list(complemento or []),
              "componentes_enunciado": list(componentes_enunciado or [])}
    for nombre, lst in listas.items():
        ajenas = [sp for sp in lst if sp not in fila]
        if ajenas:
            raise ValueError(f"Especies en `{nombre}` que no pertenecen al "
                             f"sistema: {', '.join(ajenas)} (especies: "
                             f"{', '.join(especies)})")
    rg = rango(A)
    lineas = [f"rango(A) = {rg}: se necesitan {rg} especies componentes."]

    if componentes_enunciado:
        comp = listas["componentes_enunciado"]
        lineas.append("Criterio 1: el enunciado da las componentes: "
                      + ", ".join(comp)
                      + (" (tienen prioridad sobre la alimentación)."
                         if listas["alimentacion"] else "."))
        if len(set(comp)) != len(comp) or len(comp) != rg \
                or rango([fila[sp] for sp in comp]) != rg:
            raise ValueError(
                "Las componentes del enunciado no alcanzan el rango: se "
                f"necesitan {rg} especies linealmente independientes y se "
                f"dieron {', '.join(comp)} (rango "
                f"{rango([fila[sp] for sp in comp])}).")
        lineas.append("Componentes finales: " + ", ".join(comp))
        return comp, "\n".join(lineas)

    comp: List[str] = []

    def intentar(sp: str) -> bool:
        if sp in comp:
            return False
        if rango([fila[c] for c in comp] + [fila[sp]]) > len(comp):
            comp.append(sp)
            return True
        return False

    lineas.append("Criterio 2: especies de alimentación en el orden dado: "
                  + (", ".join(listas["alimentacion"]) or "(ninguna)") + ".")
    for sp in listas["alimentacion"]:
        if len(comp) == rg or not intentar(sp):
            previas = ", ".join(comp) if comp else "(ninguna)"
            lineas.append(f"  AVISO (criterio 3): {sp} es linealmente dependiente "
                          f"de {previas}; no se usa como pivote.")
        else:
            lineas.append(f"  {sp}: pivote (alimentación).")

    if len(comp) < rg:
        faltan = rg - len(comp)
        if complemento:
            fuente = listas["complemento"]
            lineas.append(f"Criterio 4: faltan {faltan} pivote(s); se completa "
                          "con el complemento explícito: "
                          + ", ".join(fuente) + ".")
        else:
            fuente = [sp for sp in especies if sp not in comp]
            lineas.append(f"Criterio 4: faltan {faltan} pivote(s) y no hay "
                          "complemento; se completa con las demás especies "
                          "en el orden dado.")
        agregadas = []
        for sp in fuente:
            if len(comp) == rg:
                break
            if intentar(sp):
                agregadas.append(sp)
            elif complemento and sp not in comp:
                lineas.append(f"  {sp} (complemento) es linealmente dependiente "
                              f"de {', '.join(comp)}; se omite.")
        lineas.append("  Agregadas: " + (", ".join(agregadas) or "(ninguna)") + ".")
        if len(comp) < rg:
            raise ValueError(
                f"No se alcanza el rango {rg}: con la alimentación y el "
                f"complemento solo hay {len(comp)} pivote(s) independientes "
                f"({', '.join(comp) or 'ninguno'}). Agrega especies al "
                "complemento.")

    lineas.append("Componentes finales: " + ", ".join(comp))
    return comp, "\n".join(lineas)


# ------------------------------------------- reducción de Aᵀ a [I | Z]
def reducir_por_componentes(At: Sequence[Sequence[object]],
                            especies: Sequence[str],
                            componentes: Sequence[str]) -> Dict[str, object]:
    """
    Reordena las columnas de Aᵀ como [componentes | no componentes] y la
    lleva por operaciones de renglón a [I | Z], pivoteando SOLO en las
    columnas de los componentes (prefiere pivotes ±1 para evitar
    fracciones). Lanza ValueError si los componentes no son linealmente
    independientes o si no son exactamente rango(Aᵀ) especies.

    Devuelve dict con:
        columnas      : orden de columnas [componentes..., no componentes...]
        no_componentes: especies no componentes (en el orden de `especies`)
        inicial       : Aᵀ con columnas reordenadas
        pasos         : [(descripción, matriz tras la operación), ...]
        reducida      : [I | Z] (renglones nulos eliminados)
        Z             : bloque Z (componentes × no componentes)
    En las descripciones f1, f2, ... son los renglones POR POSICIÓN.
    """
    especies = list(especies)
    componentes = list(componentes)
    desconocidas = [c for c in componentes if c not in especies]
    if desconocidas:
        raise ValueError("Componentes que no están en la lista de especies: "
                         + ", ".join(desconocidas))
    if len(set(componentes)) != len(componentes):
        raise ValueError("Hay componentes repetidos: " + ", ".join(componentes))
    rg = rango(At)
    k = len(componentes)
    if k != rg:
        raise ValueError(f"Se necesitan exactamente rango(Aᵀ) = {rg} especies "
                         f"componentes y se dieron {k}: "
                         + ", ".join(componentes))

    no_comp = [sp for sp in especies if sp not in componentes]
    columnas = componentes + no_comp
    idx = [especies.index(sp) for sp in columnas]
    M = [[Fraction(fila[j]) for j in idx] for fila in At]
    inicial = [list(f) for f in M]
    n_f = len(M)
    nombres = [f"f{i + 1}" for i in range(n_f)]
    pasos: List[Tuple[str, Matriz]] = []

    def anotar(desc: str) -> None:
        pasos.append((desc, [list(f) for f in M]))

    for j in range(k):
        cand = [i for i in range(j, n_f) if M[i][j] != 0]
        if not cand:
            raise ValueError(
                "Las especies componentes NO son linealmente independientes: "
                f"la columna de '{componentes[j]}' no admite pivote después de "
                f"reducir con {', '.join(componentes[:j]) or '(ninguna)'}. "
                "Elige otro conjunto de componentes.")
        p = next((i for i in cand if abs(M[i][j]) == 1), cand[0])
        if p != j:
            M[j], M[p] = M[p], M[j]
            anotar(f"{nombres[j]} ↔ {nombres[p]}")
        piv = M[j][j]
        if piv != 1:
            M[j] = [x / piv for x in M[j]]
            anotar(f"{nombres[j]} ← ({_fmt(1 / piv)})·{nombres[j]}")
        for i in range(n_f):
            if i != j and M[i][j] != 0:
                f = M[i][j]
                M[i] = [a - f * b for a, b in zip(M[i], M[j])]
                if f == 1:
                    txt = f"{nombres[i]} ← {nombres[i]} − {nombres[j]}"
                elif f == -1:
                    txt = f"{nombres[i]} ← {nombres[i]} + {nombres[j]}"
                elif f > 0:
                    txt = f"{nombres[i]} ← {nombres[i]} − {_fmt(f)}·{nombres[j]}"
                else:
                    txt = f"{nombres[i]} ← {nombres[i]} + {_fmt(-f)}·{nombres[j]}"
                anotar(txt)

    # con k = rango, los renglones restantes quedan nulos (elementos redundantes)
    reducida = M[:k]
    Z = [fila[k:] for fila in reducida]
    return {"columnas": columnas, "no_componentes": no_comp,
            "inicial": inicial, "pasos": pasos, "reducida": reducida, "Z": Z}


# -------------------------------------------------- matriz estequiométrica
def enteros_minimos(v: Sequence[object]) -> List[int]:
    """Escala un vector racional a los enteros mínimos (mismo signo)."""
    v = [Fraction(x) for x in v]
    m = lcm(*[x.denominator for x in v]) if v else 1
    ent = [int(x * m) for x in v]
    g = gcd(*ent) if any(ent) else 1
    return [x // g for x in ent]


def matriz_estequiometrica(red: Dict[str, object],
                           especies: Sequence[str]) -> List[List[int]]:
    """ν (reacciones × especies, columnas en el orden de `especies`), una
    reacción por especie no componente, con enteros mínimos y la especie
    no componente como producto (ν > 0)."""
    especies = list(especies)
    columnas, Z = red["columnas"], red["Z"]
    k = len(Z)
    nu = []
    for j, sp in enumerate(red["no_componentes"]):
        v = [Fraction(0)] * len(especies)
        v[especies.index(sp)] = Fraction(1)
        for i in range(k):
            v[especies.index(columnas[i])] = -Z[i][j]
        nu.append(enteros_minimos(v))
    return nu


def reaccion_texto(nu_fila: Sequence[object], especies: Sequence[str],
                   flecha: str = "->") -> str:
    """[-1, -1, 1, 0, 3] -> 'CH4 + H2O -> CO + 3 H2' (orden de `especies`)."""
    def termino(c: Fraction, sp: str) -> str:
        c = abs(Fraction(c))
        return sp if c == 1 else f"{_fmt(c)} {sp}"
    izq = [termino(c, sp) for c, sp in zip(nu_fila, especies) if c < 0]
    der = [termino(c, sp) for c, sp in zip(nu_fila, especies) if c > 0]
    return f"{' + '.join(izq)} {flecha} {' + '.join(der)}"


def comprobar_conservacion(At: Sequence[Sequence[object]],
                           nu: Sequence[Sequence[object]]) -> Tuple[Matriz, bool]:
    """Aᵀν con ν en columnas (una por reacción). Devuelve (producto, es_cero)."""
    prod = [[sum((Fraction(a) * Fraction(b) for a, b in zip(fila, rx)),
                 Fraction(0)) for rx in nu] for fila in At]
    return prod, all(x == 0 for fila in prod for x in fila)


# ------------------------------- ν leída directamente de [−I | −Z] (CLAUDE.md)
def negar_renglones(M: Sequence[Sequence[object]]) -> Matriz:
    """F → −F en todos los renglones: [I | Z] -> [−I | −Z]."""
    return [[-Fraction(x) for x in fila] for fila in M]


def nu_completa(red: Dict[str, object]) -> Tuple[Matriz, List[str]]:
    """
    ν completa SIN escalar a partir de la reducción [I | Z]:
    renglones = componentes y luego no componentes (orden de red["columnas"]),
    una COLUMNA por reacción. Columna j = (−Z[:, j] ; identidad e_j), es
    decir, la reacción que forma la no componente j. Devuelve (ν, renglones).
    """
    menos = negar_renglones(red["reducida"])
    k = len(menos)
    no_comp = red["no_componentes"]
    nu = [[menos[i][k + j] for j in range(len(no_comp))] for i in range(k)]
    nu += [[Fraction(int(i == j)) for j in range(len(no_comp))]
           for i in range(len(no_comp))]
    return nu, list(red["columnas"])


def escalar_columnas(nu: Sequence[Sequence[object]]
                     ) -> Tuple[List[List[int]], List[Fraction]]:
    """Escala cada columna al mínimo entero con un factor POSITIVO (los
    signos no cambian). Devuelve (ν entera, factor de cada columna)."""
    cols = transpuesta(_copia(nu))
    enteras, factores = [], []
    for col in cols:
        ent = enteros_minimos(col)
        i = next(i for i, x in enumerate(col) if x != 0)
        factores.append(Fraction(ent[i]) / col[i])
        enteras.append(ent)
    return transpuesta(enteras), factores


def reaccion_desde_columna(columna: Sequence[object], filas: Sequence[str],
                           orden: Sequence[str] = None,
                           flecha: str = "->") -> str:
    """Lee una reacción tal cual de su columna de ν: ν < 0 reactivo,
    ν > 0 producto, ν = 0 no participa. `orden` fija el orden de los
    términos dentro de cada lado (por defecto, el de `filas`)."""
    nu = dict(zip(filas, (Fraction(x) for x in columna)))
    orden = list(orden) if orden is not None else list(filas)

    def termino(sp: str) -> str:
        c = abs(nu[sp])
        return sp if c == 1 else f"{_fmt(c)} {sp}"
    reactivos = [termino(sp) for sp in orden if nu.get(sp, 0) < 0]
    productos = [termino(sp) for sp in orden if nu.get(sp, 0) > 0]
    return f"{' + '.join(reactivos)} {flecha} {' + '.join(productos)}"


def nu_por_especie(nu: Sequence[Sequence[object]], filas: Sequence[str],
                   especies: Sequence[str]) -> List[list]:
    """ν por columnas (renglones = `filas`) -> una lista por reacción con
    columnas en el orden de `especies` (formato de matriz_estequiometrica)."""
    idx = [list(filas).index(sp) for sp in especies]
    return [[col[i] for i in idx] for col in transpuesta(nu)]


def grados_libertad(pi: int, N: int, r: int, s: int = 0) -> int:
    """Regla de las fases para sistemas reactivos: F = 2 − π + N − r − s."""
    return 2 - pi + N - r - s


# ------------------------------------------- independencia de reacciones
def _reaccion_exacta(rx: Reaccion) -> Dict[str, Fraction]:
    """Reacción (texto o dict) -> {especie: ν exacto}. El texto se valida
    con parse_reaction/verificar_balance; los coeficientes se releen como
    Fraction para no pasar por flotantes ('1/3 O2' queda 1/3)."""
    if isinstance(rx, dict):
        return {sp: Fraction(c) for sp, c in rx.items() if Fraction(c) != 0}
    nu_float = parse_reaction(rx)          # valida forma y reactivo/producto
    verificar_balance(rx, nu_float)
    for sep in _SEPARADORES:
        if sep in rx:
            izq, der = rx.split(sep, 1)
            break
    nu: Dict[str, Fraction] = {}
    for lado, signo in ((izq, -1), (der, 1)):
        for term in lado.split("+"):
            m = _RE_COEF.match(term)
            coef_txt, sp = m.group(1), m.group(2).strip()
            c = Fraction(coef_txt.replace(" ", "")) if coef_txt else Fraction(1)
            nu[sp] = nu.get(sp, Fraction(0)) + signo * c
    return {sp: c for sp, c in nu.items() if c != 0}


def matriz_reacciones(reacciones: Sequence[Reaccion],
                      especies: Sequence[str] = None
                      ) -> Tuple[Matriz, List[str]]:
    """Matriz (reacciones × especies) exacta. Sin `especies` se usan en orden
    de aparición. Acepta texto ('CH4 + H2O -> CO + 3 H2') o dicts."""
    nus = [_reaccion_exacta(rx) for rx in reacciones]
    if especies is None:
        especies = []
        for d in nus:
            for sp in d:
                if sp not in especies:
                    especies.append(sp)
    else:
        especies = list(especies)
        extra = sorted({sp for d in nus for sp in d} - set(especies))
        if extra:
            raise ValueError("Especies en las reacciones que no están en "
                             "`especies`: " + ", ".join(extra))
    return [[d.get(sp, Fraction(0)) for sp in especies] for d in nus], especies


def rango_reacciones(reacciones: Sequence[Reaccion],
                     especies: Sequence[str] = None) -> int:
    """Rango del conjunto: son independientes si rango = número de reacciones."""
    M, _ = matriz_reacciones(reacciones, especies)
    return rango(M)


def son_independientes(reacciones: Sequence[Reaccion]) -> bool:
    return rango_reacciones(reacciones) == len(reacciones)


def generan_mismo_espacio(conj1: Sequence[Reaccion],
                          conj2: Sequence[Reaccion]) -> bool:
    """True si ambos conjuntos generan el mismo espacio de reacciones:
    rango(1) = rango(2) = rango(1 ∪ 2)."""
    r1, r2 = rango_reacciones(conj1), rango_reacciones(conj2)
    return r1 == r2 == rango_reacciones(list(conj1) + list(conj2))

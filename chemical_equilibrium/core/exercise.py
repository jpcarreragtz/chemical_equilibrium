"""
exercise.py — Motor "meto datos del enunciado y sale todo" para ejercicios
de equilibrio (Método A, multirreacción, barrido de temperaturas).

La plantilla (plantilla_ejercicio.py) solo llena UN bloque de inputs y llama:

    run_exercise(titulo=..., reacciones=[...], alimentacion={...}, fase="gas",
                 P=..., P0=..., unidad="bar", temperaturas=[...],
                 modo_K="directo" | "termo", K={...} | datos_termo={...},
                 K_es_Kc=False, limitante=None)

Hace, por cada T:
  parsear reacciones -> nu/átomos/inertes, obtener K (dadas, o de ΔG°f/ΔH°f/Cp
  vía ec. 13.18 SVA, o convertidas de Kc), resolver con solve_extents,
  imprimir ξ / n_i / y_i / X del reactivo limitante (+ ΔH°, ΔG°, K si es
  modo termo) y correr selfcheck.validate (20 arranques aleatorios,
  |K_calc-K|/K < 1e-4). Al final: chequeo de Le Chatelier del barrido,
  tabla comparativa, CSV y gráfica y_i vs T (PNG).
"""

from __future__ import annotations

import csv as _csv
import re
from typing import Dict, List, Optional, Sequence

import numpy as np

from . import SALIDAS
from .equilibrium import Kc_to_Ky, conversion, format_table, solve_extents
from .reactions import build_system, parse_reaction
from .selfcheck import Problem, validate, validate_sweep
from .thermo import cargar_datos, reaction_thermo

# Paleta categórica validada (orden FIJO, seguro para daltonismo en pares
# adyacentes sobre fondo claro); la 9a serie se agrupa en "otras".
_PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
           "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
_GRIS_OTRAS = "#8a8a85"
_TINTA = "#1a1a19"


def _slug(texto: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", texto).strip("_")
    return s or "ejercicio"


def _buscar_K(K: Dict, T: float) -> Optional[Sequence[float]]:
    """Busca K[T] tolerando claves int/float (700 vs 700.0)."""
    for clave in (T, float(T), int(T) if float(T).is_integer() else None):
        if clave is not None and clave in K:
            return K[clave]
    return None


def run_exercise(*, titulo: str = "ejercicio",
                 reacciones: Sequence[str],
                 alimentacion: Dict[str, float],
                 fase: str = "gas",
                 P: float = 1.0, P0: float = 1.0, unidad: str = "bar",
                 temperaturas: Sequence[float],
                 modo_K: str = "directo",
                 K: Optional[Dict] = None,
                 datos_termo: Optional[Dict[str, dict]] = None,
                 especies_termo: Optional[Sequence[str]] = None,
                 K_es_Kc: bool = False,
                 limitante: Optional[str] = None,
                 atomos: Optional[Dict[str, Dict[str, float]]] = None,
                 csv_out: Optional[str] = None,
                 png_out: Optional[str] = None,
                 n_starts_check: int = 20,
                 k_tol_check: float = 1e-4) -> dict:
    """Resuelve el ejercicio completo. Devuelve un dict con todo
    (resultados por T, rutas de CSV/PNG, reportes de selfcheck)."""

    # ---------- 0) armar el sistema a partir de los strings ----------
    if fase not in ("gas", "liquido", "líquido"):
        raise ValueError("fase debe ser 'gas' o 'liquido'")
    es_liquido = fase != "gas"
    # `atomos` explícito = especies abstractas (A, R, S...); si es None se
    # deducen de las fórmulas de los nombres.
    especies, nu, atomos, inertes = build_system(reacciones, alimentacion,
                                                 atomos)
    nu_dicts = [parse_reaction(rx) for rx in reacciones]
    nu_arr = np.array(nu, dtype=float)
    deltas = nu_arr.sum(axis=1)
    frac_sym = "x" if es_liquido else "y"

    # reactivo limitante para reportar X (por omisión: primer reactivo de R1
    # presente en la alimentación)
    if limitante is None:
        for sp in especies:
            if nu_dicts[0].get(sp, 0.0) < 0 and alimentacion.get(sp, 0.0) > 0:
                limitante = sp
                break
    if limitante is None or alimentacion.get(limitante, 0.0) <= 0:
        raise ValueError("No pude elegir reactivo limitante: pásalo con "
                         "limitante='...' (debe estar en la alimentación)")
    n_A0 = float(alimentacion[limitante])

    # ---------- datos termo: dict a mano, o cargados del JSON ----------
    if modo_K == "termo":
        if datos_termo is None:
            if especies_termo is None:
                raise ValueError(
                    "modo_K='termo' requiere datos_termo={...} o "
                    "especies_termo=[...] para cargar data/sva_tables.json")
            datos_termo = cargar_datos(especies_termo)
        elif especies_termo is not None:
            print("  (nota: se dieron datos_termo Y especies_termo; "
                  "se usan los datos_termo explícitos)")
        participantes = sorted({sp for nud in nu_dicts for sp in nud})
        sin_datos = [sp for sp in participantes if sp not in datos_termo]
        if sin_datos:
            raise ValueError(
                "Faltan datos termodinámicos para especies que participan "
                "en las reacciones: " + ", ".join(sin_datos)
                + ". Agrégalas a datos_termo/especies_termo (y al JSON si "
                  "no están: data/sva_tables.json).")

    # en fase líquida no hay factor de presión: K = Π x_i^ν (P/P0 = 1)
    P_ef, P0_ef = (1.0, 1.0) if es_liquido else (float(P), float(P0))
    if es_liquido and K_es_Kc and np.any(np.abs(deltas) > 1e-12):
        raise ValueError("Kc en fase líquida solo es igual a Kx si δ = 0 en "
                         "todas las reacciones; convierte tu Kc a Kx a mano")

    # ---------- encabezado: lo que se parseó (revísalo vs tu enunciado) ----
    print(f"\n{'=' * 70}\n  {titulo}\n{'=' * 70}")
    print(f"  Fase: {fase}   P = {P:g} {unidad}   P0 = {P0:g} {unidad}"
          + ("   (líquido: el factor de presión no aplica)" if es_liquido
             else f"   P/P0 = {P_ef / P0_ef:g}"))
    for j, rx in enumerate(reacciones):
        print(f"  R{j + 1}: {rx}    (delta_{j + 1} = {deltas[j]:+g})")
    print(f"  Especies (orden de columnas de nu): {', '.join(especies)}")
    print(f"  Inertes: {', '.join(inertes) if inertes else '(ninguno)'}   "
          f"Reactivo limitante para X: {limitante} (n0 = {n_A0:g} mol)")
    n0_tot = sum(float(v) for v in alimentacion.values())
    terms = "".join(f" {'+' if d > 0 else '-'} {abs(d):g}·xi_{j + 1}"
                    for j, d in enumerate(deltas) if abs(d) > 1e-12)
    print(f"  n_T = {n0_tot:g}{terms if terms else '   (constante)'}")
    if modo_K == "termo":
        pendientes = [sp for sp in sorted(datos_termo)
                      if "pendiente" in str(datos_termo[sp].get("verificado",
                                                                ""))]
        if pendientes:
            print("  [OJO] datos termo AÚN NO cotejados contra el libro "
                  "para: " + ", ".join(pendientes)
                  + "  (ver data/raw/LEEME.txt)")

    # ---------- 1) resolver cada temperatura ----------
    resultados: Dict[float, dict] = {}
    K_usadas: Dict[float, List[float]] = {}
    xi_por_T: Dict[float, tuple] = {}
    reportes = []

    for T in temperaturas:
        T = float(T)
        print(f"\n{'-' * 70}\n  T = {T:g} K\n{'-' * 70}")

        # --- K de esta temperatura ---
        termo_info = None
        if modo_K == "directo":
            if K is None:
                raise ValueError("modo_K='directo' requiere el dict K={T: [...]}")
            K_lista = _buscar_K(K, T)
            if K_lista is None:
                raise ValueError(f"No hay K para T = {T:g} en el dict K")
            K_lista = [float(k) for k in np.atleast_1d(K_lista)]
        elif modo_K == "termo":
            fuera = [f"{sp} (Tmax {datos_termo[sp]['Tmax']:g} K)"
                     for sp in participantes
                     if T > datos_termo[sp].get("Tmax", float("inf"))]
            if fuera:
                print("  [OJO] Cp extrapolado fuera de su rango de ajuste: "
                      + ", ".join(fuera))
            termo_info = [reaction_thermo(nud, datos_termo, T)
                          for nud in nu_dicts]
            K_lista = [t["K"] for t in termo_info]
        else:
            raise ValueError("modo_K debe ser 'directo' o 'termo'")

        if K_es_Kc and not es_liquido:
            # K_C = C_T0^δ · Π y^ν  con C_T0 = P0/(R·T) en mol/L
            K_lista = [Kc_to_Ky(kc, d, T0=T, P0=P0, unit=unidad)
                       for kc, d in zip(K_lista, deltas)]
            print("  Kc dadas convertidas a K (fracción mol): "
                  + ", ".join(f"K{j + 1} = {k:.5g}"
                              for j, k in enumerate(K_lista)))

        if termo_info is not None:
            print("  Termodinámica por reacción (ec. 13.18 SVA):")
            print("  " + format_table(
                [[f"R{j + 1}", t["dCp_R"], t["dH"] / 1000.0,
                  t["dG"] / 1000.0, t["K"]]
                 for j, t in enumerate(termo_info)],
                ["Rxn", "dCp/R(T)", "dH° (kJ/mol)", "dG° (kJ/mol)", "K"]
            ).replace("\n", "\n  "))
        else:
            print("  K usadas: " + ", ".join(f"K{j + 1} = {k:.5g}"
                                             for j, k in enumerate(K_lista)))

        # --- resolver ---
        res = solve_extents(especies, alimentacion, nu, K_lista,
                            P=P_ef, P0=P0_ef, atoms=atomos)
        X = conversion(n_A0, res.moles[limitante])

        print("  Avances:  " + "   ".join(
            f"xi_{j + 1} = {x:.5f} mol" for j, x in enumerate(res.xi)))
        print(f"  X_{limitante} = (n0 - n)/n0 = {X:.5f}      "
              f"n_T = {res.n_total:.5f} mol")
        print("  " + format_table(
            [[sp, res.moles[sp], res.mole_fractions[sp]] for sp in especies],
            ["Especie", "n_i (mol)", f"{frac_sym}_i"]).replace("\n", "\n  "))

        # --- selfcheck automático ---
        prob = Problem(species=especies, n0=alimentacion, nu=nu, K=K_lista,
                       P=P_ef, P0=P0_ef, atoms=atomos, inerts=inertes,
                       T=T, pressure_unit=unidad, label=titulo)
        rep = validate(prob, res, n_starts=n_starts_check, k_tol=k_tol_check)
        reportes.append(rep)

        resultados[T] = {"res": res, "K": K_lista, "X": X,
                         "termo": termo_info, "reporte": rep}
        K_usadas[T] = K_lista
        xi_por_T[T] = tuple(res.xi)

    # ---------- 2) Le Chatelier sobre el barrido ----------
    if len(resultados) >= 2:
        reportes.append(validate_sweep(K_usadas, xi_por_T))

    # ---------- 3) tabla comparativa ----------
    temps = sorted(resultados)
    print(f"\n{'=' * 70}\n  Resumen del barrido — {titulo}\n{'=' * 70}")
    filas = [[T] + list(resultados[T]["res"].xi) + [resultados[T]["X"]]
             + [resultados[T]["res"].mole_fractions[sp] for sp in especies]
             for T in temps]
    print("  " + format_table(
        filas, ["T (K)"] + [f"xi_{j + 1}" for j in range(len(reacciones))]
        + [f"X_{limitante}"] + [f"{frac_sym}_{sp}" for sp in especies]
    ).replace("\n", "\n  "))

    n_fail = sum(r.n_fail for r in reportes)
    n_warn = sum(r.n_warn for r in reportes)
    veredicto = "FAIL" if n_fail else ("WARN" if n_warn else "PASS")
    print(f"\n  Autovalidación global: {veredicto} "
          f"({n_fail} FAIL, {n_warn} WARN en {len(reportes)} reportes)")

    # ---------- 4) CSV (por omisión a salidas/, nunca a la raíz) ----------
    if csv_out is None:
        SALIDAS.mkdir(exist_ok=True)
    csv_path = csv_out or str(SALIDAS / f"resultados_{_slug(titulo)}.csv")
    _exportar_csv(csv_path, temps, resultados, especies, reacciones,
                  limitante, frac_sym)
    print(f"  CSV exportado: {csv_path}")

    # ---------- 5) gráfica y_i vs T ----------
    png_path = None
    if len(temps) >= 2:
        if png_out is None:
            SALIDAS.mkdir(exist_ok=True)
        png_path = png_out or str(SALIDAS / f"grafica_{_slug(titulo)}.png")
        err = _graficar(png_path, temps, resultados, especies, titulo,
                        frac_sym)
        print(f"  Gráfica {frac_sym}_i vs T: {png_path}" if err is None
              else f"  (sin gráfica: {err})")

    return {"titulo": titulo, "especies": especies, "nu": nu,
            "inertes": inertes, "limitante": limitante,
            "resultados": resultados, "csv": csv_path, "png": png_path,
            "veredicto": veredicto}


# ----------------------------------------------------------------------------
# CSV y gráfica
# ----------------------------------------------------------------------------

def _exportar_csv(path, temps, resultados, especies, reacciones,
                  limitante, frac_sym) -> None:
    n_rxn = len(reacciones)
    con_termo = resultados[temps[0]]["termo"] is not None
    cab = (["T_K"] + [f"xi_{j + 1}" for j in range(n_rxn)]
           + [f"X_{limitante}"]
           + [f"n_{sp}" for sp in especies]
           + [f"{frac_sym}_{sp}" for sp in especies]
           + [f"K_{j + 1}" for j in range(n_rxn)])
    if con_termo:
        cab += [f"dH_R{j + 1}_J_mol" for j in range(n_rxn)]
        cab += [f"dG_R{j + 1}_J_mol" for j in range(n_rxn)]
    with open(path, "w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(cab)
        for T in temps:
            r = resultados[T]
            fila = ([T] + [f"{x:.6g}" for x in r["res"].xi]
                    + [f"{r['X']:.6g}"]
                    + [f"{r['res'].moles[sp]:.6g}" for sp in especies]
                    + [f"{r['res'].mole_fractions[sp]:.6g}" for sp in especies]
                    + [f"{k:.6g}" for k in r["K"]])
            if con_termo:
                fila += [f"{t['dH']:.6g}" for t in r["termo"]]
                fila += [f"{t['dG']:.6g}" for t in r["termo"]]
            w.writerow(fila)


def _graficar(path, temps, resultados, especies, titulo, frac_sym):
    """Líneas y_i(T), paleta categórica fija, etiquetas directas al final.
    Devuelve None si se generó, o un mensaje de por qué no."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return "instala matplotlib (pip3 install matplotlib)"

    series = {sp: [resultados[T]["res"].mole_fractions[sp] for T in temps]
              for sp in especies}
    # máximo 8 series con color propio; el resto se suma en "otras"
    orden = sorted(especies, key=lambda s: -max(series[s]))
    principales, resto = orden[:8], orden[8:]
    principales = [sp for sp in especies if sp in principales]  # orden original

    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    for k, sp in enumerate(principales):
        ax.plot(temps, series[sp], color=_PALETA[k], lw=1.8,
                marker="o", ms=4.5, label=sp)
    if resto:
        otras = [sum(series[sp][i] for sp in resto) for i in range(len(temps))]
        ax.plot(temps, otras, color=_GRIS_OTRAS, lw=1.8, ls="--",
                marker="o", ms=4.5, label=f"otras ({len(resto)})")
        series["__otras__"] = otras

    # etiquetas directas al final de cada línea, en tinta (no en color)
    y_fin = sorted(((series[sp][-1], sp) for sp in principales), reverse=True)
    rango = max(0.04, max(v for v, _ in y_fin) - min(v for v, _ in y_fin))
    sep, y_prev = 0.035 * rango, None
    dx = 0.012 * (temps[-1] - temps[0])
    for v, sp in y_fin:
        y_lab = v if y_prev is None else min(v, y_prev - sep)
        ax.annotate(sp, (temps[-1], v), xytext=(temps[-1] + dx, y_lab),
                    fontsize=8.5, color=_TINTA, va="center")
        y_prev = y_lab

    ax.set_xlabel("T (K)", color=_TINTA)
    ax.set_ylabel(f"{frac_sym}_i (fracción mol)", color=_TINTA)
    ax.set_title(titulo, color=_TINTA, fontsize=11)
    ax.set_xlim(temps[0], temps[-1] + 6.5 * dx)   # aire para las etiquetas
    ax.grid(axis="y", color="#dddddd", lw=0.7)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color("#bbbbbb")
    ax.tick_params(colors="#555555", labelsize=8.5)
    ax.legend(loc="best", frameon=False, fontsize=8.5, labelcolor=_TINTA)
    fig.tight_layout()
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    return None

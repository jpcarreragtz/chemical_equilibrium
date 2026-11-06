"""Corre los 3 problemas del 2do parcial práctico (Grupo B) y resume veredictos.

    python3 correr_examen.py        (deja la salida completa en
                                     salida_examen.txt, junto a este script)
"""
import pathlib, subprocess, sys

AQUI = pathlib.Path(__file__).resolve().parent
SCRIPTS = ["p1_etilbenceno_estireno.py", "p2_combustible_CH4_H2.py",
           "p3_oxidacion_H2S.py"]

texto, resumen = [], []
for s in SCRIPTS:
    r = subprocess.run([sys.executable, "-W", "ignore", str(AQUI / s)],
                       capture_output=True, text=True, cwd=AQUI)
    out = r.stdout + (("\n[stderr]\n" + r.stderr) if r.returncode else "")
    texto.append(f"\n{'#' * 72}\n# {s}\n{'#' * 72}\n{out}")
    ver = next((l for l in out.splitlines() if l.startswith("VEREDICTO")),
               "VEREDICTO: (error al correr)")
    resumen.append(f"  {s:32s} {ver}")

(AQUI / "salida_examen.txt").write_text("".join(texto), encoding="utf-8")
print("\n".join(resumen))
print(f"\nSalida completa: {AQUI / 'salida_examen.txt'}")

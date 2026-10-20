"""Corre los 5 problemas de práctica del 2do parcial y resume veredictos.

    python3 correr_todo.py          (deja la salida completa en
                                     salida_practica.txt, junto a este script)
"""
import pathlib, subprocess, sys

AQUI = pathlib.Path(__file__).resolve().parent
SCRIPTS = ["p1_amoniaco_presion.py", "p2_isomerizacion_margules.py",
           "p3_gas_natural_hhv.py", "p4_reformado_calor.py",
           "p5_calentador_flujo_gas.py"]

texto, resumen = [], []
for s in SCRIPTS:
    r = subprocess.run([sys.executable, "-W", "ignore", str(AQUI / s)],
                       capture_output=True, text=True, cwd=AQUI)
    out = r.stdout + (("\n[stderr]\n" + r.stderr) if r.returncode else "")
    texto.append(f"\n{'#' * 72}\n# {s}\n{'#' * 72}\n{out}")
    ver = next((l for l in out.splitlines() if l.startswith("VEREDICTO")),
               "VEREDICTO: (error al correr)")
    resumen.append(f"  {s:32s} {ver}")

(AQUI / "salida_practica.txt").write_text("".join(texto), encoding="utf-8")
print("\n".join(resumen))
print(f"\nSalida completa: {AQUI / 'salida_practica.txt'}")

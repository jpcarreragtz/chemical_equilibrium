"""Corre el simulacro completo del 2do parcial y resume veredictos.

    python3 correr_simulacro.py     (salida completa en salida_simulacro.txt)
"""
import pathlib, subprocess, sys

AQUI = pathlib.Path(__file__).resolve().parent
CORRIDAS = [("a_metanol_presion.py", i) for i in "123"] + \
           [("b_flama_etileno.py", i) for i in "abcde"] + \
           [("c_amoniaco_NO_adiabatico.py", ""), ("d_syngas_dos_reacciones.py", "")] + \
           [("e_hhv_lhv_mezclas.py", i) for i in "12"]

texto, resumen, todo_ok = [], [], True
for script, arg in CORRIDAS:
    cmd = [sys.executable, str(AQUI / script)] + ([arg] if arg else [])
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=AQUI)
    out = r.stdout + (f"\n[stderr]\n{r.stderr}" if r.stderr else "")
    out = out.replace(str(AQUI.parents[2]), "<repo>")   # sin rutas locales en el txt
    texto.append(f"\n{'#' * 72}\n# {script} {arg}\n{'#' * 72}\n{out}")
    vers = [l for l in out.splitlines() if l.startswith("VEREDICTO")]
    # OK = sin error, sin [FAIL], todos los VEREDICTO en OK y el del simulacro presente
    ok = (r.returncode == 0 and "[FAIL]" not in out and all(l.endswith("OK") for l in vers)
          and any(l.startswith("VEREDICTO SIMULACRO") for l in vers))
    todo_ok &= ok
    res = next((l for l in out.splitlines() if l.startswith("RESULTADO")), "(sin RESULTADO)")
    avisos = sum(1 for l in r.stderr.splitlines() if "Warning" in l)
    fails = out.count("[FAIL]")
    resumen.append(f"  {'OK ' if ok else 'MAL'} {script:30s} {arg:2s} FAIL={fails} avisos={avisos}  {res}")

(AQUI / "salida_simulacro.txt").write_text("".join(texto), encoding="utf-8")
print("\n".join(resumen))
print(f"\nSIMULACRO: {'TODO OK' if todo_ok else 'REVISAR'}")
print(f"Salida completa: {AQUI / 'salida_simulacro.txt'}")
sys.exit(0 if todo_ok else 1)

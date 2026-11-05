"""Sesiones 23–24 — Reacciones independientes por el método matricial.

Reformado de metano con vapor. Especies (en este orden): CH4, H2O, CO, CO2,
H2; elementos C, H, O; una fase gaseosa (π = 1); alimentación 1 mol CH4 +
4 mol H2O; restricciones estequiométricas s = 2.

Se obtiene el conjunto de reacciones independientes (NO se calculan
composiciones de equilibrio):
  1. A (especies × elementos) y Aᵀ
  2. rango y r = N − rango
  3. Aᵀ con columnas [componentes | no componentes] -> [I | Z] por
     operaciones de renglón (componentes: CH4, H2O, H2)
  4. F → −F: [−I | −Z]; −Z da los ν de las componentes
  5. ν completa: −Z sobre la identidad, una columna por reacción, cada
     columna al mínimo entero positivo; reacciones leídas de sus signos
  6. Aᵀν = 0 y regla de las fases F = 2 − π + N − r − s

Aritmética exacta (fractions) en core/independencia.py.
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # raíz del proyecto

from core.independencia import (comprobar_conservacion, elegir_componentes,
                                escalar_columnas, formatear_matriz,
                                grados_libertad, matriz_atomica,
                                negar_renglones, nu_completa,
                                num_reacciones_independientes, rango,
                                reaccion_desde_columna, reducir_por_componentes,
                                transpuesta)

# ══════════════════════════ DATOS (edítalos aquí) ═════════════════════════
ESPECIES = ["CH4", "H2O", "CO", "CO2", "H2"]
ELEMENTOS = ["C", "H", "O"]
COMPLEMENTO = ["H2"]                     # completa los pivotes tras la alimentación
COMPONENTES_ENUNCIADO = None             # el enunciado no las da
PI = 1                                   # una fase gaseosa
ALIMENTACION = {"CH4": 1.0, "H2O": 4.0}  # mol; su orden da los primeros pivotes
S = 2                                    # restricciones estequiométricas


def titulo(t: str) -> None:
    print(f"\n{'=' * 64}\n{t}\n{'=' * 64}")


# %% 1. Matriz atómica A (especies × elementos) y su transpuesta
A, elementos = matriz_atomica(ESPECIES, ELEMENTOS)
At = transpuesta(A)
N = len(ESPECIES)

titulo("1. Matriz atómica A (especies × elementos)")
print(formatear_matriz(A, filas=ESPECIES, columnas=elementos))
print("\nAᵀ (elementos × especies)")
print(formatear_matriz(At, filas=elementos, columnas=ESPECIES))

# %% 2. Rango y número de reacciones independientes
rg = rango(At)
r = num_reacciones_independientes(A)

titulo("2. Rango y reacciones independientes")
print(f"N = {N} especies, rango(A) = rango(Aᵀ) = {rg}")
print(f"r = N − rango = {N} − {rg} = {r}")

# %% 3. Elección de componentes y reducción de Aᵀ a [I | Z]
COMPONENTES, CRITERIO = elegir_componentes(ESPECIES, A, ALIMENTACION,
                                           COMPLEMENTO, COMPONENTES_ENUNCIADO)
red = reducir_por_componentes(At, ESPECIES, COMPONENTES)
k = len(COMPONENTES)

titulo("3. Componentes y operaciones de renglón sobre Aᵀ")
print(CRITERIO)
print(f"\nComponentes (pivotes): {', '.join(COMPONENTES)}; "
      f"no componentes: {', '.join(red['no_componentes'])}")
print("\nAᵀ con columnas reordenadas [componentes | no componentes]:")
print(formatear_matriz(red["inicial"], filas=[f"f{i + 1}" for i in range(len(At))],
                       columnas=red["columnas"], corte=k))
for n, (desc, M) in enumerate(red["pasos"], 1):
    print(f"\nPaso {n}: {desc}")
    print(formatear_matriz(M, filas=[f"f{i + 1}" for i in range(len(M))],
                           columnas=red["columnas"], corte=k))

print("\nMatriz reducida [I | Z]:")
print(formatear_matriz(red["reducida"], filas=COMPONENTES,
                       columnas=red["columnas"], corte=k))

# %% 4. F → −F: [−I | −Z]
menos_IZ = negar_renglones(red["reducida"])

titulo("4. F → −F en todos los renglones: [−I | −Z]")
print("−Z contiene directamente los coeficientes ν de las componentes.")
print(formatear_matriz(menos_IZ, filas=COMPONENTES, columnas=red["columnas"],
                       corte=k))

# %% 5. Matriz ν completa (una columna por reacción)
nu_frac, FILAS_NU = nu_completa(red)
NU, factores = escalar_columnas(nu_frac)
nombres_rx = [f"R{i + 1}" for i in range(len(NU[0]))]

titulo("5. Matriz ν completa: −Z sobre la identidad de las no componentes")
if any(f != 1 for f in factores):
    print(formatear_matriz(nu_frac, filas=FILAS_NU, columnas=nombres_rx))
    print("\nCada columna se escala al mínimo entero positivo: "
          + ", ".join(f"{n} ×{f}" for n, f in zip(nombres_rx, factores)) + "\n")
print(formatear_matriz(NU, filas=FILAS_NU, columnas=nombres_rx))

# %% 6. Reacciones leídas de ν
titulo("6. Reacciones leídas de ν (ν < 0 reactivo, ν > 0 producto)")
for nombre, col in zip(nombres_rx, transpuesta(NU)):
    print(f"{nombre}: {reaccion_desde_columna(col, FILAS_NU, ESPECIES)}")

# %% 7. Comprobación Aᵀν = 0
prod, ok = comprobar_conservacion(red["inicial"], transpuesta(NU))

titulo("7. Comprobación Aᵀν = 0 (columnas de Aᵀ en el orden de los renglones de ν)")
print(formatear_matriz(prod, filas=elementos, columnas=nombres_rx))
print(f"\nAᵀν = 0: {'SÍ' if ok else 'NO'} — cada reacción conserva C, H y O")

# %% 8. Regla de las fases
F = grados_libertad(PI, N, r, S)

titulo("8. Regla de las fases")
print("Alimentación: " + ", ".join(f"{v:g} mol {sp}" for sp, v in ALIMENTACION.items()))
print(f"F = 2 − π + N − r − s = 2 − {PI} + {N} − {r} − {S} = {F}")
print("Con T y P fijas, el estado de equilibrio queda determinado.")

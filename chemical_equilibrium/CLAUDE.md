# chemical_equilibrium

## Estructura de ejercicios

Cada ejercicio nuevo va en su propio folder, `ejercicios/sNN_ejN_tema/`, con su `.py` (copiado de `plantilla/plantilla_ejercicio.py`) y su `.gms` (copiado de `plantilla/plantilla_gibbs.gms`) juntos. Al crearlo se agrega una fila al índice `ejercicios/README.md`: folder, sesión, tema, método y resultado clave, tomado de la salida real.

Los scripts no usan `parents[N]` para llegar a la raíz. Buscan hacia arriba la carpeta que contiene `core/`, e `import core.entorno` va justo después:

    RAIZ = next(p for p in pathlib.Path(__file__).resolve().parents if (p / "core").is_dir())
    sys.path.insert(0, str(RAIZ))
    import core.entorno  # noqa: F401

`tests/tests.py` carga los ejercicios por nombre de archivo con `cargar_ejercicio`, que los busca en `ejercicios/*/`. Los archivos que genera GAMS quedan dentro del folder del ejercicio y no se versionan.

## Reacciones independientes: elección de especies componentes

En los ejercicios de reacciones independientes, las especies componentes se eligen así:
1. Si el enunciado da las componentes, se usan esas.
2. Si no, se usan primero las especies de alimentación, en el orden dado.
3. Una especie de alimentación linealmente dependiente de las anteriores no se usa como pivote, y se avisa en la salida.
4. Si faltan pivotes para llegar al rango, se completa con la lista explícita complemento. Sin esa lista, se completa con las demás especies en el orden dado, y se muestra cuáles se agregaron.
5. El ejercicio siempre muestra el criterio y las componentes finales.

Implementación: `core.independencia.elegir_componentes(especies, A, alimentacion, complemento=None, componentes_enunciado=None)` devuelve `(componentes, texto)`; imprime el texto en el ejercicio. En `plantilla/plantilla_ejercicio.py` se activa con `reacciones_independientes = True` y las variables `alimentacion`, `complemento` y `componentes_enunciado`.

## Reacciones independientes: lectura de ν desde la matriz

En los ejercicios de reacciones independientes, después de reducir Aᵀ a [I | Z], se aplica F → −F a todos los renglones para obtener [−I | −Z]. −Z contiene directamente los coeficientes ν de las componentes. Se agrega una identidad para las no componentes, de modo que cada columna de ν sea una reacción que forma esa no componente. Si hay fracciones, cada columna se escala al mínimo entero positivo. Las reacciones se leen directamente de esos signos y no se muestra un paso de cambio de signo.

Orden de la salida: Aᵀ → operaciones de renglón hasta [I | Z] → F → −F con [−I | −Z] → ν completa (renglones: componentes y luego no componentes; R1, R2, ... como columnas) → reacciones leídas de ν (ν < 0 reactivo, ν > 0 producto, ν = 0 no participa) → Aᵀν = 0.

Implementación en `core/independencia.py`: `negar_renglones` (F → −F), `nu_completa` (ν sin escalar y sus renglones), `escalar_columnas` (mínimo entero positivo por columna, exacto), `reaccion_desde_columna` (lee la reacción de los signos; `orden` fija el orden de los términos, en los ejercicios el de la lista de especies) y `nu_por_especie` (ν en el formato de `matriz_estequiometrica` para K y solvers).

## Datos de las diapositivas

Cuando una diapositiva o un enunciado da un link de datos, siempre se revisa y se usan esos valores. Si un enunciado no trae un dato que debería, no se llena de memoria: se busca en ese link o en la fuente verificada del proyecto, y se dice de dónde salió.

Cuando la maestra da datos (ΔG, φ, K, tablas), se usan tal cual y tienen prioridad sobre JANAF, SVA o cualquier otra fuente. No se reemplazan ni se marcan como inconsistentes en los entregables. Si algo parece raro, se le pregunta al usuario. JANAF y SVA solo se usan cuando la maestra no da el dato.

## JANAF en cada ejercicio

Si un ejercicio necesita ΔG°f y JANAF aplica, los datos se consultan en el momento y nunca se llenan de memoria:

1. Buscar la especie en `https://janaf.nist.gov/tables/<Elemento>-index.html` (por ejemplo `C-index.html`) y tomar el ID de la tabla en la fase correcta (por ejemplo `C-067`).
2. Leer `https://janaf.nist.gov/tables/<ID>.txt`. Es texto separado por tabuladores: la columna `delta-f G` está en kJ/mol y T va cada 100 K.
3. Agregar el valor a `data/janaf_dgf.json` en J/mol: `"<especie>": {"tabla": "<ID>", "<T>": <valor>}`. Si la especie ya existe, solo se agrega la T.
4. Si T no está en la malla de 100 K, avisar antes de interpolar.
5. Si la especie no aparece en el índice (C2H6, C3H8), agregarla a `no_en_janaf` con la razón, usar `dGf_sva` y decirlo en el resultado.

## Minimización de la energía de Gibbs

Implementación en `core/gibbs.py`: `minimizar_gibbs` (RAND en ln n_i, sin valor inicial), `validar_gibbs` (7 puntos, más el 8 y el 9 cuando aplican), `reacciones_independientes`, `prueba_de_fase` y `cargar_dgf`.

ΔG°f: `cargar_dgf` usa `data/janaf_dgf.json`, capturado de https://janaf.nist.gov/. Una especie o T nueva se agrega ahí con su tabla. La ec. 13.18 de SVA (`dGf_sva`) solo se usa para especies declaradas en `no_en_janaf`, como C2H6.

Alcance: una fase gaseosa ideal con φ_i constante, más restricciones lineales extra y fases condensadas puras. No aplica con φ_i dependiente de y, ni con soluciones líquidas.

Restricciones extra (S27): `restricciones=[({especie: coef}, rhs), ...]` se agregan como filas a A n = b dentro de RAND, cada una con su multiplicador (`res.lam_restr`). Ejemplo, ASF: n_Ck − α^(k−1)·n_C1 = 0. Con restricciones, `validar_gibbs` incluye sus multiplicadores en el punto 3 (Lagrange), marca el punto 4 (K) como n/a porque el equilibrio está restringido, y verifica cada restricción en el punto 8. `validar_gibbs(..., nota_carbono="el modelo de la maestra no incluye C(s)")` deja el punto 7 (a_C) como n/a con esa nota (S27); sin el parámetro se comporta igual que siempre.

Fases condensadas puras y prueba de fase (S28): `condensadas=["C(s)", "Si(s)", "S(ref)", "H2O(l)", ...]`, con la fase en el nombre y su ΔG°f en `dGf` (las de referencia valen 0: JANAF C-002, S-001, Si-001). El procedimiento: (1) se resuelve solo el gas; (2) con sus λ, Δ_k = g_k/RT + Σ_e A_ek λ_e/RT, con la convención del Ejemplo 13.14 (λ_C/RT = 0.7635); (3) si Δ_k < 0, la fase entra con g_k/RT + Σ A λ = 0 y n_k ≥ 0, y se recalcula; (4) si n_k sale negativo, se quita. Sirve para fases de uno o varios elementos, siempre que sus elementos estén en alguna especie gaseosa. `res.delta_fase_gas` y `res.delta_fase` guardan Δ_k (solo gas y final), `res.n_cond` los moles presentes, y el punto 9 de `validar_gibbs` reporta todo. `prueba_de_fase(res, especie, dGf)` devuelve Δ_k sin recalcular. Para el grafito, Δ_C = λ_C/RT = −ln a_C.

Límite conocido del Newton RAND: desde algunos arranques al azar se estanca si una especie mayoritaria arranca como traza (por ejemplo el O2 en H2S:O2 = 1:2). Desde el reparto uniforme converge. `validar_gibbs` lo reporta como WARN en el punto 5, no como FAIL. Tampoco converge si dos elementos solo se distinguen por trazas menores a ~1e-16 (p. ej. H y O cuando solo hay H2O y trazas de H2 y O2).

GAMS: minimizar no basta para las trazas, porque casi no pesan en G y el solver se detiene aunque reporte "Locally Optimal". Se usa ln n_i como variable y una segunda etapa CNS con las condiciones de Lagrange, como en `ejercicios/s25_ej2_steam_cracker/s25_ej2_gibbs.gms`. No se deja el límite inferior en 0, porque ln(0) detiene el solve. Por la misma razón, una plantilla con y.lo = 1e-8 (como la de la maestra) deja en 1e-8 toda traza cuyo valor real sea menor: el solver para en la cota. Se usa `lnn.lo = -700`. Los `.gms` de S26–S28 escriben un CSV solo si se pasa `--salida=archivo.csv`; en el `put` va `.nz = 0`, porque si no, los valores menores a 1e-5 salen como 0.

Salida de los `.gms` (regla para TODO `.gms`, nuevo o existente; la usan S25–S28): todo ejercicio nuevo en GAMS se arma copiando `plantilla/plantilla_gibbs.gms` y debe dar resultados limpios como los `.py`. Los `.gms` empiezan con `$offListing` y `$offSymXRef offSymList`, y usan `option limrow = 0, limcol = 0, solprint = silent` (con `off` se queda el resumen de cada solve). IPOPT se calla con `ipopt.opt` (`print_level 0`, `sb yes`, `print_eval_error no`), que el mismo `.gms` escribe con `$onEcho` y activa con `<modelo>.optfile = 1`. Al final escriben `<nombre>_resultados.txt` junto al `.gms`, en el folder del ejercicio, con `put` (títulos del `.py`, ceros explícitos y una línea de estado con los solves que salieron bien) y lo copian a la terminal con `execute 'cat ...'`. Ese bloque va dentro de `$ifThen "%system.incParent%" == ""`, para que la revisión automática de la extensión de VS Code, que incluye el archivo desde otro `.gms`, no lo escriba. `ipopt.opt`, los `.lst` y `*_resultados.txt` están en `.gitignore` para cualquier subfolder de `ejercicios/`. Las líneas `---` de GAMS por cada solve no se pueden quitar desde el `.gms`; solo con `lo=0`, y la extensión de VS Code siempre pasa `lo=3`. Antes de entregar un `.gms` nuevo se comprueba que el `.lst` no pase de unos cientos de líneas, que `<nombre>_resultados.txt` tenga la tabla con los mismos títulos que el `.py` y que la línea de estado diga cuántos solves salieron bien.

Licencias: los `.gms` corren con la licencia Free Personal, vigente hasta el 30-nov-2027. Se usa IPOPT en NLP y CNS porque corre sin restricción; CONOPT y otros solvers comerciales corren ahí en modo demo. CONOPT sin modo demo requiere GAMSPy con la licencia Free Academic GAMSPy, que no ejecuta `.gms`.

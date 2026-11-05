# chemical_equilibrium

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

## JANAF en cada ejercicio

Si un ejercicio necesita ΔG°f y JANAF aplica, los datos se consultan en el momento y nunca se llenan de memoria:

1. Buscar la especie en `https://janaf.nist.gov/tables/<Elemento>-index.html` (por ejemplo `C-index.html`) y tomar el ID de la tabla en la fase correcta (por ejemplo `C-067`).
2. Leer `https://janaf.nist.gov/tables/<ID>.txt`. Es texto separado por tabuladores: la columna `delta-f G` está en kJ/mol y T va cada 100 K.
3. Agregar el valor a `data/janaf_dgf.json` en J/mol: `"<especie>": {"tabla": "<ID>", "<T>": <valor>}`. Si la especie ya existe, solo se agrega la T.
4. Si T no está en la malla de 100 K, avisar antes de interpolar.
5. Si la especie no aparece en el índice (C2H6, C3H8), agregarla a `no_en_janaf` con la razón, usar `dGf_sva` y decirlo en el resultado.

## Minimización de la energía de Gibbs

Implementación en `core/gibbs.py`: `minimizar_gibbs` (RAND en ln n_i, sin valor inicial), `validar_gibbs` (7 puntos), `reacciones_independientes` y `cargar_dgf`.

ΔG°f: `cargar_dgf` usa `data/janaf_dgf.json`, capturado de https://janaf.nist.gov/. Una especie o T nueva se agrega ahí con su tabla. La ec. 13.18 de SVA (`dGf_sva`) solo se usa para especies declaradas en `no_en_janaf`, como C2H6.

Alcance: una fase gaseosa con φ_i constante. No aplica con φ_i dependiente de y, ni con sólidos o líquidos.

GAMS: minimizar no basta para las trazas, porque casi no pesan en G y el solver se detiene aunque reporte "Locally Optimal". Se usa ln n_i como variable y una segunda etapa CNS con las condiciones de Lagrange, como en `ejercicios/s25_ej2_gibbs.gms`. No se deja el límite inferior en 0, porque ln(0) detiene el solve.

Licencias: los `.gms` corren con la licencia Free Personal, vigente hasta el 30-nov-2027. Se usa IPOPT en NLP y CNS porque corre sin restricción; CONOPT y otros solvers comerciales corren ahí en modo demo. CONOPT sin modo demo requiere GAMSPy con la licencia Free Academic GAMSPy, que no ejecuta `.gms`.

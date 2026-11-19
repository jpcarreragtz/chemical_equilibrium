# Equilibrio Químico — IBERO CDMX (Otoño 2026)

Trabajo del curso de **Equilibrio Químico** (Ingeniería Química, Universidad
Iberoamericana, Grupo B, Otoño 2026): un solver propio en Python con
autovalidación, modelos en GAMS para la minimización de Gibbs, las tablas de
Smith–Van Ness transcritas y verificadas, y los ejercicios y parciales
resueltos con ellos.

## Por dónde entrar

| carpeta | qué hay |
|---|---|
| [chemical_equilibrium/](chemical_equilibrium/) | El proyecto: solver (`core/`), datos, tests y su [README](chemical_equilibrium/README.md) con el mapa completo y cómo correr todo |
| [chemical_equilibrium/ejercicios/](chemical_equilibrium/ejercicios/README.md) | Índice de ejercicios: un folder por ejercicio, con sesión, tema, método y resultado clave |
| [chemical_equilibrium/examenes/](chemical_equilibrium/examenes/) | Parciales resueltos, prácticas y simulacros |
| [chemical_equilibrium/plantilla/](chemical_equilibrium/plantilla/) | Plantillas para empezar un ejercicio nuevo en Python o en GAMS |
| [Tablas y Valores/](Tablas%20y%20Valores/) | Tablas de Smith–Van Ness (Apéndices B y C) en Excel y JSON, en [tablas_sva/](Tablas%20y%20Valores/tablas_sva/) |
| [excel/](excel/) | Libros de Excel de ejercicios resueltos |

## Para empezar

```bash
pip3 install -r requirements.txt
cd chemical_equilibrium
python3 tests/tests.py      # exit 0 = todos los casos PASS
```

Los PDF del libro y las fotos de tablas no están en el repositorio por
derechos de autor (ver `.gitignore`).

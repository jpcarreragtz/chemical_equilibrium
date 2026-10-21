# Simulacro del 2do parcial (Otoño 2026)

Un problema nuevo por tipo de examen. Cada script es una **copia literal de una
plantilla** (mismo cuerpo, solo cambia el bloque DATOS/INPUTS) con el
planteamiento a mano en el docstring y una autovalidación independiente al final.

```bash
cd examenes/parcial2_otono2026/simulacro
python3 correr_simulacro.py                 # todo + resumen; salida en salida_simulacro.txt
python3 b_flama_etileno.py c                # o un inciso suelto
```

| Script | Plantilla | Incisos | Resultado |
|---|---|---|---|
| `a_metanol_presion.py` | `plantilla_presion_fugacidad.py` | 1, 2, 3 | y = 0.9291 (300 K, 1 bar); T = 364.5 K (1 bar); T = 504.2 K (100 bar, Lewis-Randall "maestra"; gas ideal 516.5 K) |
| `b_flama_etileno.py` | `plantilla_balance_energia.py` ("T_ad") | a–e | 2533.5, 2198.6, 1950.9, 1609.2, 2282.5 K (= Excel s18) |
| `c_amoniaco_NO_adiabatico.py` | `plantilla_balance_energia.py` ("T_ad") | — | T salida = 978.8 K |
| `d_syngas_dos_reacciones.py` | `plantilla_balance_energia.py` ("Q") | — | Q = 328.0 kJ/mol CH4 |
| `e_hhv_lhv_mezclas.py` | `plantilla_balance_energia.py` ("HHV") | 1, 2 | S17: HHV 932.95 / LHV 843.16 kJ/mol (LHV 46 830 kJ/kg); S16: HHV 858.36 / LHV 773.85 kJ/mol (53.24 / 48.00 MJ/kg), O2 191 mol, aire 909.5 mol por 100 mol |

Avisos esperados: Cp extrapolado (T > 2000 K) en b) a, b, e; en a) inciso 3,
φ de CO y CH3OH con Vr < 2 (virial de Pitzer extrapolada, Fig. 3.14 SVA).

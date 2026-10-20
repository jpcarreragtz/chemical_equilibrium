# 2do parcial (Otoño 2026) — ejercicios adicionales de práctica

Enunciados: slides 216–217 de `Notas/EQ. Químico.pdf`.
Planteamientos a mano y comparación con la clave: **`GUIA_planteamientos.md`**.

```bash
cd examenes/parcial2_otono2026/practica
python3 correr_todo.py            # los 5 problemas + resumen de veredictos
python3 p1_amoniaco_presion.py    # o uno por uno
```

| Script | Tema | Motor usado |
|---|---|---|
| `p1_amoniaco_presion.py` | K(T) ec. 13.18, efecto de P, Lewis-Randall (φ Pitzer) | `core.thermo`, `core.nonideal.K_phi_lewis` |
| `p2_isomerizacion_margules.py` | equilibrio líquido con γ (Margules) | `core.nonideal.margules_gamma` |
| `p3_gas_natural_hhv.py` | ΔH°c / HHV de mezclas | `core.energia.combustion_mezcla` |
| `p4_reformado_calor.py` | Q con 2 reacciones, salida conocida | `core.energia.balance_energia` |
| `p5_calentador_flujo_gas.py` | combustión incompleta, Q al agua, V̇ | `core.energia` + `h_liquido_sat` |

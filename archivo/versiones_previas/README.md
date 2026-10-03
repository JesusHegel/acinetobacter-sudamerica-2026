# Versiones previas

Scripts reemplazados durante el trabajo. Se conservan para que las decisiones
metodológicas queden documentadas.

| Script | Por qué se reemplazó |
|---|---|
| `curar_metadatos.py` | Sustituido por `pipeline/03_tipificacion/01_curar_metadatos.py` |
| `agrupamiento_andino.py` | Emplea enlace simple, descartado tras comprobar que produce agrupamientos cuyo diámetro interno supera el umbral declarado. El criterio adoptado es el enlace completo |
| `figura_arbol_st2.py` | Dibujaba el árbol como dendrograma. Sustituido por `pipeline/05_poblacion/03_red_expansion_minima_st2.py`, que lo representa como red de expansión mínima |

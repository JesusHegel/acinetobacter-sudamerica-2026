# Análisis descartado

Análisis completados y verificados que no forman parte del manuscrito. Se
conservan porque el trabajo se realizó y porque pueden sustentar una
publicación independiente.

## Partición de varianza

`C1_datos_modelo.py`, `C3_modelo_completo.R`, `C4_figura2.R`

Modelo multinivel que estima qué fracción de la variación observada entre
genomas corresponde al país de origen y cuál al proyecto que los depositó.

## Bootstrap del cociente de varianzas

`bootstrap_var.py`

Dos mil réplicas sobre el cociente de varianzas intranacional e
internacional. **El estimador resultó inestable y no se reportó.** El script
se conserva porque documenta esa decisión.

## Validación en *Klebsiella pneumoniae*

`klebsiella/`

Réplica del análisis de partición de varianza sobre 3447 genomas de otra
especie, para comprobar que el fenómeno no era exclusivo de
*A. baumannii*. Incluye su propio `pipeline.sh` y su documentación.

# Validación en *Klebsiella pneumoniae*

Aplicación del mismo procedimiento de estratificación por proyecto de origen
a un segundo organismo, como comprobación de que el fenómeno descrito para
*Acinetobacter baumannii* no es específico de esa especie.

## Conjunto

| | |
|---|---|
| Registros recuperados (NCBI Pathogen Detection, 9 países) | 4151 |
| Con ensamblado depositado | 3943 |
| Clínicos con ensamblado | 3639 |
| Confirmados como *K. pneumoniae* por ANI | 3447 |
| BioProjects con n ≥ 8 | 55 |
| Países con ≥ 2 proyectos | 7 |

De los 3639 analizados, 192 (5,3 %) resultaron pertenecer a otras especies
del complejo: *K. quasipneumoniae* (93), *K. aerogenes* (59) y
*K. variicola* (40).

## Variables

Las cuatro binarias equivalentes a las modeladas en *A. baumannii*:

| *A. baumannii* | *K. pneumoniae* |
|---|---|
| Sin carbapenemasa adquirida | Sin carbapenemasa adquirida |
| Complejo clonal 2 (ST2, ST2724) | Complejo clonal 258 (ST11, ST258, ST437, ST340) |
| bla<sub>OXA-23</sub> | bla<sub>KPC</sub> |
| bla<sub>OXA-72</sub> | bla<sub>NDM</sub> |

## Ficheros

| Fichero | Contenido |
|---|---|
| `pipeline.sh` | ANI (skani), MLST y AMRFinderPlus sobre los 3639 genomas |
| `modelo_kleb.R` | Modelo binomial multinivel en brms |
| `metadatos_kleb.tsv` | Metadatos descargados de NCBI Pathogen Detection |
| `ani_asignacion.tsv` | Especie asignada por ANI a cada genoma |
| `mlst_kleb.tsv` | Salida de mlst |
| `tabla_kleb.tsv` | Tabla curada: 3447 genomas con país, proyecto, ST y genes |
| `modelo_entrada_kleb.tsv` | 55 proyectos × 4 variables |
| `varianzas_kleb.tsv` | Coeficientes de partición de la varianza |

## Resultado

| Variable | VPC | P(VPC > 0,5) |
|---|---|---|
| Sin carbapenemasa | 0,879 | 0,984 |
| bla<sub>NDM</sub> | 0,790 | 0,949 |
| bla<sub>KPC</sub> | 0,695 | 0,876 |
| CC258 | 0,426 | 0,353 |

Los modelos convergieron sin transiciones divergentes (R̂ ≤ 1,003).


## Figura

`FigS1_klebsiella.pdf` — proporciones estimadas por proyecto para las cuatro
variables, agrupadas por país. Equivalente a la Figura 2 del manuscrito
principal. Generada con `figS_kleb.R`.

# Estratificación por proyecto de origen en genómica pública de *Acinetobacter baumannii*

Código y datos derivados del análisis de 900 genomas clínicos de *A. baumannii*
de nueve países sudamericanos depositados en NCBI Pathogen Detection.

## Contenido

| Ruta | Descripción |
|---|---|
| `EJECUTAR.md` | **Orden de ejecución, duración de cada paso y cifras de verificación** |
| `pipeline/` | Procedimiento completo, numerado por etapas |
| `archivo/` | Material que documenta el proceso pero no forma parte del análisis |
| `datos/` | Tablas de entrada: accesión→BioProject, longitudes de contig, metadatos |
| `resultados/` | Tablas de salida, matrices de distancia y figuras generadas |
| `entornos/` | Especificación conda de los entornos empleados |
| `figuras/` | Figuras del manuscrito en formato SVG y PDF |
| `metadatos/` | Tablas de metadatos curadas de NCBI Pathogen Detection |

### Etapas del `pipeline/`

| Carpeta | Qué hace |
|---|---|
| `01_datos/` | Descarga de los 900 genomas y de las referencias del complejo ACB |
| `02_ensamblado/` | Ensamblado de los registros disponibles solo como lecturas crudas |
| `03_tipificacion/` | Especie, MLST, cápsula, resistoma y replicones plasmídicos |
| `04_contexto/` | Contexto genético de las carbapenemasas: ISAba1, co-localización, sitios pdif |
| `05_poblacion/` | cgMLST, redes de expansión mínima y agrupamientos |
| `06_estadistica/` | Contrastes sobre eventos epidemiológicos independientes |
| `07_figuras/` | Figuras del manuscrito |
| `08_tablas/` | Tabla S1 y tablas suplementarias |

Cada script indica en su cabecera qué recibe, qué produce, qué entorno
necesita y cuánto tarda.

## Entornos

Tres entornos conda. Reconstruir con:

```bash
conda env create -f entornos/abaumannii.yml   # tipificación, anotación y figuras
conda env create -f entornos/ensamblaje.yml   # ensamblado y cgMLST
conda env create -f entornos/qc.yml           # control de calidad y ANI
```

| Herramienta | Versión | Entorno |
|---|---|---|
| Shovill | 1.4.2 | ensamblaje |
| SPAdes | 3.15.5 | ensamblaje |
| chewBBACA | 3.5.3 | ensamblaje |
| Prodigal | 2.6.3 | ensamblaje |
| CheckM2 | 1.1.0 | qc |
| skani | 0.3.2 | qc |
| BLAST+ | 2.17.0 | abaumannii |
| mlst | 2.35.0 | abaumannii |
| Kaptive | 3.2.2 | abaumannii |
| AMRFinderPlus | 4.2.7 (BD 2026-05-15.1) | abaumannii |
| SeqKit | 2.13.0 | abaumannii |
| ISEScan | 1.7.3 | abaumannii |
| GrapeTree | 2.2 | abaumannii |
| networkx | 3.6.1 | abaumannii |
| matplotlib | 3.11.1 | abaumannii |

## Datos de origen

Los genomas públicos se descargan de NCBI Pathogen Detection; las accesiones
están en `resultados/tabla_clinica_900.tsv`. Los 58 ensamblados generados aquí
proceden de lecturas depositadas por otros grupos en cinco BioProjects
(PRJEB107069, PRJNA1015678, PRJEB27899, PRJEB39593, PRJNA1322038) y se
depositaron enlazados a sus BioSamples originales; la correspondencia completa
está en `datos/metadatos_58_ensamblados.tsv`.

## Cita

[Referencia del artículo, pendiente]

## Licencia

Código bajo licencia MIT. Datos bajo CC BY 4.0.

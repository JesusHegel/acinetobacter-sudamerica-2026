# Estratificación por proyecto de origen en genómica pública de *Acinetobacter baumannii*

Código y datos derivados del análisis de 900 genomas clínicos de *A. baumannii*
de nueve países sudamericanos depositados en NCBI Pathogen Detection.

## Contenido

| Ruta | Descripción |
|---|---|
| `scripts/` | Scripts de análisis (Python, bash, awk) |
| `datos/` | Tablas de entrada: accesión→BioProject, longitudes de contig, metadatos de los 58 ensamblados |
| `resultados/` | Tablas de salida: tipificación de los 900 genomas, co-localización gen–replicón, matrices de distancia cgMLST |
| `entornos/` | Especificación conda de los entornos empleados |
| `figuras/` | Las siete figuras del manuscrito en formato SVG y PDF |
| `metadatos/` | Tablas de metadatos curadas de NCBI Pathogen Detection |
| `ensamblados_58.tar.gz` | Los 58 ensamblados generados en este trabajo (solo en el archivo de Zenodo) |

## Entornos

Tres entornos conda. Reconstruir con:

```bash
conda env create -f entornos/abaumannii.yml   # tipificación y anotación
conda env create -f entornos/ensamblaje.yml   # ensamblado y cgMLST
conda env create -f entornos/qc.yml           # control de calidad y ANI
```

| Herramienta | Versión | Entorno |
|---|---|---|
| Shovill | 1.4.2 | ensamblaje |
| SPAdes | 3.15.5 | ensamblaje |
| chewBBACA | 3.5.4 | ensamblaje |
| CheckM2 | 1.1.0 | qc |
| skani | 0.3.2 | qc |
| BLAST+ | 2.17.0 | abaumannii |
| mlst | 2.35.0 | abaumannii |
| Kaptive | 3.2.2 | abaumannii |
| AMRFinderPlus | 4.2.7 (BD 2026-05-15.1) | abaumannii |
| SeqKit | 2.13.0 | abaumannii |
| ISEScan | 1.7.3 | abaumannii |

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

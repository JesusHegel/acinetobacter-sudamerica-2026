# Estratificación por proyecto de origen en genómica pública de *Acinetobacter baumannii*

Código y datos derivados del análisis de 900 genomas clínicos de *A. baumannii*
de nueve países sudamericanos depositados en NCBI Pathogen Detection.

## Contenido

| Ruta | Descripción |
|---|---|
| `scripts/` | Scripts de análisis (Python, R, bash, awk) |
| `datos/` | Tablas de entrada: accesión→BioProject, longitudes de contig, metadatos de los 58 ensamblados |
| `resultados/` | Tablas de salida: tipificación de los 900 genomas, co-localización gen–replicón, contexto de ISAba1, matrices de distancia cgMLST |
| `resultados/grapetree/` | Red de expansión mínima de los 213 genomas ST2: perfiles, árbol en formato Newick y figuras coloreadas por país, año, BioProject y replicón portador |
| `resultados/grapetree_900/` | Red de expansión mínima de los 900 genomas sobre el esquema cgMLST de 2133 loci (umbral de presencia 0,95) |
| `resultados/plasmidos/` | Anotación de los elementos plasmídicos portadores de blaOXA-72 y del armazón chileno cerrado, con la figura comparativa |
| `resultados/figuras_rawgraphs/` | Gráficos de círculos anidados (país → BioProject → ST) generados con RawGraphs |
| `resultados/figuras_grapetree_web/` | Árboles exportados desde la interfaz de GrapeTree; las versiones reproducibles están en `resultados/grapetree/` |
| `entornos/` | Especificación conda de los entornos empleados |
| `figuras/` | Figuras del manuscrito en formato SVG y PDF |
| `metadatos/` | Tablas de metadatos curadas de NCBI Pathogen Detection |
| `klebsiella/` | Validación del análisis de partición de varianza en *K. pneumoniae* |
| `ensamblados_58.tar.gz` | Los 58 ensamblados generados en este trabajo (solo en el archivo de Zenodo) |

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

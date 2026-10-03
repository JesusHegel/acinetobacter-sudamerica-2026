# Resultados

Salidas del análisis descrito en `EJECUTAR.md`. Permiten comparar una
reproducción independiente sin necesidad de rehacer los pasos previos.

## Qué hay

| Ruta | Contenido |
|---|---|
| `CIFRAS_MANUSCRITO.txt` | Valores reportados en el texto, con su procedencia |
| `TablaS1_900_genomas.tsv` | Los 900 genomas con metadatos y tipificación completa |
| `tipificacion/` | ANI, especie asignada, CheckM2, Kaptive (KL y OCL) y AMRFinderPlus |
| `cgmlst_900_eval/` | Perfiles alélicos de los 2390 loci y matriz al umbral de 0,95 |
| `cgmlst_st2_eval/`, `cgmlst_st79_eval/` | Matrices de distancia por linaje |
| `colocalizacion.tsv`, `apt_hits_v2.tsv` | Replicones detectados y co-localización con los genes de resistencia |
| `isaba1_*.tsv` | Contexto de ISAba1 por familia de carbapenemasa |
| `plasmidos/` | Anotación de los elementos portadores y del armazón chileno |
| `grapetree/`, `grapetree_900/` | Árboles en formato Newick, metadatos y figuras |
| `figuras_rawgraphs/`, `figuras_grapetree_web/` | Figuras exportadas desde herramientas interactivas |

## Qué no hay

Se omiten las carpetas de trabajo intermedias, que pesan en conjunto más de
2 GB y se regeneran al ejecutar el procedimiento:

- Ficheros intermedios de CheckM2; solo se incluye su `quality_report`
- Coordenadas de CDS y datos por contig de chewBBACA; solo se incluyen los perfiles
- Registro completo del ensamblado con SPAdes

## Cómo comparar

Las cifras de verificación están en `EJECUTAR.md`. Para una comparación
fichero a fichero, los dos que concentran el resultado son
`TablaS1_900_genomas.tsv` y `cgmlst_900_eval/cgMLST95.tsv`.

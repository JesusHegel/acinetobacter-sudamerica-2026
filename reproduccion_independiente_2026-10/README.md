# Reproducción independiente — octubre de 2026

Dos reejecuciones completas del pipeline (pasos 01 a 08) sobre una máquina
limpia (Windows 11 + WSL2 + Debian 13), realizadas por una persona que no
participó en el análisis original, partiendo únicamente de este repositorio
y de `EJECUTAR.md`.

| Pasada | Fecha | Sobre | Resultado |
|---|---|---|---|
| Primera | 11 de octubre | el repositorio tal como estaba | reprodujo las cifras y destapó los defectos de empaquetado |
| Segunda | 27 al 29 de octubre | commit `a3d4f8e`, ya corregido | certificación: todo corre sin intervención |

La primera pasada está en la raíz de esta carpeta. La segunda, en
`certificacion_a3d4f8e/`.

## Qué encontró la primera

Defectos que solo aparecen en una máquina donde las herramientas no conviven
en un mismo entorno y donde no hay ficheros sueltos de ejecuciones anteriores.
Entre ellos, dos ficheros intermedios del paso 05 que ningún script producía
aunque uno anunciara haberlos escrito, y varias rutas que apuntaban fuera del
repositorio. Todos están corregidos.

## Qué certifica la segunda

Pipeline completo desde un clon nuevo, en un directorio base distinto
(`ABAUMANNII_BASE`), con los 842 genomas públicos y el esquema cgMLST
descargados otra vez. Nada reutilizado de la primera.

Siete tablas comparadas fila a fila con las publicadas, **las siete
idénticas**: TablaS1 (900 filas), TablaS2 (22), TablaS3 (33), TablaS4 (7),
TablaS5 (58), `tabla_clinica_900` (900) y el anexo de ISAba1 (44). La matriz
de distancias cgMLST de los ST2 coincide en sus **22 791 pares**. El árbol de
los 213 ST2 comparte **las 71 particiones**. El de los 900 comparte 303 de
307, y ninguna de las ocho discrepantes contiene un genoma de Perú o Ecuador.

Cifras del hallazgo reproducidas: 73 portadores de blaOXA-72, 63
co-localizados; CC2 andino de 24 genomas, los 24 con r3-T18; Brasil 35,
ninguno con r3-T18; nueve eventos andinos; hueco KL2/KL9 de 21 a 30;
Fisher p = 3 × 10⁻⁷ sobre eventos.

## Entorno

Esquema cgMLST de 2390 loci de Chewie-NS
(`http://chewbbaca.online/NS/api/species/2/schemas/1`, versión
`2025-07-31T14:24:53`). AMRFinderPlus 4.2.7 con base 2026-08-07.1 —tres meses
posterior a la del análisis original—, chewBBACA 3.5.4, mlst 2.35.0, Kaptive
3.2.2, BLAST 2.17.0+, SeqKit 2.13.0, skani 0.3.2. Detalle en los dos
`entorno.txt`.

El cambio de catálogo de AMRFinderPlus no movió ninguna cifra: las 68 líneas
adicionales que aporta no son carbapenemasas. Y los nueve falsos negativos de
blaOXA-72 siguen sin detectarse con él, lo que confirma que su causa es la
linealización del elemento circular y no un catálogo desactualizado.

## Verificación de integridad

Cada pasada trae un `MANIFEST.md5` con la huella de sus ficheros, calculada en
la máquina donde se produjeron:

    md5sum -c MANIFEST.md5

El `.gitattributes` de esta carpeta declara `* -text` para que git no convierta
los finales de línea, aquí y en las subcarpetas. Sin esa línea, una copia
clonada en Windows recibiría los `.tsv` y los `.txt` alterados y los
manifiestos dejarían de cuadrar. No conviene retirarlo.

## Commits posteriores a la certificación

La certificación cubre `a3d4f8e`. Los commits posteriores no alteran el
análisis:

| Commit | Qué cambia | Efecto sobre una reproducción |
|---|---|---|
| `a801916` | Etiquetas impresas del bloque D del `04/10` y entrada del `03/02` | Las etiquetas son texto en pantalla |
| `a593635` | Filtro por los 900 en el cruce del `03/02`; dos ficheros de `resultados/` regenerados | Ver abajo |
| `5cb9b04` | Duraciones declaradas en esta guía | Ninguno |

Los dos ficheros regenerados (`apt_hits_v2.tsv` y `colocalizacion.tsv`) viven
en `repo/resultados/`, que ningún paso del pipeline lee: los scripts leen de
`$B/resultados/`, el directorio de trabajo.

El filtro añadido al `03/02` descarta los genomas que no están entre los 900.
En una ejecución limpia `amrfinder_todos.tsv` contiene solo esos 900, de modo
que no hay nada que descartar: la reejecución independiente obtuvo 1354 y 878
filas con el script anterior al filtro, que son las cifras correctas. Los 87 y
31 sobrantes de los ficheros publicados procedían del directorio de trabajo de
la máquina de desarrollo, que arrastraba salidas del conjunto anterior de 944
genomas, no del script.

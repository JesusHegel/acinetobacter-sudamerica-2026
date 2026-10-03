# Reproducción del análisis

Este documento indica el orden de ejecución, lo que produce cada paso y
cuánto tarda. Los números de los ficheros marcan el orden dentro de cada
carpeta.

## Requisitos previos

**Espacio en disco:** 50 GB libres recomendados.

**Ubicación del proyecto.** Los scripts emplean rutas absolutas bajo
`~/abaumannii/`. El repositorio debe clonarse en `~/abaumannii/repo` y los
datos se escribirán en `~/abaumannii/datos` y `~/abaumannii/resultados`:

```bash
mkdir -p ~/abaumannii && cd ~/abaumannii
git clone https://github.com/JesusHegel/acinetobacter-sudamerica-2026.git repo
```

**La herramienta `datasets` del NCBI** vive en el entorno `abaumannii`. Los
pasos que la emplean lo indican en su cabecera.

**Entornos conda.** Tres entornos, reconstruibles desde `entornos/`:

```bash
conda env create -f entornos/abaumannii.yml   # tipificación, anotación, figuras
conda env create -f entornos/ensamblaje.yml   # ensamblado y cgMLST
conda env create -f entornos/qc.yml           # control de calidad y ANI
```

Cada paso indica en su cabecera qué entorno necesita.

---

## 01 · Datos de partida

| Paso | Qué hace | Entorno | Duración |
|---|---|---|---|
| `01_descargar_genomas.sh` | Descarga los 842 ensamblados públicos a partir de la lista fija de accesiones | abaumannii | 30-60 min |
| `02_descargar_referencias_acb.sh` | Las cinco referencias del complejo ACB para la confirmación de especie | abaumannii | 1 min |
| `03_descargar_esquemas.sh` | Catálogo de AMRFinderPlus, esquema APT y esquema cgMLST | abaumannii + ensamblaje | 20-40 min |

La lista `accesiones_900.txt` contiene los 900 identificadores del conjunto:
842 ensamblados públicos (GCA) y 58 registros del SRA que se ensamblan en
el paso 02.

**No se repite la consulta original a NCBI Pathogen Detection.** Las bases
públicas crecen con el tiempo, de modo que una consulta nueva recuperaría
un conjunto distinto. La lista fija garantiza el mismo punto de partida.

---

## 02 · Ensamblado

| Paso | Qué hace | Entorno | Duración |
|---|---|---|---|
| `00_descargar_lecturas.sh` | Descarga del SRA las lecturas crudas de los 58 registros | ensamblaje | 2-4 h |
| `01_ensamblar_lecturas.sh` | Ensambla esas lecturas con SPAdes | ensamblaje | 6-10 h |

Se procesaron 63 registros, de los cuales 58 superaron el control de calidad
del paso 03 e integran el conjunto final. La lista de accesiones está en
`accesiones_sra_58.txt`.

### Este paso es opcional

El ensamblado de los 63 registros añade entre 6 y 10 horas de cómputo y
verifica el funcionamiento de SPAdes, no las conclusiones del trabajo. Hay
dos formas de proceder:

**Opción A, recomendada: omitir el ensamblado.** Los 58 ensamblados que
resultaron de este paso se han remitido a GenBank bajo el BioProject
PRJNA1505778, **cuya liberación queda condicionada a la publicación del
trabajo**. Hasta entonces no son descargables desde el NCBI, y deben
solicitarse a los autores. Se entregan como un único archivo comprimido que
se extrae en `datos/ensamblados_63/`:

```bash
mkdir -p ~/abaumannii/datos/ensamblados_63
tar -xzf ensamblados_58.tar.gz -C ~/abaumannii/datos/ensamblados_63 --strip-components=1
ls ~/abaumannii/datos/ensamblados_63/*.fna | wc -l   # debe dar 58
```

**Opción B: reensamblar.** Solo si el objetivo es comprobar también el
procedimiento de ensamblado. En ese caso, los ensamblados obtenidos pueden
diferir mínimamente de los originales en el número de contigs, sin que ello
afecte a la tipificación.

**En ambos casos los 900 genomas deben estar presentes** antes de continuar
con el paso 03. Si se omite el ensamblado sin descargar los 58, el conjunto
quedará incompleto y ninguna de las cifras de verificación coincidirá.

---

## 03 · Tipificación

| Paso | Qué hace | Entorno | Duración |
|---|---|---|---|
| `00_tipificar.sh` | Especie (skani), MLST, cápsula (Kaptive) y resistoma (AMRFinderPlus) | qc + abaumannii | 2-4 h |
| `01_curar_metadatos.py` | Normaliza país, año y fuente de aislamiento | base | min |
| `02_tipificar_plasmidos.sh` | Replicones del esquema APT y co-localización con los genes de resistencia | abaumannii | 1-2 h |

---

## 04 · Contexto genético

| Paso | Qué hace | Entorno | Duración |
|---|---|---|---|
| `01_coordenadas_carbapenemasas.py` | Extrae las coordenadas de cada gen desde la salida de AMRFinderPlus | base | min |
| `02_cruce_isaba1.awk` | Cruza esas coordenadas con las posiciones de ISAba1 | base | min |
| `03_anexo_isaba1_oxa23.py` | Anexo de los genomas con ISAba1 río arriba de blaOXA-23 | base | min |
| `04_rastreo_dirigido.sh` | Busca blaOXA-72 sin umbral de cobertura, para recuperar los partidos por el punto de linealización | abaumannii | 20 min |
| `05_informe_rastreo.py` | Compara lo recuperado frente a lo notificado | base | min |
| `06_recalcular_colocalizacion.py` | Recalcula el contexto exigiendo replicón en el mismo contig | base | min |
| `07_sitios_pdif.py` | Localiza los sitios pdif en los elementos | base | min |
| `08_comparar_plasmidos_cerrados.py` | Compara el elemento andino con plásmidos circulares cerrados | abaumannii | min |
| `09_armazon_con_y_sin_gen.py` | Contrasta el armazón r3-T18 vacío frente al portador | base | min |
| `10_r3t18_fuera_del_cc2.py` | Distribución de r3-T18 en los linajes no andinos | base | min |

---

## 05 · Estructura poblacional

| Paso | Qué hace | Entorno | Duración |
|---|---|---|---|
| `00_cgmlst.sh` | AlleleCall y ExtractCgMLST sobre los 900 genomas | ensamblaje | 20 min |
| `01_arbol_grapetree.sh` | Red de expansión mínima, algoritmo MSTreeV2 | abaumannii | 2 min |
| `02_agrupamiento_enlace_completo.py` | Agrupamiento jerárquico y barrido de umbrales | abaumannii | min |
| `03_red_expansion_minima_st2.py` | Red de resolución fina sobre ST2 | abaumannii | min |

---

## 06 · Estadística

| Paso | Qué hace | Entorno | Duración |
|---|---|---|---|
| `01_fisher_exclusividad.py` | Prueba exacta de Fisher sobre eventos independientes | base | seg |

---

## 07 · Figuras y 08 · Tablas

Ambos requieren el entorno `abaumannii`. Se ejecutan en cualquier orden, una
vez completados los pasos anteriores.
Las tablas se generan en el orden numerado: la S1 primero y después sus
correcciones.

---

## Cifras de verificación

Si la reproducción es correcta, estos valores deben coincidir:

| Resultado | Valor |
|---|---|
| Genomas clínicos tras los filtros | 900 |
| Eventos epidemiológicos independientes | 360 |
| Genomas ST2 | 213 |
| Loci cgMLST retenidos al umbral 0,95 | 2133 de 2390 |
| Portadores de blaOXA-72 | 73 |
| Portadores con replicón tipificado | 63 (86,3 %) |
| blaOXA-23 co-localizado con replicón | 5 de 600 (0,83 %) |
| Genomas con ISAba1 río arriba de blaOXA-23 | 44 |
| Replicón r3-T18 en el conjunto | 46 genomas, 6 países |
| Contraste principal (Fisher) | p = 3 × 10⁻⁷ |

Una discrepancia en cualquiera de estos valores indica un problema que
conviene localizar antes de continuar.

---

## Qué hay en `archivo/`

Material que documenta el proceso pero no forma parte del análisis
publicado. Se conserva por trazabilidad.

- `auditoria/` — verificaciones de consistencia realizadas durante el trabajo
- `analisis_descartado/` — partición de varianza y validación en *K. pneumoniae*, no incluidas en el manuscrito
- `versiones_previas/` — scripts reemplazados por versiones posteriores
- `deposito/` — preparación del envío de los ensamblados a GenBank

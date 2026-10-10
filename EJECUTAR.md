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
| `02_unificar_conjunto.sh` | **Obligatorio en ambas opciones.** Reúne los 842 públicos y los 58 propios en un único directorio | cualquiera | segundos |

El paso `02_unificar_conjunto.sh` debe ejecutarse siempre, tanto si se
reensambla como si se extrae el archivo comprimido. Los pasos 03 y 05 operan
sobre `datos/genomas_900/`, que es el directorio que crea. Si se omite, el
paso 03 falla en su primera orden.

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
| `01_coordenadas_carbapenemasas.py` | Extrae las coordenadas de cada gen desde las salidas por genoma de AMRFinderPlus | base | min |
| `02_buscar_isaba1.sh` | Localiza ISAba1 en los 900 genomas por BLASTN, y calcula las longitudes de contig | abaumannii | 15-30 min |
| `03_cruce_isaba1.sh` | Cruza las posiciones de ISAba1 con las coordenadas de cada gen y clasifica cada copia | base | seg |
| `04_anexo_isaba1_oxa23.py` | Anexo de los genomas con ISAba1 rio arriba de blaOXA-23 | base | min |
| `05_rastreo_dirigido.sh` | Busca blaOXA-72 sin umbral de cobertura | abaumannii | 20 min |
| `06_informe_rastreo.py` | Compara lo recuperado frente a lo notificado. Debe indicar 9 portadores nuevos | base | min |
| `06b_propagar_rastreo.py` | Escribe los 9 recuperados en la tabla clinica y en la Tabla S1. El `06` solo informa; sin este paso ambas quedan en 64 portadores | base | seg |
| `07_recalcular_colocalizacion.py` | Recalcula el contexto exigiendo replicon en el mismo contig | base | min |
| `08_sitios_pdif.py` | Localiza los sitios pdif en los elementos | base | min |
| `09_comparar_plasmidos_cerrados.py` | Compara el elemento andino con plasmidos circulares cerrados | abaumannii | min |
| `10_armazon_con_y_sin_gen.py` | Contrasta el armazon r3-T18 vacio frente al portador | base | min |
| `11_r3t18_fuera_del_cc2.py` | Distribucion de r3-T18 en los linajes no andinos | base | min |

El paso `03_cruce_isaba1.sh` invoca el fichero `03_cruce_isaba1.awk`, que no se ejecuta por si solo.

**Orden de ejecución.** Los pasos `07` y `11` leen la Tabla S1, que produce `08_tablas/01_tabla_S1.py`. El orden ejecutable es:

```
03  ->  08/01  ->  04  ->  08/02  ->  05/00b  ->  05/00c  ->  08/04  ->  06  ->  07 y 08/03
```

**`05/00b` y `05/00c` van despues de la Tabla S1.** El `00b` selecciona los 213
ST2 leyendo la columna `st_pasteur` de la Tabla S1, de modo que no puede
correrse antes de `08/01`. El `00c` necesita la matriz que produce el `00b`, y
`08/04` necesita esa misma matriz para el agrupamiento del CC2 andino.

**`08/04` va antes que `06`.** `04_recuento_final.py` escribe el numero de eventos
andinos en `resultados/verif/eventos_andinos.txt`, que lee el Fisher. Si falta, el
Fisher avisa y usa el valor publicado (9).

**Dentro del paso 04 el orden numerado no es el de ejecucion.** El `01` lee `datos/contig_len_900.tsv`, que produce el `02`, de modo que el orden real es `02 → 01 → 03 → 04 → …`. Los pasos `01`, `04`, `07` y `11` leen ademas la Tabla S1, que produce `08_tablas/01_tabla_S1.py`.

## 05 · Estructura poblacional

| Paso | Qué hace | Entorno | Duración |
|---|---|---|---|
| `00_cgmlst.sh` | AlleleCall y ExtractCgMLST sobre los 900 genomas | ensamblaje | 20 min |
| `00b_cgmlst_st2.sh` | AlleleCall y ExtractCgMLST restringidos a los 213 ST2. Produce la matriz de distancias del linaje | ensamblaje | 5 min |
| `00c_arbol_y_metadatos_st2.sh` | Red de expansión mínima de los ST2 y sus metadatos de anotación | abaumannii | 1 min |
| `01_arbol_grapetree.sh` | Red de expansión mínima, algoritmo MSTreeV2 | abaumannii | 2 min |
| `02_agrupamiento_enlace_completo.py` | Agrupamiento jerárquico y barrido de umbrales | abaumannii | min |
| `03_red_expansion_minima_st2.py` | Red de resolución fina sobre ST2 | abaumannii | min |

Los pasos `00b` y `00c` no existían hasta el 10 de octubre de 2026: la matriz
de ST2, su árbol y sus metadatos se habían generado a mano y no eran
reproducibles. Los consumen `02`, `03`, `07_figuras/fig5_dendrograma_andino.py`
y `08_tablas/04_recuento_final.py`. El `00c` reproduce el árbol y los metadatos
publicados byte a byte.

---

## 06 · Estadística

| Paso | Qué hace | Entorno | Duración |
|---|---|---|---|
| `01_fisher_exclusividad.py` | Prueba exacta de Fisher sobre eventos independientes | base | seg |

El numero de eventos andinos ya no va fijo en el codigo: lo toma de
`resultados/verif/eventos_andinos.txt`, que produce `08_tablas/04_recuento_final.py`.

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
| Genomas sin carbapenemasa adquirida | 133 |
| Genomas con carbapenemasa adquirida | 767 |
| Copias de blaOXA-23 evaluadas para ISAba1 | 696, evaluabilidad 7,3 % |
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

---

## Nota sobre una correccion del analisis

Las cifras de esta guia corresponden al analisis corregido el 4 de octubre
de 2026. Una version anterior reportaba **176 genomas sin carbapenemasa
adquirida**; la cifra correcta es **139**.

**Que ocurrio.** El procedimiento que construia la tabla de tipificacion
buscaba las carbapenemasas por nombre literal (`blaOXA-23`, `blaOXA-72`,
`blaOXA-58`, `blaOXA-143`), de modo que las variantes de esas mismas
familias registradas con otro nombre no se contabilizaban. El apartado 2.4
de metodos declara que el recuento se hace por familia y no por variante,
pero ese criterio solo se habia aplicado a la familia intrinseca.

**Alcance.** Cuarenta y cuatro genomas pasan de figurar sin carbapenemasa a
portador, todos con coincidencia exacta y 100 % de identidad y cobertura:
38 con `blaOXA-253`, 4 con `blaOXA-366` y 2 con `blaOXA-657`. La familia
OXA-143 pasa de 1 a 41 genomas. Se corrigen ademas tres asignaciones
erroneas: `blaOXA-407` y `blaOXA-241` pertenecen a la familia intrinseca y
`blaGES-11` es una cefalosporinasa, no una carbapenemasa.

**Lo que no cambia.** El hallazgo central no se altera. blaOXA-72 sigue
presente en 73 genomas, 63 con replicon co-localizado, y la separacion entre
el replicon andino y el brasileno mantiene su significacion (p = 3,2e-07).

---

## Segunda correccion: propagacion del rastreo dirigido (10 de octubre de 2026)

Las cifras de genomas con y sin carbapenemasa adquirida vuelven a moverse:
**139 -> 133** sin carbapenemasa y **761 -> 767** con carbapenemasa.

**Que ocurrio.** El rastreo dirigido de blaOXA-72 (`04_contexto/05` y `06`)
recupera 9 portadores que AMRFinderPlus etiqueta como `blaOXA` generico por
tratarse de alineamientos partidos entre fragmentos de un contig
(`PARTIAL_CONTIG_ENDX`). El paso `06` solo imprimia el resultado: nadie lo
escribia. La columna `oxa72_contexto` de la Tabla S1 si recogia los 9, pero
`genes_carbapenemasa` y `carbapenemasa_adquirida` se quedaban en el recuento
de AMRFinderPlus, de modo que la propia tabla se contradecia (73 frente a 64).

**Alcance.** Nueve genomas suman `blaOXA-72` en `genes_carbapenemasa`
(4 brasilenos y 5 peruanos); de ellos, seis pasan ademas de `no` a `si` en
`carbapenemasa_adquirida`. Todos con identidad BLASTN de 99,5 a 100 % frente
a la referencia, por encima del umbral de 95 % del rastreo. Lo escribe el
paso nuevo `04_contexto/06b_propagar_rastreo.py`, que es idempotente.

**Lo que no cambia.** Los 73 portadores de blaOXA-72, los 63 con replicon
co-localizado (86,3 %) y el contraste principal (p = 3,2e-07) se mantienen.
Las dos tablas son ahora coherentes entre si.

**Tercer arreglo de la misma tanda.** El numero de eventos andinos ya no va
fijo en el codigo del Fisher: lo escribe `08_tablas/04_recuento_final.py` en
`resultados/verif/eventos_andinos.txt`. Y los seis escritores CSV del
pipeline llevan ya `lineterminator="\n"`, de modo que las tablas no salen
con fin de linea CRLF (afectaba a la Tabla S1 y a la tabla clinica, y hacia
que la ultima columna arrastrase un retorno de carro al leerla con awk).

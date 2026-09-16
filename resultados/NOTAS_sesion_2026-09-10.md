# Resultados nuevos — sesión del 10 de septiembre de 2026

Trabajo derivado de la segunda reunión con la Dra. Pons.

---

## 1. Concentración del depósito: PRJNA613847

**Un solo BioProject brasileño aporta 89 de los 213 genomas ST2 del
continente: el 42 %.**

Al colorear el árbol por año de colección aparece un pico marcado en 2021,
con 95 de los 213 genomas ST2. **De esos 95, setenta y siete proceden de
PRJNA613847.** Los 18 restantes se reparten entre once proyectos distintos,
la mayoría con un solo genoma.

El pico no refleja una expansión del linaje: refleja un depósito.

Es la demostración concreta, sobre datos propios, de por qué el año de
colección se emplea como variable descriptiva y no como serie temporal.

Figuras: `grapetree/red_st2_anio.png` y `grapetree/red_st2_bioproject.png`.

---

## 2. ISAba1 en las carbapenemasas adquiridas

Mismo procedimiento aplicado a blaOXA-51-like (ventana de 500 pb río arriba,
BLASTN de la referencia de 1180 pb con identidad >= 95 % y alineamiento >= 300 pb).

| Determinante            | Copias | Presente | Ausente | No eval. | Evaluab. |
|-------------------------|-------:|---------:|--------:|---------:|---------:|
| blaOXA-51-like (intr.)  |   1002 |        7 |     834 |      161 |   83,9 % |
| blaOXA-23               |    696 |       51 |       0 |      645 |    7,3 % |
| blaOXA-24/40            |     89 |        0 |      69 |       20 |   77,5 % |
| blaOXA-58               |     49 |        0 |      13 |       36 |   26,5 % |
| blaOXA-143              |     47 |        0 |      23 |       24 |   48,9 % |
| blaNDM                  |     44 |        0 |      30 |       14 |   68,2 % |

**ISAba1 se detecta río arriba únicamente de blaOXA-23 y del gen intrínseco.
Cero casos en las otras cuatro familias, sobre 135 copias evaluables.**

Los 51 hits de blaOXA-23 corresponden a 44 genomas de 8 países, 10 linajes y
22 BioProjects. Ningún genoma peruano: los cinco peruanos con blaOXA-23
quedaron como no evaluables.

La evaluabilidad del 7,3 % en blaOXA-23 no es un fallo del análisis. Refleja
que el gen se sitúa habitualmente a menos de 500 pb del borde de su contig,
consecuencia de viajar en transposones flanqueados por ISAba1 que el
ensamblador de lecturas cortas no resuelve.

Anexo por genoma: `anexo_isaba1_oxa23.tsv`.

---

## 3. Anotación del elemento andino r3-T18

Prodigal en modo `meta` (Prokka falla con secuencias < 20 000 pb).
Doce CDS sobre 8111 pb.

- **rep r3-T18** (596-1546, 317 aa) — 100 % de identidad con la referencia
  `r3-T18_CP042561.1_pE47_005_c53` del catálogo APT
- **mobQ** (2225-3880, 552 aa) — 92,6 % con `mobQ-T36`
- **blaOXA-72** (7150-7977, 276 aa) — 100 % de identidad y cobertura contra
  el catálogo de AMRFinderPlus. **Único CDS del elemento que corresponde a un
  gen de resistencia.**
- Seis proteínas hipotéticas
- Una proteína partida por el punto de linealización (CDS 1 y 12)
- **Sin parA**: el elemento carece de sistema de partición

---

## 4. Comparación con el armazón chileno

p3UC20804 (CP076810.1), 8229 pb, circular cerrado, mismo tipo r3-T18, sin el gen.
Rotado 595 posiciones para situar su rep en la misma fase que el andino.

**Nueve de los doce CDS andinos tienen equivalente en el chileno. Seis son
idénticos al 100 %. El único CDS ausente en el chileno es blaOXA-72.**

| Andino    | Chileno | Identidad |
|-----------|---------|----------:|
| andino_2  | _2      |    97,9 % |
| andino_3  | _1      |    99,4 % |
| andino_4  | _11     |    88,2 % |
| andino_5  | _10     |     100 % |
| andino_6  | _9      |     100 % |
| andino_7  | _8      |     100 % |
| andino_8  | _7      |     100 % |
| andino_9  | _6      |     100 % |
| andino_10 | _5      |     100 % |

Esta comparación gen por gen sustituye con ventaja al argumento previo basado
en porcentajes de cobertura recíproca (83-85 %).

Figura: `plasmidos/figura_plasmidos.png`.

---

## Pendiente

- Figura del panorama de los 900 genomas por MLST (equivalente a la Figura 6
  de Hummel et al., Microorganisms 2024).
- Decidir con la Dra. Pons el alcance del primer artículo.

---

## 5. Panorama de linajes (Figura A)

goeBURST sobre los 900 genomas clínicos, esquema Pasteur de 7 loci.
Equivalente a la Figura 6 de Hummel et al. (Microorganisms 2024), que emplea
PHYLOViZ Online sobre 452 genomas europeos.

- 73 tipos de secuencia distintos; 72 con perfil alélico completo
- 27 conexiones de variante de locus único
- **ST2724 aparece unido a ST2 por una arista**, lo que confirma por
  perfil alélico su pertenencia al complejo clonal 2

Figura: `grapetree/goeburst_900.png`.
Script: `scripts/goeburst_900.py`.

---

## 5. Figuras de panorama

### Círculos anidados (RawGraphs)

Equivalente a la Figura 7 de Hummel et al., pero con **tres niveles en lugar de
dos**: país → BioProject → tipo de secuencia, con el área proporcional al número
de genomas.

La versión coloreada por BioProject muestra que dentro de Brasil un único
proyecto concentra una fracción desproporcionada del total, lo que hace visible
el problema de agregación que el trabajo cuantifica.

Datos: `circulos_pais_bioproject_st.csv` (263 filas, 900 genomas).
Figuras: `figuras_rawgraphs/circulos_por_pais.svg` y
`figuras_rawgraphs/circulos_por_bioproject.svg`.

### Árboles de ST2 (GrapeTree)

Exportados desde la interfaz de GrapeTree sobre el mismo árbol MSTreeV2,
coloreados por país, año y plásmido portador. Las versiones reproducibles,
generadas por script, están en `grapetree/red_st2_*.png`.

Figuras: `figuras_grapetree_web/MSTree_ST2_*.svg`.

### Descartado

Se probó un goeBURST de 7 loci sobre los 900 genomas (nodo = ST, arista =
variante de locus único). Se descartó: solo 27 de los 73 STs tienen algún
vecino a un locus de distancia, de modo que la figura resulta poco informativa
y no corresponde a ninguna de las figuras del trabajo de referencia.

La Figura 6 de Hummel et al. se basa en perfiles cgMLST de todos los aislados,
no en MLST de 7 loci. Replicarla exigiría ejecutar chewBBACA sobre los 900
genomas. Queda pendiente de decidir con la Dra. Pons si resulta necesaria.

---

## 6. cgMLST sobre los 900 genomas (Figura A definitiva)

Se ejecutó chewBBACA AlleleCall sobre los 900 genomas clínicos con el esquema
completo de 2390 loci (BSR 0,6, `--no-inferred`). Duración: 13 min 12 s.
Se clasificaron 2 163 586 CDS, de los cuales 2 070 439 fueron coincidencias
exactas.

ExtractCgMLST sobre la matriz resultante:

| Umbral de presencia | Loci retenidos |
|---------------------|---------------:|
| 0,95                |   **2133**/2390 |
| 0,99                |       1588/2390 |
| 1,00                |        312/2390 |

**Se adopta el umbral de 0,95, con 2133 loci.** Coincide con el tamaño de
esquema empleado por Hummel et al. (2024) en su Figura 6.

Sobre esa matriz se construyó una red de expansión mínima con GrapeTree
(MSTreeV2): 900 hojas, 1208 nodos totales. Es el equivalente sudamericano de
la Figura 6 de Hummel et al., que reúne 452 aislados europeos.

Nota: el conjunto de ST2 analizado previamente retiene 1719 loci sobre 213
genomas. Que los 900 retengan más (2133) no es contradictorio: el umbral de
0,95 tolera ausencia en un 5 % de los genomas por locus, mientras que la
matriz de ST2 se calculó sobre loci presentes en la totalidad del subconjunto.

Figura: `grapetree_900/red_900_por_st.png`.
Perfiles: `cgmlst_900_eval/cgMLST95.tsv` · esquema: `cgMLSTschema95.txt`.

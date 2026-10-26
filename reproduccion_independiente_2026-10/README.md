# Reproducción independiente — octubre de 2026

Reejecución completa del pipeline (pasos 01 a 08) sobre una máquina limpia
(Windows 11 + WSL2 + Debian 13), realizada por una persona que no participó en
el análisis original, partiendo únicamente de este repositorio y de
`EJECUTAR.md`.

## Contenido

- `registros/` — salida en pantalla de siete pasos, tal como se produjo: sitios
  pdif (04/08), armazón con y sin gen (04/10), r3-T18 fuera del CC2 (04/11),
  agrupamiento por enlace completo (05/02), prueba de Fisher (06/01), tablas
  suplementarias (08/03) y recuento final (08/04).
- `scripts/` — los siete scripts de verificación escritos durante la
  reejecución, más dos reescrituras de trabajo (paso 03 reanudable y rastreo de
  ISAba1 sobre disco en lugar de tmpfs).
- `salidas/` — ficheros regenerados: tablas S1 a S5, tabla clínica, matriz y
  perfiles cgMLST de los 213 ST2, árboles de los 213 y de los 900, cgMLST95 de
  los 900 y figuras.
- `entorno.txt` — sistema, entornos conda con versiones, base de datos de
  AMRFinderPlus y esquema cgMLST empleados.

## Entorno

Esquema cgMLST de 2390 loci obtenido de Chewie-NS
(`http://chewbbaca.online/NS/api/species/2/schemas/1`, versión
`2025-07-31T14:24:53`), descargado el 4 de octubre de 2026. chewBBACA 3.5.4,
AMRFinderPlus 4.2.7 con base 2026-08-07.1, mlst 2.35.0, Kaptive 3.2.2, BLAST
2.17.0+, SeqKit 2.13.0, skani 0.3.2. Detalle completo en `entorno.txt`.

## Verificación de integridad

`MANIFEST.md5` contiene la huella de los 49 ficheros, calculada en la máquina
donde se produjeron. Para comprobar que no han cambiado:

    md5sum -c MANIFEST.md5

El `.gitattributes` de esta carpeta declara `* -text` para que git no convierta
los finales de línea. Sin esa línea, una copia clonada en Windows recibiría los
`.tsv` y los `.txt` alterados y el manifiesto dejaría de cuadrar. No conviene
retirarlo.

## Resultado

Las cifras que sostienen el trabajo se reprodujeron. Comparaciones directas con
los ficheros publicados: la matriz de distancias de los 213 ST2 coincide en sus
22 791 pares, el árbol de los ST2 comparte las 71 particiones, y las tablas S4 y
S5 son idénticas.

Diferencias encontradas, documentadas en `EJECUTAR.md`: los perfiles cgMLST de
los 900 genomas difieren en 128 celdas de 1 919 700 (0,007 %), concentradas en
7 loci, lo que desplaza 15 847 de 404 550 pares con un máximo de 3 loci sobre
2133; y el árbol de expansión mínima de los 900 comparte 302 de sus 307
particiones, ninguna de las discrepantes con genomas de Perú o Ecuador.

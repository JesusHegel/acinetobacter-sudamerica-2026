#!/usr/bin/env bash
# =====================================================================
# Tipificacion de genoma central sobre los 900 genomas.
#
# Entrada:  datos/genomas_900/
#           datos/cgmlst_abaumannii/   (esquema de 2390 loci, Chewie-NS)
# Produce:  resultados/cgmlst_900/      perfiles alelicos
#           resultados/cgmlst_900_eval/ matrices por umbral de presencia
# Entorno:  ensamblaje
#
# Este script NO activa el entorno por si mismo: un script de bash no puede
# cambiar el entorno conda del proceso que lo invoca. Ejecutelo asi:
#
#   conda activate ensamblaje && bash 00_cgmlst.sh
# Duracion: 13 min (AlleleCall) + 5 min (ExtractCgMLST)
#
# El esquema se descarga de Chewie-NS:
#   chewBBACA.py DownloadSchema -sp 2 -sc 1 -o datos/cgmlst_abaumannii
# =====================================================================
set -euo pipefail
BASE=~/abaumannii

echo "=== 1/2  AlleleCall ==="
chewBBACA.py AlleleCall \
  -i "$BASE/datos/genomas_900" \
  -g "$BASE/datos/cgmlst_abaumannii/Acinetobacter_baumannii_cgMLST" \
  -o "$BASE/resultados/cgmlst_900" \
  --cpu 4 --bsr 0.6 --no-inferred

echo "=== 2/2  ExtractCgMLST ==="
# Retiene los loci presentes en >= 95 % de los genomas (2133 de 2390).
# Produce tambien las matrices a 0,99 y 1,00 para comprobar el efecto del umbral.
chewBBACA.py ExtractCgMLST \
  -i "$BASE/resultados/cgmlst_900/results_alleles.tsv" \
  -o "$BASE/resultados/cgmlst_900_eval"

echo "listo"

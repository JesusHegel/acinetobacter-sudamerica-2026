#!/usr/bin/env bash
# =====================================================================
# Tipificacion de los 900 genomas: especie, MLST, capsula y resistoma.
#
# Entrada:  datos/genomas_900/        (pasos 01 y 02)
# Produce:  resultados/tipificacion/
# Duracion: 2 a 4 h, dominada por AMRFinderPlus
#
# IMPORTANTE: este paso emplea VARIOS entornos conda y un unico script de
# bash no puede cambiar entre ellos. Ejecutelo en tandas:
#
#   conda activate qc          && bash 00_tipificar.sh especie
#   conda activate abaumannii  && bash 00_tipificar.sh resto
#
# Sin argumento corre los cuatro bloques, que es lo correcto cuando un
# mismo entorno tiene todas las herramientas.
#
# Si reconstruyo los entornos desde entornos/minimos/, las dos llamadas a
# mlst del bloque 2 necesitan ademas el entorno mlst, porque mlst 2.35.0
# arrastra blast 2.16.0 y el resto del analisis usa blast 2.17.0:
#
#   conda activate abaumannii  && bash 00_tipificar.sh sinmlst
#   conda activate mlst        && bash 00_tipificar.sh mlst
#
# En la maquina de desarrollo ambas versiones conviven en abaumannii, de
# modo que alli el bloque 2 corre con el resto.
#
# Argumentos admitidos:
#   todo      los cuatro bloques (es lo que hace sin argumento)
#   especie   1. confirmacion de especie con skani
#   mlst      2. tipificacion de secuencia multilocus
#   resto     2, 3 y 4
#   sinmlst   3 y 4
#
# El bloque de AMRFinderPlus omite los genomas ya procesados, de modo que
# una interrupcion no obliga a repetir el trabajo hecho.
# =====================================================================
set -euo pipefail
BASE=~/abaumannii
IN=$BASE/datos/genomas_900
OUT=$BASE/resultados/tipificacion
REF=$BASE/datos/referencias_acb
mkdir -p "$OUT"

# Seleccion de bloques. Sin argumento corre los cuatro, que es lo correcto
# cuando un unico entorno tiene todas las herramientas.
QUE="${1:-todo}"
case "$QUE" in
  todo)    HAZ_ESPECIE=1; HAZ_MLST=1; HAZ_KAPTIVE=1; HAZ_AMR=1 ;;
  especie) HAZ_ESPECIE=1; HAZ_MLST=0; HAZ_KAPTIVE=0; HAZ_AMR=0 ;;
  mlst)    HAZ_ESPECIE=0; HAZ_MLST=1; HAZ_KAPTIVE=0; HAZ_AMR=0 ;;
  resto)   HAZ_ESPECIE=0; HAZ_MLST=1; HAZ_KAPTIVE=1; HAZ_AMR=1 ;;
  sinmlst) HAZ_ESPECIE=0; HAZ_MLST=0; HAZ_KAPTIVE=1; HAZ_AMR=1 ;;
  *) echo "uso: bash 00_tipificar.sh [todo|especie|mlst|resto|sinmlst]"; exit 1 ;;
esac
echo "bloques solicitados: $QUE"

# ---------------------------------------------------------------------
# 0. Bases de datos            (solo la primera vez)
# ---------------------------------------------------------------------
# AMRFinderPlus requiere descargar su catalogo antes del primer uso.
# Este trabajo empleo la version 2026-05-15.1; una version posterior
# puede notificar variantes adicionales.
#   amrfinder -u
#
# El esquema cgMLST del paso 05 se descarga de Chewie-NS:
#   chewBBACA.py DownloadSchema -sp 2 -sc 1 -o ~/abaumannii/datos/cgmlst_abaumannii

# ---------------------------------------------------------------------
# 1. Confirmacion de especie            entorno: qc
# ---------------------------------------------------------------------
# ANI frente a las cinco referencias del complejo A. calcoaceticus-baumannii.
# Umbral de 95 % para delimitacion de especie; se asigna la referencia con
# mayor identidad. Las referencias estan en 01_datos/referencias_acb/
if [ "$HAZ_ESPECIE" -eq 1 ]; then
echo "=== 1/4  Confirmacion de especie (skani) ==="
ls "$IN"/*.fna > "$OUT/lista_genomas.txt"
ls "$REF"/*.fna > "$OUT/lista_referencias.txt"
skani dist --ql "$OUT/lista_genomas.txt" --rl "$OUT/lista_referencias.txt" \
  -t 4 -o "$OUT/ani_skani.tsv"

fi

# ---------------------------------------------------------------------
# 2. Tipificacion de secuencia multilocus          entorno: abaumannii
#                                           (o mlst, ver la cabecera)
# ---------------------------------------------------------------------
# abaumannii_2 = esquema Pasteur (adoptado para el analisis)
# abaumannii   = esquema Oxford  (solo para comparar cobertura)
if [ "$HAZ_MLST" -eq 1 ]; then
echo "=== 2/4  MLST ==="
mlst --scheme abaumannii_2 "$IN"/*.fna > "$OUT/mlst_pasteur.tsv"
mlst --scheme abaumannii   "$IN"/*.fna > "$OUT/mlst_oxford.tsv"

fi

# ---------------------------------------------------------------------
# 3. Tipificacion capsular                         entorno: abaumannii
# ---------------------------------------------------------------------
# Las bases vienen incluidas en el propio paquete de Kaptive.
KDB=$(python3 -c "import kaptive,os;print(os.path.join(os.path.dirname(kaptive.__file__),'data'))")
if [ "$HAZ_KAPTIVE" -eq 1 ]; then
echo "=== 3/4  Kaptive  (bases en $KDB) ==="
kaptive assembly "$KDB/Acinetobacter_baumannii_k_locus_primary_reference.gbk" \
  "$IN"/*.fna -o "$OUT/kaptive_kl.tsv"
kaptive assembly "$KDB/Acinetobacter_baumannii_OC_locus_primary_reference.gbk" \
  "$IN"/*.fna -o "$OUT/kaptive_ocl.tsv"

fi

# ---------------------------------------------------------------------
# 4. Genes de resistencia                          entorno: abaumannii
# ---------------------------------------------------------------------
if [ "$HAZ_AMR" -eq 1 ]; then
echo "=== 4/4  AMRFinderPlus ==="
# Se produce un fichero por genoma en analisis_amr4/, ademas del fichero
# concatenado. El paso 04 localiza las coordenadas de cada gen por genoma,
# de modo que necesita las salidas individuales.
mkdir -p "$BASE/analisis_amr4"
for f in "$IN"/*.fna; do
  n=$(basename "$f" .fna)
  n=$(echo "$n" | sed -E 's/^(GC[AF]_[0-9]+\.[0-9]+).*/\1/')
  # Se omite lo ya calculado, de modo que una interrupcion no obligue a
  # repetir el trabajo hecho. El fichero por genoma se escribe SIN --name,
  # para que conserve el formato de columnas que espera el paso 04.
  if [ ! -s "$BASE/analisis_amr4/$n.tsv" ]; then
    amrfinder -n "$f" --organism Acinetobacter_baumannii --plus \
      > "$BASE/analisis_amr4/$n.tsv"
  fi
done

# El fichero concatenado se reconstruye a partir de las salidas por genoma,
# anteponiendo la accesion a cada linea. Asi queda completo aunque la
# ejecucion se haya reanudado.
: > "$OUT/amrfinder_todos.tsv"
for g in "$BASE"/analisis_amr4/*.tsv; do
  a=$(basename "$g" .tsv)
  tail -n +2 "$g" | awk -v A="$a" -F'\t' 'BEGIN{OFS="\t"}{print A, $0}' \
    >> "$OUT/amrfinder_todos.tsv"
done
echo "genomas con salida de AMRFinderPlus: $(ls "$BASE"/analisis_amr4/*.tsv | wc -l)"

fi

echo
echo "Salidas en $OUT"
ls -la "$OUT"

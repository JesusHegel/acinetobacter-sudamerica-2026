#!/usr/bin/env bash
# =====================================================================
# Tipificacion de los 900 genomas: especie, MLST, capsula y resistoma.
#
# Entrada:  datos/genomas_900/        (pasos 01 y 02)
# Produce:  resultados/tipificacion/
# Duracion: 2 a 4 h, dominada por AMRFinderPlus
#
# IMPORTANTE: este paso emplea DOS entornos conda y un unico script de bash
# no puede cambiar entre ellos. Ejecutelo en dos tandas:
#
#   conda activate qc          && bash 00_tipificar.sh especie
#   conda activate abaumannii  && bash 00_tipificar.sh resto
#
# Sin argumento ejecuta solo la parte que corresponde al entorno activo.
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

# ---------------------------------------------------------------------
# 0. Bases de datos            (solo la primera vez)
# ---------------------------------------------------------------------
# AMRFinderPlus requiere descargar su catalogo antes del primer uso.
# Este trabajo empleo la version 2026-05-15.1; una version posterior
# puede notificar variantes adicionales.
#   amrfinder -u
#
# El esquema cgMLST del paso 05 se descarga de Chewie-NS:
#   chewBBACA.py DownloadSchema -sp 5 -sc 1 -o ~/abaumannii/datos/cgmlst_abaumannii

# ---------------------------------------------------------------------
# 1. Confirmacion de especie            entorno: qc
# ---------------------------------------------------------------------
# ANI frente a las cinco referencias del complejo A. calcoaceticus-baumannii.
# Umbral de 95 % para delimitacion de especie; se asigna la referencia con
# mayor identidad. Las referencias estan en 01_datos/referencias_acb/
echo "=== 1/4  Confirmacion de especie (skani) ==="
ls "$IN"/*.fna > "$OUT/lista_genomas.txt"
ls "$REF"/*.fna > "$OUT/lista_referencias.txt"
skani dist --ql "$OUT/lista_genomas.txt" --rl "$OUT/lista_referencias.txt" \
  -t 4 -o "$OUT/ani_skani.tsv"

# ---------------------------------------------------------------------
# 2. Tipificacion de secuencia multilocus          entorno: abaumannii
# ---------------------------------------------------------------------
# abaumannii_2 = esquema Pasteur (adoptado para el analisis)
# abaumannii   = esquema Oxford  (solo para comparar cobertura)
echo "=== 2/4  MLST ==="
mlst --scheme abaumannii_2 "$IN"/*.fna > "$OUT/mlst_pasteur.tsv"
mlst --scheme abaumannii   "$IN"/*.fna > "$OUT/mlst_oxford.tsv"

# ---------------------------------------------------------------------
# 3. Tipificacion capsular                         entorno: abaumannii
# ---------------------------------------------------------------------
# Las bases vienen incluidas en el propio paquete de Kaptive.
KDB=$(python3 -c "import kaptive,os;print(os.path.join(os.path.dirname(kaptive.__file__),'data'))")
echo "=== 3/4  Kaptive  (bases en $KDB) ==="
kaptive assembly "$KDB/Acinetobacter_baumannii_k_locus_primary_reference.gbk" \
  "$IN"/*.fna -o "$OUT/kaptive_kl.tsv"
kaptive assembly "$KDB/Acinetobacter_baumannii_OC_locus_primary_reference.gbk" \
  "$IN"/*.fna -o "$OUT/kaptive_ocl.tsv"

# ---------------------------------------------------------------------
# 4. Genes de resistencia                          entorno: abaumannii
# ---------------------------------------------------------------------
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

echo
echo "Salidas en $OUT"
ls -la "$OUT"

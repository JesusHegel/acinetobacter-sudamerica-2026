#!/usr/bin/env bash
# =====================================================================
# Tipificacion de genoma central restringida a los 213 genomas ST2.
#
# El analisis de estructura fina del CC2 andino (seccion 3.9) necesita una
# matriz de distancias calculada solo sobre el linaje: con los 900 genomas
# el umbral de presencia del 95 % retiene loci distintos y las distancias
# intra-linaje no son comparables.
#
# Entrada:  datos/genomas_900/
#           datos/cgmlst_abaumannii/   (esquema de 2390 loci, Chewie-NS)
#           resultados/TablaS1_900_genomas.tsv   (para seleccionar los ST2)
# Produce:  resultados/cgmlst_st2/       perfiles alelicos
#           resultados/cgmlst_st2_eval/  distance_matrix_symmetric.tsv
# Entorno:  ensamblaje
#
#   conda activate ensamblaje && bash 00b_cgmlst_st2.sh
# Duracion: 4 min (AlleleCall) + 1 min (ExtractCgMLST) + segundos (matriz)
#
# El paso 4 lo hace 00b2_matriz_distancias_st2.py, en este mismo
# directorio: ExtractCgMLST no calcula distancias y los dos ficheros que
# consumen el 00c, el 05/02, fig5_dendrograma_andino.py y el 08/04 no los
# producia ningun paso.
# =====================================================================
set -euo pipefail
BASE=${ABAUMANNII_BASE:-$HOME/abaumannii}
LISTA="$BASE/resultados/genomas_st2.txt"
DIR="$BASE/datos/genomas_st2"

echo "=== 1/4  seleccionando los ST2 de la Tabla S1 ==="
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="st_pasteur") c=i; next}
            $c=="2"{print $1}' "$BASE/resultados/TablaS1_900_genomas.tsv" > "$LISTA"
echo "    genomas ST2: $(wc -l < "$LISTA")"

# Los ficheros de datos/genomas_900 llevan la accesion como nombre, sin
# sufijo de ensamblado; se enlazan en vez de copiarse.
rm -rf "$DIR"; mkdir -p "$DIR"
n=0; faltan=0
while read -r acc; do
  f="$BASE/datos/genomas_900/${acc}.fna"
  if [ -f "$f" ]; then ln -s "$f" "$DIR/${acc}.fna"; n=$((n+1))
  else echo "    FALTA: $acc"; faltan=$((faltan+1)); fi
done < "$LISTA"
echo "    enlazados: $n   faltantes: $faltan"
[ "$faltan" -eq 0 ] || { echo "ERROR: faltan genomas en datos/genomas_900"; exit 1; }

echo "=== 2/4  AlleleCall ==="
chewBBACA.py AlleleCall \
  -i "$DIR" \
  -g "$BASE/datos/cgmlst_abaumannii/Acinetobacter_baumannii_cgMLST" \
  -o "$BASE/resultados/cgmlst_st2" \
  --cpu 4 --bsr 0.6 --no-inferred

echo "=== 3/4  ExtractCgMLST ==="
chewBBACA.py ExtractCgMLST \
  -i "$BASE/resultados/cgmlst_st2/results_alleles.tsv" \
  -o "$BASE/resultados/cgmlst_st2_eval"

echo "=== 4/4  perfiles y matriz de distancias ==="
# ExtractCgMLST no calcula distancias: hay que derivarlas del cgMLST100.tsv.
python3 "$(dirname "$0")/00b2_matriz_distancias_st2.py" \
  "$BASE/resultados/cgmlst_st2_eval"

for f in cgMLST_profiles.tsv distance_matrix_symmetric.tsv; do
  r="$BASE/resultados/cgmlst_st2_eval/$f"
  [ -s "$r" ] || { echo "ERROR: no se produjo $f"; exit 1; }
done
echo "listo"

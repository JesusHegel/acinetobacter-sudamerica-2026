#!/usr/bin/env bash
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate abaumannii
set -u

BASE=~/abaumannii
IN=$BASE/datos/genomas_900
OUT=$BASE/resultados/tipificacion
REF=$BASE/datos/referencias_acb
AMR=$BASE/analisis_amr4
mkdir -p "$OUT" "$AMR"

if [ -s "$OUT/ani_skani.tsv" ]; then
  echo "=== 1/4  skani: ya estaba hecho, se omite ==="
else
  echo "=== 1/4  Confirmacion de especie, skani  ($(date +%H:%M)) ==="
  ls "$IN"/*.fna > "$OUT/lista_genomas.txt"
  ls "$REF"/*.fna > "$OUT/lista_referencias.txt"
  conda run --no-capture-output -n qc skani dist \
    --ql "$OUT/lista_genomas.txt" --rl "$OUT/lista_referencias.txt" \
    -t 4 -o "$OUT/ani_skani.tsv"
fi

if [ -s "$OUT/mlst_pasteur.tsv" ]; then
  echo "=== 2/4  MLST: ya estaba hecho, se omite ==="
else
  echo "=== 2/4  MLST  ($(date +%H:%M)) ==="
  mlst --scheme abaumannii_2 "$IN"/*.fna > "$OUT/mlst_pasteur.tsv"
  mlst --scheme abaumannii   "$IN"/*.fna > "$OUT/mlst_oxford.tsv"
fi

if [ -s "$OUT/kaptive_kl.tsv" ]; then
  echo "=== 3/4  Kaptive: ya estaba hecho, se omite ==="
else
  echo "=== 3/4  Kaptive  ($(date +%H:%M)) ==="
  KDB=$(python3 -c "import kaptive,os;print(os.path.join(os.path.dirname(kaptive.__file__),'data'))")
  kaptive assembly "$KDB/Acinetobacter_baumannii_k_locus_primary_reference.gbk" \
    "$IN"/*.fna -o "$OUT/kaptive_kl.tsv"
  kaptive assembly "$KDB/Acinetobacter_baumannii_OC_locus_primary_reference.gbk" \
    "$IN"/*.fna -o "$OUT/kaptive_ocl.tsv"
fi

echo "=== 4/4  AMRFinderPlus  ($(date +%H:%M)) ==="
total=$(ls "$IN"/*.fna | wc -l); i=0; nuevos=0
for f in "$IN"/*.fna; do
  i=$((i+1))
  n=$(basename "$f" .fna)
  n=$(echo "$n" | sed -E 's/^(GC[AF]_[0-9]+\.[0-9]+).*/\1/')
  if [ -s "$AMR/$n.tsv" ]; then continue; fi
  amrfinder -n "$f" --organism Acinetobacter_baumannii --plus --name "$n" > "$AMR/$n.tsv"
  nuevos=$((nuevos+1))
  if [ $((nuevos % 25)) -eq 0 ]; then echo "   $i de $total  ($(date +%H:%M))"; fi
done
echo "   procesados en esta tanda: $nuevos"

echo "=== reconstruyendo el fichero concatenado ==="
: > "$OUT/amrfinder_todos.tsv"
for t in "$AMR"/*.tsv; do tail -n +2 "$t" >> "$OUT/amrfinder_todos.tsv"; done

echo
echo "ficheros individuales de AMRFinder: $(ls "$AMR"/*.tsv | wc -l)  (deben ser $total)"
echo "lineas en el concatenado:           $(wc -l < "$OUT/amrfinder_todos.tsv")"
echo "TERMINADO  ($(date +%H:%M))"

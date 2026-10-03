#!/bin/bash
READS=~/abaumannii/datos/reads
OUT=~/abaumannii/datos/ensamblados_63
TMP=~/abaumannii/datos/tmp_shovill
mkdir -p "$OUT" "$TMP"

for r1 in "$READS"/*_1.fastq.gz; do
  run=$(basename "$r1" _1.fastq.gz)
  r2="$READS/${run}_2.fastq.gz"
  [ -f "$OUT/${run}.fna" ] && { echo "ya: $run"; continue; }
  [ -f "$r2" ] || { echo "sin par: $run"; continue; }

  echo "=== $run  $(date +%H:%M) ==="
  shovill --R1 "$r1" --R2 "$r2" \
          --outdir "$OUT/tmp_${run}" \
          --tmpdir "$TMP" \
          --cpus 6 --ram 12 --minlen 200 --force

  if [ -s "$OUT/tmp_${run}/contigs.fa" ]; then
    mv "$OUT/tmp_${run}/contigs.fa" "$OUT/${run}.fna"
  else
    echo "$run" >> "$OUT/fallos_ensamblaje.txt"
    echo ">>> FALLO: $run"
  fi
  rm -rf "$OUT/tmp_${run}" "$TMP"/*
done

echo "--- ensamblados: $(ls "$OUT"/*.fna 2>/dev/null | wc -l) ---"

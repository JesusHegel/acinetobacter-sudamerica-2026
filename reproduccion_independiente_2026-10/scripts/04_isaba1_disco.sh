#!/usr/bin/env bash
# Version de 02_buscar_isaba1.sh que no agota la memoria:
#   - el temporal va en disco, no en /tmp (que es tmpfs)
#   - base BLAST indexada en vez de -subject
#   - no silencia los errores
set -eo pipefail
B=~/abaumannii
REF=$B/repo/datos/isaba1_ref.fa
OUT=$B/datos/isaba1_hits900.tsv
TMP=$B/tmp_isaba1
mkdir -p "$TMP"

if [ "$(wc -l < "$B/datos/contig_len_900.tsv" 2>/dev/null || echo 0)" -ge 100000 ]; then
  echo "longitudes de contig: ya estaban ($(wc -l < "$B/datos/contig_len_900.tsv") contigs)"
else
  echo "calculando longitudes de contig..."
  : > "$B/datos/contig_len_900.tsv"
  for f in "$B"/datos/genomas_900/*.fna; do
    acc=$(basename "$f" .fna)
    seqkit fx2tab -nl "$f" | awk -v A="$acc" -F'\t' '{print A"\t"$1"\t"$NF}' \
      >> "$B/datos/contig_len_900.tsv"
  done
  echo "contigs: $(wc -l < "$B/datos/contig_len_900.tsv")"
fi

if [ ! -s "$TMP/all.fna" ]; then
  echo "concatenando genomas en disco..."
  : > "$TMP/all.fna"
  for f in "$B"/datos/genomas_900/*.fna; do
    acc=$(basename "$f" .fna)
    awk -v A="$acc" '/^>/{sub(/^>/,">"A"|");print;next}{print}' "$f" >> "$TMP/all.fna"
  done
fi
echo "secuencias: $(grep -c '^>' "$TMP/all.fna")"

if [ ! -s "$TMP/db.nin" ] && [ ! -s "$TMP/db.nal" ]; then
  echo "construyendo base BLAST indexada..."
  makeblastdb -in "$TMP/all.fna" -dbtype nucl -out "$TMP/db" > /dev/null
fi

echo "BLAST de ISAba1 ($(date +%H:%M))..."
blastn -query "$REF" -db "$TMP/db" \
       -outfmt "6 sseqid sstart send pident length" \
       -evalue 1e-5 -max_target_seqs 100000 -num_threads 4 \
  | awk -F'\t' '$4>=95 && $5>=300' \
  | awk -F'\t' '{split($1,a,"|"); print a[1]"\t"a[2]"\t"$2"\t"$3"\t"$4"\t"$5}' > "$OUT"

echo "hits de ISAba1: $(wc -l < "$OUT")"
echo "genomas con al menos un hit: $(cut -f1 "$OUT" | sort -u | wc -l)"
echo "TERMINADO ($(date +%H:%M))"

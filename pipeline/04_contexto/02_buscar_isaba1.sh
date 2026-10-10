#!/usr/bin/env bash
# =====================================================================
# Localiza la secuencia de insercion ISAba1 en los 900 genomas.
#
# BLASTN de la referencia de ISAba1 (1180 pb) contra el conjunto completo,
# conservando los alineamientos con identidad >= 95 % y longitud >= 300 pb.
#
# Entrada:  datos/genomas_900/        (paso 02)
#           datos/isaba1_ref.fa       referencia de ISAba1
# Produce:  datos/isaba1_hits900.tsv
#           datos/contig_len_900.tsv
# Entorno:  conda activate abaumannii
# Duracion: 15 a 30 min
# =====================================================================
set -euo pipefail
B=~/abaumannii
REF=$B/datos/isaba1_ref.fa
OUT=$B/datos/isaba1_hits900.tsv
# El temporal va bajo el directorio del proyecto y no en /tmp, que en WSL
# y en otras instalaciones es un sistema de ficheros en memoria: aqui se
# concatenan 3,5 GB de secuencia.
TMP=$(mktemp -d "$B/tmp_isaba1.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# longitudes de contig, que el cruce necesita para decidir evaluabilidad
: > "$B/datos/contig_len_900.tsv"
for f in "$B"/datos/genomas_900/*.fna; do
  acc=$(basename "$f" .fna)
  seqkit fx2tab -nil "$f" | awk -v A="$acc" -F'\t' '{print A"\t"$1"\t"$NF}' \
    >> "$B/datos/contig_len_900.tsv"
done
echo "contigs: $(wc -l < "$B/datos/contig_len_900.tsv")"

# una sola base con todos los contigs, prefijados por accesion
: > "$TMP/all.fna"
for f in "$B"/datos/genomas_900/*.fna; do
  acc=$(basename "$f" .fna)
  awk -v A="$acc" '/^>/{sub(/^>/,">"A"|");print;next}{print}' "$f" >> "$TMP/all.fna"
done

# Se indexa el sujeto en lugar de pasarlo con -subject: ese modo carga la
# secuencia entera en memoria y agota la RAM en maquinas de 8 GB.
makeblastdb -in "$TMP/all.fna" -dbtype nucl -out "$TMP/db" >/dev/null

blastn -query "$REF" -db "$TMP/db" \
       -outfmt "6 sseqid sstart send pident length" \
       -evalue 1e-5 -max_target_seqs 100000 \
  | awk -F'\t' '$4>=95 && $5>=300' \
  | awk -F'\t' '{split($1,a,"|"); print a[1]"\t"a[2]"\t"$2"\t"$3"\t"$4"\t"$5}' > "$OUT"

echo "hits de ISAba1: $(wc -l < "$OUT")"
echo "genomas con al menos un hit: $(cut -f1 "$OUT" | sort -u | wc -l)"

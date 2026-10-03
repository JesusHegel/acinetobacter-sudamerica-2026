#!/usr/bin/env bash
# =====================================================================
# Descarga las lecturas crudas de los 58 registros que el conjunto
# incorpora como ensamblados propios.
#
# Requiere: sra-tools (prefetch, fasterq-dump)
# Entorno:  conda activate ensamblaje
# Produce:  datos/reads/<accesion>_1.fastq.gz y _2.fastq.gz
# Duracion: 2 a 4 h segun conexion
# Espacio:  unos 25 GB
#
# Este paso solo es necesario si se opta por reensamblar (Opcion B del
# documento EJECUTAR.md).
# =====================================================================
set -euo pipefail
BASE=~/abaumannii
LISTA=$BASE/repo/pipeline/02_ensamblado/accesiones_sra_58.txt
OUT=$BASE/datos/reads
mkdir -p "$OUT" "$BASE/tmp_sra"
cd "$BASE/tmp_sra"

n=0
total=$(wc -l < "$LISTA")
while read -r acc; do
  [ -z "$acc" ] && continue
  n=$((n+1))
  if [ -f "$OUT/${acc}_1.fastq.gz" ]; then
    echo "[$n/$total] $acc ya descargado"
    continue
  fi
  echo "[$n/$total] $acc"
  prefetch "$acc" --max-size 50G
  fasterq-dump "$acc" --split-files --threads 4 --outdir "$OUT"
  gzip -f "$OUT/${acc}_1.fastq" "$OUT/${acc}_2.fastq"
  rm -rf "$acc"
done < "$LISTA"

echo
echo "pares descargados: $(ls "$OUT"/*_1.fastq.gz 2>/dev/null | wc -l)"

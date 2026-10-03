#!/usr/bin/env bash
# =====================================================================
# Descarga las cinco referencias del complejo A. calcoaceticus-baumannii
# empleadas en la confirmacion de especie por ANI (paso 03).
#
# Accesiones verificadas contra el NCBI el 3 de octubre de 2026.
#
# Requiere: datasets (NCBI), unzip, seqkit
# Entorno:  conda activate abaumannii
# Produce:  datos/referencias_acb/   (~29 MB)
# Duracion: 1 min
# =====================================================================
set -euo pipefail
DEST=~/abaumannii/datos/referencias_acb
TMP=/tmp/ref_acb
mkdir -p "$DEST" "$TMP"

declare -A REF=(
  [baumannii]=GCF_009035845.1       # ATCC 19606,   NZ_CP045110.1
  [nosocomialis]=GCF_041021905.1    # XH1679,       NZ_CP157432.1
  [pittii]=GCF_000191145.1          # PHEA-2,       NC_016603.1
  [seifertii]=GCF_016064815.1       # S21,          NZ_CP065820.1
  [calcoaceticus]=GCF_900444805.1   # NCTC 12983,   NZ_UFSJ01000002.1
)

for sp in "${!REF[@]}"; do
  echo "descargando $sp (${REF[$sp]})"
  datasets download genome accession "${REF[$sp]}" --include genome \
    --filename "$TMP/$sp.zip" >/dev/null
  unzip -q -o "$TMP/$sp.zip" -d "$TMP/$sp"
  # Se descartan los plasmidos de la cepa, que no formaban parte de las
  # referencias empleadas. El filtro actua sobre la descripcion de cada
  # secuencia, de modo que respeta los genomas en borrador cuyos contigs
  # no estan anotados como plasmido.
  find "$TMP/$sp" -name "*.fna" -exec cat {} \; \
    | seqkit grep -nrvip "plasmid" > "$DEST/$sp.fna"
done

rm -rf "$TMP"
echo
echo "referencias en $DEST"
for f in "$DEST"/*.fna; do
  echo "  $(basename "$f")  $(head -1 "$f" | cut -c2-60)"
done

# Nota: el fichero resultante puede diferir del original en el ancho de
# linea del FASTA. Las secuencias son identicas y el calculo de ANI no
# se ve afectado.

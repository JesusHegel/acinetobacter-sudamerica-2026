#!/usr/bin/env bash
# =====================================================================
# Descarga los 842 ensamblados publicos del conjunto, a partir de la
# lista fija de accesiones.
#
# No se repite la consulta original a NCBI Pathogen Detection: las bases
# publicas crecen con el tiempo y una consulta nueva recuperaria un
# conjunto distinto. La lista fija garantiza que se obtengan exactamente
# los mismos genomas analizados.
#
# Los 58 restantes del conjunto son ensamblados propios y se obtienen en
# el paso 02.
#
# Requiere: datasets (NCBI), unzip
# Entorno:  conda activate abaumannii
# Produce:  datos/genomas_945/ncbi_dataset/data/<accesion>/*.fna
# Duracion: 30 a 60 min segun conexion
#
# Nota sobre el nombre del directorio: se conserva "genomas_945" porque es
# el que emplean los pasos posteriores. Corresponde al conjunto descargado
# antes de los filtros de calidad, de los que resultaron 842 publicos.
# =====================================================================
set -euo pipefail
BASE=~/abaumannii
LISTA=$BASE/repo/pipeline/01_datos/accesiones_900.txt
DEST=$BASE/datos/genomas_945
TMP=$BASE/tmp_descarga
mkdir -p "$DEST" "$TMP"

grep '^GCA' "$LISTA" > "$TMP/gca.txt"
echo "ensamblados publicos a descargar: $(wc -l < "$TMP/gca.txt")"

datasets download genome accession --inputfile "$TMP/gca.txt" \
  --include genome --filename "$TMP/genomas.zip"

# Se conserva la jerarquia ncbi_dataset/data/<accesion>/ que produce
# datasets, porque los pasos 03, 04, 07 y 08 localizan cada genoma por
# su accesion dentro de esa estructura.
unzip -q -o "$TMP/genomas.zip" -d "$DEST"
rm -rf "$TMP"

n=$(find "$DEST/ncbi_dataset/data" -name "*.fna" | wc -l)
echo "descargados: $n"
[ "$n" -eq 842 ] || echo "AVISO: se esperaban 842 ficheros"

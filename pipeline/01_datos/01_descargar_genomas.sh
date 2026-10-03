#!/usr/bin/env bash
# =====================================================================
# Descarga los 900 genomas del conjunto a partir de su lista de accesiones.
#
# No se repite la consulta original a NCBI Pathogen Detection: las bases
# publicas crecen con el tiempo y una consulta nueva recuperaria un
# conjunto distinto. La lista fija garantiza que se obtengan exactamente
# los mismos 900 genomas analizados.
#
# Requiere: datasets (NCBI), unzip
# Entorno:  conda activate abaumannii
# Produce:  datos/genomas_900/   (~3.5 GB)
# Duracion: 30 a 60 min segun conexion
# =====================================================================
set -euo pipefail
BASE=~/abaumannii
LISTA=$BASE/repo/pipeline/01_datos/accesiones_900.txt
DEST=$BASE/datos/genomas_900
mkdir -p "$DEST" "$BASE/tmp_descarga"

echo "Genomas a descargar: $(grep -c GCA "$LISTA" || true)"
echo "Nota: los registros que no empiezan por GCA son ensamblados propios,"
echo "      generados en el paso 02 a partir de lecturas del SRA."

grep '^GCA' "$LISTA" > "$BASE/tmp_descarga/gca.txt"
datasets download genome accession --inputfile "$BASE/tmp_descarga/gca.txt" \
  --include genome --filename "$BASE/tmp_descarga/genomas.zip"
unzip -q -o "$BASE/tmp_descarga/genomas.zip" -d "$BASE/tmp_descarga/x"
find "$BASE/tmp_descarga/x" -name "*.fna" -exec cp {} "$DEST/" \;

echo "descargados: $(ls "$DEST" | wc -l)"

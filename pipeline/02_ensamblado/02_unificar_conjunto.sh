#!/usr/bin/env bash
# =====================================================================
# Reune en un solo directorio los 900 genomas del conjunto.
#
# Los pasos anteriores dejan los genomas en dos ubicaciones distintas:
#   - los 842 publicos en datos/genomas_945/ncbi_dataset/data/<accesion>/
#   - los 58 propios en datos/ensamblados_63/
#
# Los pasos 03 y 05 operan sobre el conjunto completo, de modo que se
# crean enlaces a ambos bloques en un unico directorio. Se emplean
# enlaces simbolicos y no copias para no duplicar 3,5 GB en disco.
#
# Entrada:  datos/genomas_945/ y datos/ensamblados_63/
# Produce:  datos/genomas_900/<accesion>.fna
# Entorno:  cualquiera
# Duracion: segundos
# =====================================================================
set -euo pipefail
BASE=~/abaumannii
DEST=$BASE/datos/genomas_900
mkdir -p "$DEST"

# Se enlaza unicamente lo que figura en la lista de accesiones, de modo que
# el conjunto queda definido por esa lista y no por lo que haya en disco.
LISTA=$BASE/repo/pipeline/01_datos/accesiones_900.txt
n=0; m=0; falta=0
while read -r acc; do
  [ -z "$acc" ] && continue
  if [[ "$acc" == GCA_* ]]; then
    f=$(ls "$BASE"/datos/genomas_945/ncbi_dataset/data/"$acc"/*.fna 2>/dev/null | head -1 || true)
    if [ -n "$f" ]; then ln -sf "$f" "$DEST/$acc.fna"; n=$((n+1))
    else echo "  ausente (publico): $acc"; falta=$((falta+1)); fi
  else
    f=$BASE/datos/ensamblados_63/$acc.fna
    if [ -f "$f" ]; then ln -sf "$f" "$DEST/$acc.fna"; m=$((m+1))
    else echo "  ausente (propio): $acc"; falta=$((falta+1)); fi
  fi
done < "$LISTA"

echo "publicos enlazados:  $n"
echo "propios enlazados:   $m"
echo "ausentes:            $falta"

t=$(ls "$DEST"/*.fna 2>/dev/null | wc -l)
echo "total en el conjunto: $t"
if [ "$t" -ne 900 ]; then
  echo
  echo "AVISO: se esperaban 900 genomas. Revise los pasos 01 y 02 antes de continuar."
  exit 1
fi

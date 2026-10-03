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

# publicos: se normaliza el nombre a <accesion>.fna
n=0
for f in "$BASE"/datos/genomas_945/ncbi_dataset/data/*/*.fna; do
  [ -e "$f" ] || continue
  acc=$(basename "$(dirname "$f")")
  ln -sf "$f" "$DEST/$acc.fna"
  n=$((n+1))
done
echo "publicos enlazados:  $n"

# propios
m=0
for f in "$BASE"/datos/ensamblados_63/*.fna; do
  [ -e "$f" ] || continue
  ln -sf "$f" "$DEST/$(basename "$f")"
  m=$((m+1))
done
echo "propios enlazados:   $m"

t=$(ls "$DEST"/*.fna 2>/dev/null | wc -l)
echo "total en el conjunto: $t"
[ "$t" -eq 900 ] || echo "AVISO: se esperaban 900 genomas"

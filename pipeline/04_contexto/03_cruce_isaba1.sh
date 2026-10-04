#!/usr/bin/env bash
# =====================================================================
# Cruza las posiciones de ISAba1 con las coordenadas de cada gen para
# determinar si el elemento se situa en la ventana de 500 pb rio arriba.
#
# Clasifica cada copia en tres categorias:
#   presente      ISAba1 se solapa con la ventana, respetando orientacion
#   ausente       la ventana esta contenida en el contig y no hay solape
#   no_evaluable  el gen esta a menos de 500 pb del extremo del contig,
#                 de modo que la ausencia no se distingue de la
#                 fragmentacion del ensamblado
#
# Entrada:  datos/contig_len_900.tsv       (paso 04/02)
#           datos/isaba1_hits900.tsv       (paso 04/02)
#           datos/<gen>_coords_900.tsv     (paso 04/01)
# Produce:  resultados/isaba1_<gen>_900.tsv
# Entorno:  cualquiera
# Duracion: segundos
# =====================================================================
set -euo pipefail
B=~/abaumannii
AWK=$B/repo/pipeline/04_contexto/03_cruce_isaba1.awk

for gen in oxa51 oxa23 oxa24 oxa58 oxa143 ndm; do
  COORD=$B/datos/${gen}_coords_900.tsv
  [ -f "$COORD" ] || { echo "  sin coordenadas: $gen"; continue; }
  OUT=$B/resultados/isaba1_${gen}_900.tsv
  awk -f "$AWK" "$B/datos/contig_len_900.tsv" "$B/datos/isaba1_hits900.tsv" "$COORD" > "$OUT"
  printf "%-10s" "$gen"
  awk -F'\t' '{t++; if($7=="presente")p++; else if($7=="ausente")a++} 
    END{printf "copias=%-6d presente=%-5d ausente=%-5d no_eval=%-5d evaluabilidad=%.1f %%\n",
        t, p, a, t-p-a, 100*(p+a)/t}' "$OUT"
done

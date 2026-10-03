#!/usr/bin/env bash
# =====================================================================
# Obtiene los tres catalogos externos que el analisis necesita.
#
# Requiere: git, amrfinder, chewBBACA.py
# Produce:  AcinetobacterPlasmidTyping-main/   esquema APT
#           datos/cgmlst_abaumannii/           esquema cgMLST, 2390 loci
#           la base de AMRFinderPlus en su ubicacion por defecto
# Duracion: 20 a 40 min
# =====================================================================
set -euo pipefail
BASE=~/abaumannii
cd "$BASE"

echo "=== 1/3  Esquema Acinetobacter Plasmid Typing ==="
if [ ! -d AcinetobacterPlasmidTyping-main ]; then
  git clone https://github.com/MehradHamidian/AcinetobacterPlasmidTyping.git \
    AcinetobacterPlasmidTyping-main
fi
ls AcinetobacterPlasmidTyping-main/*rep_DNA-seqs*.fasta

echo
echo "=== 2/3  Base de datos de AMRFinderPlus ==="
# Este trabajo empleo la version 2026-05-15.1. La orden descarga siempre
# la mas reciente; una base posterior puede notificar variantes adicionales
# y mover los recuentos de carbapenemasas.
conda run -n abaumannii amrfinder -u
conda run -n abaumannii amrfinder --database_version 2>/dev/null || true

echo
echo "=== 3/3  Esquema cgMLST (Chewie-NS) ==="
if [ ! -d "$BASE/datos/cgmlst_abaumannii" ]; then
  conda run -n ensamblaje chewBBACA.py DownloadSchema \
    -sp 5 -sc 1 -o "$BASE/datos/cgmlst_abaumannii"
fi
find "$BASE/datos/cgmlst_abaumannii" -name "*.fasta" | wc -l

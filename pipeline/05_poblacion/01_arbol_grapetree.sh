#!/usr/bin/env bash
# =====================================================================
# Red de expansion minima sobre los perfiles cgMLST.
#
# Entrada:  resultados/cgmlst_900_eval/cgMLST95.tsv
# Produce:  resultados/grapetree_900/mst_900.nwk
# Entorno:  conda activate abaumannii
# Duracion: menos de 2 min
#
# GrapeTree 2.2 presenta un fallo de importacion en su ejecutable, por lo
# que se invoca su funcion interna directamente.
# =====================================================================
set -euo pipefail
BASE=~/abaumannii
mkdir -p "$BASE/resultados/grapetree_900"

python3 -c "
from grapetree.module.MSTrees import backend
nwk = backend(profile='$BASE/resultados/cgmlst_900_eval/cgMLST95.tsv', method='MSTreeV2')
open('$BASE/resultados/grapetree_900/mst_900.nwk','w').write(nwk)
print('newick:', len(nwk), 'caracteres')
"

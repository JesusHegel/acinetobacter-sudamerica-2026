#!/usr/bin/env bash
# =====================================================================
# Red de expansion minima de los 213 ST2 y sus metadatos de anotacion.
#
# Produce los dos ficheros que consume 03_red_expansion_minima_st2.py y
# que hasta ahora no generaba ningun paso del pipeline.
#
# Entrada:  resultados/cgmlst_st2_eval/cgMLST_profiles.tsv  (paso 00b)
#           resultados/TablaS1_900_genomas.tsv
# Produce:  resultados/grapetree/st2_msa.nwk
#           resultados/grapetree/metadatos_grapetree_st2.tsv
# Entorno:  conda activate abaumannii
# Duracion: menos de 1 min
#
# GrapeTree 2.2 presenta un fallo de importacion en su ejecutable, por lo
# que se invoca su funcion interna directamente, igual que en el paso 01.
# =====================================================================
set -euo pipefail
BASE=~/abaumannii
mkdir -p "$BASE/resultados/grapetree"

echo "=== 1/2  red de expansion minima ==="
python3 -c "
from grapetree.module.MSTrees import backend
nwk = backend(profile='$BASE/resultados/cgmlst_st2_eval/cgMLST_profiles.tsv',
              method='MSTreeV2')
open('$BASE/resultados/grapetree/st2_msa.nwk','w').write(nwk)
print('    newick:', len(nwk), 'caracteres')
"

echo "=== 2/2  metadatos de anotacion ==="
python3 - <<'PY'
import os, csv, re
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")

s1 = {r["accesion"].strip(): r for r in csv.DictReader(
      open(f"{B}/resultados/TablaS1_900_genomas.tsv"), delimiter="\t")}

# Los identificadores de la matriz llevan el sufijo del ensamblado que
# chewBBACA toma de las cabeceras FASTA; se cruzan por la accesion.
hdr = open(f"{B}/resultados/cgmlst_st2_eval/distance_matrix_symmetric.tsv"
          ).readline().rstrip("\n").split("\t")[1:]

cols = ["ID","pais","anio","ST","KL","plasmido_oxa72","portador","bioproject"]
filas, sin_cruce = [], []
for ident in hdr:
    m = ACC.match(ident)
    d = s1.get(m.group(1)) if m else None
    if d is None:
        sin_cruce.append(ident); continue
    ctx = (d["oxa72_contexto"] or "").strip()
    filas.append([ident, d["pais"], d["anio"], "ST" + d["st_pasteur"], d["kl"],
                  ctx if ctx else "-", "si" if ctx else "no", d["bioproject"]])

sal = f"{B}/resultados/grapetree/metadatos_grapetree_st2.tsv"
with open(sal, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols); w.writerows(filas)

print(f"    filas: {len(filas)}   sin cruce: {len(sin_cruce)}")
if sin_cruce: print("    AVISO, no cruzan:", sin_cruce[:5])
PY

echo "listo"

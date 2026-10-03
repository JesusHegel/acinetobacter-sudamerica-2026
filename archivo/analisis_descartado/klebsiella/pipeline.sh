#!/bin/bash
set -e
cd ~/abaumannii/tmp_kleb
source ~/miniforge3/etc/profile.d/conda.sh

echo "=== [1/4] descomprimiendo ==="
mkdir -p genomas
unzip -qo datos/kleb.zip -d datos/unz
find datos/unz -name "*.fna" | while read f; do
  a=$(basename $(dirname "$f"))
  [ -f "genomas/$a.fna" ] || cp "$f" "genomas/$a.fna"
done
echo "  genomas: $(ls genomas/*.fna | wc -l)"

echo "=== [2/4] ANI con skani ==="
conda activate qc
ls ref/*.fna > ref_list.txt
ls genomas/*.fna > gen_list.txt
skani dist --ql gen_list.txt --rl ref_list.txt -t 6 --medium -o ani_raw.tsv 2>/dev/null
python3 - << 'PY'
best={}
for i,l in enumerate(open('ani_raw.tsv')):
    if i==0: continue
    f=l.split('\t')
    if len(f)<5: continue
    q=f[1].split('/')[-1].replace('.fna',''); r=f[0].split('/')[-1].replace('.fna','')
    ani=float(f[2])
    if q not in best or ani>best[q][1]: best[q]=(r,ani)
NOM={'GCF_000240185.1':'K.pneumoniae','GCF_020525545.1':'K.variicola',
     'GCF_020099175.1':'K.quasipneumoniae','GCF_000751755.1':'K.quasipneumoniae_sub',
     'GCF_020525665.1':'K.quasivariicola','GCF_019048125.1':'K.aerogenes'}
with open('ani_asignacion.tsv','w') as o:
    o.write('genoma\tespecie\tani\n')
    for q,(r,a) in sorted(best.items()):
        o.write(f'{q}\t{NOM.get(r,r)}\t{a:.2f}\n')
from collections import Counter
print('  asignacion:', dict(Counter(NOM.get(v[0],v[0]) for v in best.values())))
PY

echo "=== [3/4] MLST ==="
conda activate abaumannii
mlst --scheme klebsiella --quiet --threads 6 genomas/*.fna > mlst_kleb.tsv 2>/dev/null
echo "  con ST: $(awk -F'\t' '$3!="-"' mlst_kleb.tsv | wc -l) de $(wc -l < mlst_kleb.tsv)"

echo "=== [4/4] AMRFinderPlus ==="
mkdir -p amr
n=0
for f in genomas/*.fna; do
  a=$(basename "$f" .fna)
  [ -s "amr/$a.tsv" ] && continue
  amrfinder -n "$f" --organism Klebsiella_pneumoniae --plus -o "amr/$a.tsv" --threads 6 >/dev/null 2>&1 || true
  n=$((n+1)); [ $((n % 200)) -eq 0 ] && echo "    $n procesados"
done
echo "  ficheros AMR: $(ls amr/*.tsv | wc -l)"
echo "=== TERMINADO ==="

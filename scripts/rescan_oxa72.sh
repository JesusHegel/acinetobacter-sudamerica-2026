#!/bin/bash
set -e
B=~/abaumannii; cd $B
mkdir -p tmp_rescan
echo "[1/4] concatenando genomas con prefijo de accesion..."
: > tmp_rescan/all.fna
for d in datos/genomas_945/ncbi_dataset/data/*/; do
  a=$(basename "$d"); f=$(ls "$d"/*.fna 2>/dev/null | head -1)
  [ -n "$f" ] && awk -v A="$a" '/^>/{sub(/^>/,">"A"|");print;next}{print}' "$f" >> tmp_rescan/all.fna
done
for f in datos/ensamblados_63/*.fna; do
  a=$(basename "$f" .fna)
  awk -v A="$a" '/^>/{sub(/^>/,">"A"|");print;next}{print}' "$f" >> tmp_rescan/all.fna
done
echo "    secuencias: $(grep -c '>' tmp_rescan/all.fna)"
echo "[2/4] extrayendo referencia blaOXA-72..."
AMR=$(ls ~/miniforge3/envs/abaumannii/share/amrfinderplus/data/*/AMR_CDS.fa | head -1)
python3 - "$AMR" << 'PY'
import sys,re
k=None;buf=[];out=[]
for l in open(sys.argv[1]):
    if l.startswith(">"):
        if k and re.search(r"blaOXA-72\b",k): out.append((k,"".join(buf)))
        k=l[1:].strip();buf=[]
    else: buf.append(l.strip())
if k and re.search(r"blaOXA-72\b",k): out.append((k,"".join(buf)))
with open("tmp_rescan/oxa72.fa","w") as f:
    for h,s in out[:1]: f.write(f">blaOXA-72\n{s}\n"); print("    len:",len(s))
PY
echo "[3/4] makeblastdb..."
makeblastdb -in tmp_rescan/all.fna -dbtype nucl -out tmp_rescan/db > /dev/null
echo "[4/4] blastn..."
blastn -query tmp_rescan/oxa72.fa -db tmp_rescan/db -evalue 1e-5 -perc_identity 95 \
  -outfmt '6 sseqid pident length qlen sstart send slen' -num_threads 4 \
  > tmp_rescan/hits.tsv
echo "    hits: $(wc -l < tmp_rescan/hits.tsv)"

#!/bin/bash
cd ~/abaumannii
echo "=== N5: publicacion asociada a cada BioProject con n>=8 ==="
echo
awk -F'\t' 'NR>1{print $2"\t"$1"\t"$3}' resultados/TablaS2_proyectos.tsv | sort -u | \
while IFS=$'\t' read -r proj pais n; do
  id=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=${proj}&retmode=json" \
       | python3 -c "import json,sys; d=json.load(sys.stdin); print(d['esearchresult']['idlist'][0] if d['esearchresult']['idlist'] else '')" 2>/dev/null)
  if [ -z "$id" ]; then echo "  $proj ($pais, n=$n): no hallado"; sleep 1; continue; fi
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=bioproject&id=${id}&retmode=json" \
    | python3 -c "
import json,sys
try:
    d=json.load(sys.stdin); u=d['result'][d['result']['uids'][0]]
    t=u.get('project_title','')[:78]
    o=u.get('submitter_organization','')[:48]
    print(f'  $proj  $pais  n=$n')
    print(f'      titulo: {t}')
    print(f'      envia : {o}')
except Exception as e:
    print(f'  $proj ($pais): error')
" 2>/dev/null
  sleep 1
done

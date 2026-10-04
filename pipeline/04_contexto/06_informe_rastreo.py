#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Informe del rastreo dirigido de blaOXA-72: compara los portadores detectados
por BLAST frente a los declarados por AMRFinderPlus (seccion 2.3.3)."""
#!/usr/bin/env python3
import re,os
from collections import defaultdict,Counter
B=os.path.expanduser("~/abaumannii")
ACC=re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m=ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]
rows=tsv(f"{B}/resultados/tabla_clinica_900.tsv")
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)",rows[0][0]): rows=rows[1:]
C=["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
clin={norm(r[0]):dict(zip(C,(r+[""]*10)[:10])) for r in rows}
h=defaultdict(list)
for r in tsv(f"{B}/tmp_rescan/hits.tsv"):
    g,ctg=r[0].split("|",1)
    h[norm(g)].append((float(r[1]),int(r[2]),int(r[3]),int(r[4]),int(r[5]),int(r[6])))
det={g for g in h if g in clin}
cur={g for g,d in clin.items() if "blaOXA-72" in d["genes"]}
print(f"detectados por BLAST: {len(det)}   declarados: {len(cur)}   NUEVOS: {len(det-cur)}")
print(f"declarados sin hit BLAST: {sorted(cur-det)}\n")
print("=== genomas NUEVOS (falsos negativos) ===")
for g in sorted(det-cur):
    d=clin[g]; tot=sum(x[1] for x in h[g]); best=max(h[g],key=lambda x:x[1])
    edge=min(min(best[3],best[4])-1, best[5]-max(best[3],best[4]))
    print(f"  {g:20s} {d['pais']:10s} ST{d['st']:6s} {d['kl']:6s} "
          f"aln_tot={tot}/{best[2]} id={best[0]:.1f} borde={edge} pb  carb=[{d['carb']}] reps=[{d['rep']}]")
print(f"\n  por pais: {dict(Counter(clin[g]['pais'] for g in det-cur))}")
print(f"  por ST:   {dict(Counter('ST'+clin[g]['st'] for g in det-cur))}")
print(f"\n=== recuento corregido ===")
print(f"  portadores de blaOXA-72: {len(cur)} -> {len(det)}")
for p in sorted({clin[g]['pais'] for g in det}):
    a=sum(1 for g in cur if clin[g]['pais']==p); b=sum(1 for g in det if clin[g]['pais']==p)
    print(f"    {p:12s} {a:3d} -> {b:3d}")

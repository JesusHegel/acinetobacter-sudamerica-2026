#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Contenido de replicones por contig, diametros de grupo, conflicto Pons."""
import re, os
from collections import Counter, defaultdict
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p):
    return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]

rows = tsv(f"{B}/resultados/tabla_clinica_900.tsv")
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)", rows[0][0]): rows = rows[1:]
C = ["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
clin = {norm(r[0]): dict(zip(C,(r+[""]*10)[:10])) for r in rows}

cr = tsv(f"{B}/resultados/apt_v2/colocalizacion.tsv")
ci = {n:i for i,n in enumerate(cr[0])}
O = [r for r in cr[1:] if r[ci["gen"]]=="blaOXA-72"]

print("=== A. genomas con blaOXA-72 AUSENTES de la tabla clinica ===")
for r in O:
    g = norm(r[ci["genoma"]])
    if g not in clin: print(f"  {g}  estado={r[ci['estado']]}  contig={r[ci['gen_contig']]}")

COL = [r for r in O if r[ci["reps_mismo_contig"]] not in ("-","","NA") and norm(r[ci["genoma"]]) in clin]
print(f"\n=== B. contigs co-localizados clinicos: {len(COL)} filas, "
      f"{len(set(norm(r[ci['genoma']]) for r in COL))} genomas ===")

AND = {"Peru","Perú","Ecuador"}
print("\n=== C. ANDINOS CC2: replicones EN EL MISMO CONTIG ===")
for r in sorted(COL, key=lambda x: norm(x[ci["genoma"]])):
    g = norm(r[ci["genoma"]]); d = clin[g]
    if d["pais"] in AND and d["st"] in ("2","2724"):
        print(f"  {g:20s} {d['pais']:8s} {d['kl']:5s} len={r[ci['contig_len']]:>7s} "
              f"MISMO=[{r[ci['reps_mismo_contig']]}]  otros=[{r[ci['reps_otros_contigs']]}]")

print("\n=== D. reps_mismo_contig, recuento por region ===")
agg = defaultdict(Counter)
for r in COL:
    d = clin[norm(r[ci["genoma"]])]
    reg = "ANDES-CC2" if (d["pais"] in AND and d["st"] in ("2","2724")) else d["pais"]
    agg[reg][r[ci["reps_mismo_contig"]]] += 1
for reg in sorted(agg):
    print(f"  {reg}:")
    for k,v in agg[reg].most_common(): print(f"      {v:3d} x [{k}]")

print("\n=== E. contigs duplicados por genoma (regla 17) ===")
cnt = Counter(norm(r[ci["genoma"]]) for r in COL)
for g,n in cnt.most_common():
    if n>1: print(f"  {g}  {n} contigs  ({clin[g]['pais']})")

print("\n=== F. hits APT crudos de los dos genomas de Pons (punto 11) ===")
H = f"{B}/resultados/apt_v2/apt_hits_v2.tsv"
hr = tsv(H); print("  cabecera:", "\t".join(hr[0]))
for l in hr[1:]:
    if any(a in l[0] for a in ("GCA_051941745","GCA_051942065")): print("  "+"\t".join(l))

# --- diametros ---
with open(f"{B}/resultados/cgmlst_st2_eval/distance_matrix_symmetric.tsv") as f:
    names = f.readline().rstrip("\n").split("\t")[1:]
    D = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in f}
idx = {n:i for i,n in enumerate(names)}; m2n = {norm(n):n for n in names}
def dist(a,b): return int(float(D[a][idx[b]]))
cc2 = sorted({norm(r[ci["genoma"]]) for r in COL
              if clin[norm(r[ci["genoma"]])]["pais"] in AND
              and clin[norm(r[ci["genoma"]])]["st"] in ("2","2724")} & set(m2n))
N = [m2n[g] for g in cc2]; acc = {m2n[g]:g for g in cc2}
def clusters(nodes,t):
    par={n:n for n in nodes}
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i,a in enumerate(nodes):
        for b in nodes[i+1:]:
            if dist(a,b)<=t: par[f(a)]=f(b)
    g=defaultdict(list)
    for n in nodes: g[f(n)].append(n)
    return sorted(g.values(), key=len, reverse=True)
print("\n=== G. diametro de cada grupo (encadenamiento por enlace simple) ===")
for t in (8,10,12,25):
    print(f"  umbral {t}: {len(clusters(N,t))} grupos")
    for i,gr in enumerate(clusters(N,t),1):
        dia = max((dist(a,b) for a in gr for b in gr if a!=b), default=0)
        kls = sorted({clin[acc[n]]['kl'] for n in gr})
        pj  = len({clin[acc[n]]['anio'] for n in gr})
        print(f"     grupo {i}: n={len(gr):2d}  diametro={dia:3d}  KL={','.join(kls)}  anios={pj}")

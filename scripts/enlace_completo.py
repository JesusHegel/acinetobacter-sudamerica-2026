#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V42: enlace simple vs completo en el conjunto andino CC2 + hueco KL2/KL9."""
import re, os, sys
from collections import defaultdict

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
AND = {"Peru","Perú","Ecuador"}
cc2 = sorted({norm(r[ci["genoma"]]) for r in cr[1:]
              if r[ci["gen"]] == "blaOXA-72"
              and r[ci["reps_mismo_contig"]] not in ("-","","NA")
              and norm(r[ci["genoma"]]) in clin
              and clin[norm(r[ci["genoma"]])]["pais"] in AND
              and clin[norm(r[ci["genoma"]])]["st"] in ("2","2724")})

MAT = f"{B}/resultados/cgmlst_st2_eval/distance_matrix_symmetric.tsv"
with open(MAT) as f:
    names = f.readline().rstrip("\n").split("\t")[1:]
    D = {}
    for l in f:
        p = l.rstrip("\n").split("\t"); D[p[0]] = p[1:]
idx = {n:i for i,n in enumerate(names)}
m2n = {norm(n):n for n in names}
def dist(a,b): return int(float(D[a][idx[b]]))

N   = [m2n[g] for g in cc2 if g in m2n]
acc = {m2n[g]:g for g in cc2 if g in m2n}
print(f"subconjunto CC2 con blaOXA-72 co-localizado: {len(cc2)} genomas")
print(f"evaluados en la matriz: {len(N)}   fuera: {[g for g in cc2 if g not in m2n]}\n")
if len(N) < 2: sys.exit("ERROR: menos de 2 genomas evaluables")

def single(nodes, t):
    par = {n:n for n in nodes}
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i,a in enumerate(nodes):
        for b in nodes[i+1:]:
            if dist(a,b) <= t: par[f(a)] = f(b)
    g = defaultdict(list)
    for n in nodes: g[f(n)].append(n)
    return sorted(g.values(), key=len, reverse=True)

def complete(nodes, t):
    cl = [[n] for n in nodes]
    while True:
        best = None; bd = None
        for i in range(len(cl)):
            for j in range(i+1, len(cl)):
                d = max(dist(a,b) for a in cl[i] for b in cl[j])
                if d <= t and (bd is None or d < bd): bd = d; best = (i,j)
        if best is None: break
        i,j = best; cl[i] += cl[j]; del cl[j]
    return sorted(cl, key=len, reverse=True)

print("=== A. numero de grupos: enlace SIMPLE vs COMPLETO ===")
print("umbral  simple  completo")
for t in range(0,26):
    print(f"{t:5d}  {len(single(N,t)):6d}  {len(complete(N,t)):8d}")

for t in (5,8,10):
    print(f"\n=== B. enlace COMPLETO a umbral {t} ===")
    for i,gr in enumerate(complete(N,t),1):
        dia = max((dist(a,b) for a in gr for b in gr if a!=b), default=0)
        print(f"  grupo {i}: n={len(gr)}  diametro={dia}")
        for n in sorted(gr, key=lambda x: acc[x]):
            d = clin[acc[n]]
            print(f"     {acc[n]:20s} {d['pais']:8s} {d['anio']:5s} {d['kl']}")

print("\n=== C. hueco de la biparticion KL2 / KL9 ===")
kl2 = [n for n in N if clin[acc[n]]["kl"]=="KL2"]
kl9 = [n for n in N if clin[acc[n]]["kl"]=="KL9"]
d2 = [dist(a,b) for i,a in enumerate(kl2) for b in kl2[i+1:]]
d9 = [dist(a,b) for i,a in enumerate(kl9) for b in kl9[i+1:]]
en = [dist(a,b) for a in kl2 for b in kl9]
print(f"  KL2 n={len(kl2)}: intra min={min(d2)} max={max(d2)}")
print(f"  KL9 n={len(kl9)}: intra min={min(d9)} max={max(d9)}")
print(f"  entre grupos:      min={min(en)} max={max(en)}")
print(f"  HUECO: intra maxima={max(max(d2),max(d9))}  ->  inter minima={min(en)}")

print("\n=== D. pares que exceden el umbral dentro de cada grupo (simple, t=8) ===")
for i,gr in enumerate(single(N,8),1):
    ex = [(acc[a],acc[b],dist(a,b)) for j,a in enumerate(gr) for b in gr[j+1:] if dist(a,b)>8]
    print(f"  grupo {i} (n={len(gr)}): {len(ex)} pares >8")
    for a,b,d in sorted(ex, key=lambda x:-x[2])[:6]:
        print(f"     {a} - {b}  d={d}")

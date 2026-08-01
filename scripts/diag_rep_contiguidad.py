#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Punto 43: la deteccion de genes rep depende de la contiguidad del ensamblado?"""
import re, os
from collections import Counter, defaultdict
from statistics import median
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]

rows = tsv(f"{B}/resultados/tabla_clinica_900.tsv")
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)", rows[0][0]): rows = rows[1:]
C = ["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
clin = {norm(r[0]): dict(zip(C,(r+[""]*10)[:10])) for r in rows}
bp = {norm(r[0]): r[1] for r in tsv(f"{B}/datos/acc_bioproject_944.tsv")
      if len(r)>=2 and r[1].startswith("PRJ")}

ncont = Counter()
for r in tsv(f"{B}/datos/contig_len.tsv"):
    if len(r)>=2: ncont[norm(r[0])] += 1
print(f"contig_len.tsv: {len(ncont)} genomas con recuento de contigs")

hr = tsv(f"{B}/resultados/apt_v2/apt_hits_v2.tsv"); hi={n:i for i,n in enumerate(hr[0])}
reps = defaultdict(set); edge = []
for r in hr[1:]:
    g = norm(r[hi["genoma"]]); reps[g].add(r[hi["tipo_rep"]])
    try:
        ini,fin,cl = int(r[hi["rep_ini"]]), int(r[hi["rep_fin"]]), int(r[hi["contig_len"]])
        edge.append(min(min(ini,fin)-1, cl-max(ini,fin)))
    except (ValueError, KeyError): pass

ev = [g for g in clin if g in ncont]
print(f"genomas clinicos con ambos datos: {len(ev)}\n")

print("=== A. riqueza de tipos rep por cuartil de n contigs ===")
srt = sorted(ev, key=lambda g: ncont[g]); q = len(srt)//4
for i,(lo,hi_) in enumerate([(0,q),(q,2*q),(2*q,3*q),(3*q,len(srt))],1):
    sub = srt[lo:hi_]
    nr = [len(reps.get(g,())) for g in sub]
    print(f"  Q{i}: n={len(sub):3d}  contigs {ncont[sub[0]]:4d}-{ncont[sub[-1]]:4d}  "
          f"mediana tipos rep={median(nr):4.1f}  sin ningun rep={sum(1 for x in nr if x==0):3d} "
          f"({100*sum(1 for x in nr if x==0)/len(sub):.1f} %)")

print("\n=== B. correlacion de rangos (Spearman, calculo directo) ===")
def rank(v):
    s=sorted(range(len(v)), key=lambda i:v[i]); r=[0]*len(v); i=0
    while i<len(s):
        j=i
        while j+1<len(s) and v[s[j+1]]==v[s[i]]: j+=1
        avg=(i+j)/2+1
        for k in range(i,j+1): r[s[k]]=avg
        i=j+1
    return r
x=[ncont[g] for g in ev]; y=[len(reps.get(g,())) for g in ev]
rx,ry=rank(x),rank(y); n=len(x)
mx,my=sum(rx)/n,sum(ry)/n
num=sum((a-mx)*(b-my) for a,b in zip(rx,ry))
den=(sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry))**.5
print(f"  rho = {num/den:+.3f}   (n={n})")

print("\n=== C. proximidad de los hits al borde del contig ===")
edge.sort()
print(f"  hits con coordenadas validas: {len(edge)}")
for p in (0,1,5,10,25,50):
    print(f"    p{p:<2d} = {edge[min(int(len(edge)*p/100), len(edge)-1)]:>8d} pb del borde")
print(f"  hits a <100 pb del borde: {sum(1 for e in edge if e<100)} "
      f"({100*sum(1 for e in edge if e<100)/len(edge):.2f} %)")

print("\n=== D. por BioProject con n>=8: contigs vs deteccion de rep ===")
byp = defaultdict(list)
for g in ev: byp[bp.get(g,"?")].append(g)
print(f"  {'BioProject':16s} {'pais':10s} {'n':>4s} {'med.contigs':>12s} {'med.reps':>9s} {'%sin rep':>9s} {'dominante':>12s}")
for p,gs in sorted(byp.items(), key=lambda kv: -len(kv[1])):
    if len(gs)<8 or p=="?": continue
    nr=[len(reps.get(g,())) for g in gs]
    dom=Counter(t for g in gs for t in reps.get(g,())).most_common(1)
    print(f"  {p:16s} {clin[gs[0]]['pais']:10s} {len(gs):4d} {median([ncont[g] for g in gs]):12.1f} "
          f"{median(nr):9.1f} {100*sum(1 for x in nr if x==0)/len(gs):8.1f}% "
          f"{(dom[0][0] if dom else '-'):>12s}")

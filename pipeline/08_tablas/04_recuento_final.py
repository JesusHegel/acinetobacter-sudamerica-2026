#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recuento definitivo con los 73 portadores."""
import re, os
from collections import defaultdict, Counter
B=os.path.expanduser("~/abaumannii")
ACC=re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m=ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]
rows=tsv(f"{B}/resultados/tabla_clinica_900.tsv")
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)",rows[0][0]): rows=rows[1:]
C=["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
clin={norm(r[0]):dict(zip(C,(r+[""]*10)[:10])) for r in rows}
bp={norm(r[0]):r[1] for r in tsv(f"{B}/datos/acc_bioproject_944.tsv")
    if len(r)>=2 and r[1].startswith("PRJ")}

o72=defaultdict(set)
for r in tsv(f"{B}/tmp_rescan/hits.tsv"):
    g,ctg=r[0].split("|",1)
    if norm(g) in clin: o72[norm(g)].add(ctg)
apt=defaultdict(lambda: defaultdict(set)); clen={}
for r in tsv(f"{B}/resultados/apt_v2/apt_hits_v2.tsv")[1:]:
    apt[norm(r[0])][r[2]].add(r[1]); clen[(norm(r[0]),r[2])]=r[5]

print(f"=== A. co-localizacion recalculada (por genoma) ===")
col={}; noev=[]
for g,ctgs in o72.items():
    hit={c:apt[g][c] for c in ctgs if apt[g].get(c)}
    if hit: col[g]=hit
    else: noev.append(g)
print(f"  portadores={len(o72)}  co-localizados={len(col)} ({100*len(col)/len(o72):.1f} %)  "
      f"no evaluables={len(noev)}   [antes 54/64 = 84,4 %]")
print(f"  por pais: {dict(Counter(clin[g]['pais'] for g in col))}")

AND={"Peru","Perú","Ecuador"}
cc2=sorted(g for g in col if clin[g]['pais'] in AND and clin[g]['st'] in ("2","2724"))
print(f"\n=== B. CC2 andino: {len(cc2)} genomas [antes 19] ===")
for g in cc2:
    d=clin[g]; reps=sorted({t for s in col[g].values() for t in s})
    L=[clen.get((g,c),'?') for c in col[g]]
    print(f"  {g:20s} {d['pais']:8s} ST{d['st']:5s} {d['anio']:5s} {d['kl']:6s} "
          f"{bp.get(g,'?'):14s} {','.join(reps):12s} len={','.join(L)}")
print(f"  replicones: {dict(Counter(t for g in cc2 for s in col[g].values() for t in s))}")
br=[g for g in col if clin[g]['pais']=='Brasil']
print(f"  Brasil n={len(br)} replicones: "
      f"{dict(Counter(','.join(sorted(s)) for g in br for s in col[g].values()))}")

with open(f"{B}/resultados/cgmlst_st2_eval/distance_matrix_symmetric.tsv") as f:
    names=f.readline().rstrip("\n").split("\t")[1:]
    D={l.split("\t")[0]:l.rstrip("\n").split("\t")[1:] for l in f}
idx={n:i for i,n in enumerate(names)}; m2n={norm(n):n for n in names}
def dist(a,b): return int(float(D[a][idx[b]]))
N=[m2n[g] for g in cc2 if g in m2n]; acc={m2n[g]:g for g in cc2 if g in m2n}
print(f"\n=== C. agrupamiento: {len(N)} en matriz, fuera={[g for g in cc2 if g not in m2n]} ===")
def single(nd,t):
    par={n:n for n in nd}
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i,a in enumerate(nd):
        for b in nd[i+1:]:
            if dist(a,b)<=t: par[f(a)]=f(b)
    g=defaultdict(list)
    for n in nd: g[f(n)].append(n)
    return sorted(g.values(),key=len,reverse=True)
def comp(nd,t):
    cl=[[n] for n in nd]
    while True:
        bst=None;bd=None
        for i in range(len(cl)):
            for j in range(i+1,len(cl)):
                d=max(dist(a,b) for a in cl[i] for b in cl[j])
                if d<=t and (bd is None or d<bd): bd=d;bst=(i,j)
        if bst is None: break
        i,j=bst; cl[i]+=cl[j]; del cl[j]
    return sorted(cl,key=len,reverse=True)
print("  umbral  simple  completo")
for t in range(0,26): print(f"  {t:5d}  {len(single(N,t)):6d}  {len(comp(N,t)):8d}")
print(f"\n  --- enlace COMPLETO, umbral 8: {len(comp(N,8))} eventos ---")
for i,gr in enumerate(comp(N,8),1):
    dia=max((dist(a,b) for a in gr for b in gr if a!=b),default=0)
    print(f"   grupo {i}: n={len(gr)} diam={dia} "
          f"KL={sorted({clin[acc[n]]['kl'] for n in gr})} "
          f"paises={sorted({clin[acc[n]]['pais'] for n in gr})}")
    for n in sorted(gr,key=lambda x:acc[x]): print(f"      {acc[n]} {clin[acc[n]]['anio']}")

print("\n=== D. hueco KL2/KL9 ===")
k2=[n for n in N if clin[acc[n]]['kl']=='KL2']; k9=[n for n in N if clin[acc[n]]['kl']=='KL9']
o=[n for n in N if clin[acc[n]]['kl'] not in ('KL2','KL9')]
print(f"  KL2 n={len(k2)}  KL9 n={len(k9)}  otros={[(acc[n],clin[acc[n]]['kl']) for n in o]}")
if k2 and k9:
    d2=[dist(a,b) for i,a in enumerate(k2) for b in k2[i+1:]]
    d9=[dist(a,b) for i,a in enumerate(k9) for b in k9[i+1:]]
    en=[dist(a,b) for a in k2 for b in k9]
    print(f"  KL2 intra 0-{max(d2)}   KL9 intra 0-{max(d9)}   inter {min(en)}-{max(en)}")
    print(f"  HUECO: {max(max(d2),max(d9))} -> {min(en)}")

print("\n=== E. Tabla 2 corregida (regla del representante) ===")
for g in o72:
    if "blaOXA-72" not in clin[g]["genes"]:
        clin[g]["genes"]=(clin[g]["genes"].replace("ninguna","").strip(",")+",blaOXA-72").strip(",")
        clin[g]["carb"]="si"
ev=defaultdict(list)
for g,d in clin.items(): ev[(d["pais"],d["st"],bp.get(g,"?"),d["anio"])].append(g)
print(f"  {'pais':11s} {'ev':>4s} {'OXA-23':>8s} {'OXA-72':>8s} {'OXA-58':>8s} {'NDM':>8s} {'sinCarb':>8s}")
for p,_ in Counter(d["pais"] for d in clin.values()).most_common():
    ks=[k for k in ev if k[0]==p]; n=len(ks)
    f=lambda t: 100*sum(1 for k in ks if t in clin[sorted(ev[k])[0]]["genes"])/n
    sc=100*sum(1 for k in ks if clin[sorted(ev[k])[0]]["carb"].lower() in ("no","-",""))/n
    print(f"  {p:11s} {n:4d} {f('blaOXA-23'):7.1f}% {f('blaOXA-72'):7.1f}% "
          f"{f('blaOXA-58'):7.1f}% {f('blaNDM'):7.1f}% {sc:7.1f}%")
print("\n  sin desduplicar:")
for p in ("Peru","Ecuador","Brasil"):
    gs=[g for g,d in clin.items() if d["pais"]==p]
    k=sum(1 for g in gs if "blaOXA-72" in clin[g]["genes"])
    print(f"    {p}: {k}/{len(gs)} = {100*k/len(gs):.1f} %")

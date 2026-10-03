#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura 5: dendrograma del conjunto andino ST2 y distribucion de distancias."""
import os, re
from collections import defaultdict
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m=ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]

S1=tsv(f"{B}/resultados/TablaS1_900_genomas.tsv"); h=S1[0]
clin={r[0]:dict(zip(h,r+[""]*len(h))) for r in S1[1:]}
with open(f"{B}/resultados/cgmlst_st2_eval/distance_matrix_symmetric.tsv") as f:
    names=f.readline().rstrip("\n").split("\t")[1:]
    D={l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in f}
idx={n:i for i,n in enumerate(names)}; m2n={norm(n):n for n in names}
def dist(a,b): return int(float(D[a][idx[b]]))

AND={"Peru","Perú","Ecuador"}
port={g for g,d in clin.items() if d["oxa72_contexto"] not in ("","no_evaluable")
      and d["pais"] in AND and d["st_pasteur"] in ("2","2724") and g in m2n}
N=sorted(port, key=lambda g:(clin[g]["kl"], clin[g]["pais"], g))
print(f"genomas andinos portadores en la matriz: {len(N)}")

# UPGMA
cl=[[g] for g in N]; hist=[]
def dmed(a,b): return sum(dist(m2n[x],m2n[y]) for x in a for y in b)/(len(a)*len(b))
pos={g:i for i,g in enumerate(N)}; xpos={tuple(c):pos[c[0]] for c in cl}; hgt={}
while len(cl)>1:
    best=None; bd=None
    for i in range(len(cl)):
        for j in range(i+1,len(cl)):
            d=dmed(cl[i],cl[j])
            if bd is None or d<bd: bd=d; best=(i,j)
    i,j=best; a,b=cl[i],cl[j]
    new=a+b
    hist.append((tuple(a),tuple(b),bd))
    xpos[tuple(new)]=(xpos[tuple(a)]+xpos[tuple(b)])/2
    hgt[tuple(new)]=bd
    cl=[c for k,c in enumerate(cl) if k not in (i,j)]+[new]

mx=max(d for _,_,d in hist)
W,H=980,620; X0,Y0=170,60; PW=560; PH=330
sc=lambda d: Y0+PH-(d/mx)*PH
KLC={"KL2":"#2166ac","KL9":"#b2182b"}
s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
   f'<rect width="{W}" height="{H}" fill="white"/>',
   f'<text x="20" y="26" font-family="Helvetica,Arial" font-size="14" font-weight="700">A</text>',
   f'<text x="40" y="26" font-family="Helvetica,Arial" font-size="12">Agrupamiento de los {len(N)} genomas andinos portadores (UPGMA sobre distancias alelicas cgMLST)</text>']
step=PW/max(len(N)-1,1)
gx=lambda g: X0+pos[g]*step
for a,b,d in hist:
    xa=X0+xpos[a]*step; xb=X0+xpos[b]*step; y=sc(d)
    ya=sc(hgt.get(a,0)) if a in hgt else Y0+PH
    yb=sc(hgt.get(b,0)) if b in hgt else Y0+PH
    s.append(f'<path d="M{xa:.1f},{ya:.1f} L{xa:.1f},{y:.1f} L{xb:.1f},{y:.1f} L{xb:.1f},{yb:.1f}" fill="none" stroke="#444" stroke-width="1.1"/>')
for g in N:
    d=clin[g]; x=gx(g)
    s.append(f'<circle cx="{x:.1f}" cy="{Y0+PH}" r="4" fill="{KLC.get(d["kl"],"#888")}"/>')
    s.append(f'<text x="{x:.1f}" y="{Y0+PH+16}" text-anchor="end" font-family="Helvetica,Arial" font-size="7.5" '
             f'transform="rotate(-70 {x:.1f} {Y0+PH+16})">{g} {d["pais"][:2]} {d["anio"]}</text>')
for v,lb in [(8,"8 alelos (rango de brote)"),(21,"21 (max. intragrupo)"),(30,"30 (min. intergrupo)")]:
    if v<=mx:
        y=sc(v)
        s.append(f'<line x1="{X0-12}" y1="{y:.1f}" x2="{X0+PW+10}" y2="{y:.1f}" stroke="#bbb" stroke-dasharray="4,3" stroke-width="1"/>')
        s.append(f'<text x="{X0+PW+14}" y="{y+3:.1f}" font-family="Helvetica,Arial" font-size="8.5" fill="#666">{lb}</text>')
s.append(f'<line x1="{X0-12}" y1="{Y0}" x2="{X0-12}" y2="{Y0+PH}" stroke="#333"/>')
for v in [0,10,20,30,40,50]:
    if v<=mx:
        y=sc(v)
        s.append(f'<line x1="{X0-17}" y1="{y:.1f}" x2="{X0-12}" y2="{y:.1f}" stroke="#333"/>')
        s.append(f'<text x="{X0-21}" y="{y+3:.1f}" text-anchor="end" font-family="Helvetica,Arial" font-size="9">{v}</text>')
s.append(f'<text x="{X0-52}" y="{Y0+PH/2}" text-anchor="middle" font-family="Helvetica,Arial" font-size="10" '
         f'transform="rotate(-90 {X0-52} {Y0+PH/2})">Diferencias alelicas</text>')
lx=X0
for kl,c in KLC.items():
    s.append(f'<circle cx="{lx}" cy="{Y0-14}" r="4" fill="{c}"/>')
    s.append(f'<text x="{lx+9}" y="{Y0-10}" font-family="Helvetica,Arial" font-size="9.5">Tipo capsular {kl}</text>')
    lx+=125

y2=520
s.append(f'<text x="20" y="{y2-40}" font-family="Helvetica,Arial" font-size="14" font-weight="700">B</text>')
s.append(f'<text x="40" y="{y2-40}" font-family="Helvetica,Arial" font-size="12">Distancias por pares (mediana y rango)</text>')
kl2=[g for g in N if clin[g]["kl"]=="KL2"]; kl9=[g for g in N if clin[g]["kl"]=="KL9"]
AN2={g for g,d in clin.items() if d["pais"] in AND and d["st_pasteur"] in ("2","2724") and g in m2n}
otros=[n for n in names if norm(n) not in AN2]
def stat(a,b,mismo=False):
    v=[dist(m2n[x] if x in m2n else x, m2n[y] if y in m2n else y)
       for i,x in enumerate(a) for y in (b[i+1:] if mismo else b)]
    return (min(v), sorted(v)[len(v)//2], max(v), len(v)) if v else (0,0,0,0)
FIL=[("KL2 entre si",stat(kl2,kl2,True),"#2166ac"),
     ("KL9 entre si",stat(kl9,kl9,True),"#b2182b"),
     ("KL2 frente a KL9",stat(kl2,kl9),"#7b3294"),
     ("Andinos frente a ST2 no andinos",stat([m2n[g] for g in N],otros),"#666")]
MX=max(f[1][2] for f in FIL); BX=280; BW=520
for i,(lb,(lo,md,hi,np_),c) in enumerate(FIL):
    y=y2+i*24
    s.append(f'<text x="{BX-10}" y="{y+4}" text-anchor="end" font-family="Helvetica,Arial" font-size="10">{lb}</text>')
    x1=BX+BW*lo/MX; x2=BX+BW*hi/MX; xm=BX+BW*md/MX
    s.append(f'<line x1="{x1:.1f}" y1="{y}" x2="{x2:.1f}" y2="{y}" stroke="{c}" stroke-width="3"/>')
    s.append(f'<circle cx="{xm:.1f}" cy="{y}" r="4.5" fill="white" stroke="{c}" stroke-width="2"/>')
    s.append(f'<text x="{BX+BW+10}" y="{y+4}" font-family="Helvetica,Arial" font-size="9" fill="#555">{lo}-{hi}, mediana {md}</text>')
    print(f"  {lb:34s} min={lo:3d} mediana={md:3d} max={hi:3d}  ({np_} pares)")
s.append('</svg>')
open(f"{B}/resultados/figuras/Fig5_cgmlst.svg","w").write("\n".join(s))
print(f"\nescrita: resultados/figuras/Fig5_cgmlst.svg")

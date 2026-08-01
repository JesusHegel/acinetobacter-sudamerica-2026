#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cuatro verificaciones pendientes."""
import os, glob, subprocess
from collections import defaultdict, Counter
B = os.path.expanduser("~/abaumannii")
def fasta(p):
    d={};k=None;buf=[]
    for l in open(p):
        if l.startswith(">"):
            if k: d[k]="".join(buf)
            k=l[1:].split()[0]; buf=[]
        else: buf.append(l.strip())
    if k: d[k]="".join(buf)
    return d
def fpath(a):
    g=glob.glob(f"{B}/datos/genomas_945/ncbi_dataset/data/{a}/*.fna")
    if g: return g[0]
    p=f"{B}/datos/ensamblados_63/{a}.fna"
    return p if os.path.exists(p) else None
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p,encoding="utf-8",errors="replace") if l.strip()]

print("="*66)
print("V1. Que contienen los 253-305 pb extra del contig andino de 8111 pb")
print("="*66)
A=fasta(fpath("GCA_051942065.1"))["JBFECW010000069.1"]
Bs=fasta(fpath("GCA_053320315.1"))["JBROGJ010000074.1"]
os.makedirs(f"{B}/tmp_ver", exist_ok=True)
open(f"{B}/tmp_ver/a.fa","w").write(f">A_8111\n{A}\n")
open(f"{B}/tmp_ver/b.fa","w").write(f">B_7851\n{Bs}\n")
r=subprocess.run(["blastn","-query",f"{B}/tmp_ver/a.fa","-subject",f"{B}/tmp_ver/b.fa",
                  "-evalue","1e-5","-outfmt","6 qstart qend sstart send pident length"],
                 capture_output=True,text=True)
iv=[]
for l in r.stdout.strip().split("\n"):
    if l:
        f=l.split("\t"); iv.append((min(int(f[0]),int(f[1])),max(int(f[0]),int(f[1]))))
        print(f"  bloque A:{f[0]}-{f[1]}  B:{f[2]}-{f[3]}  id={f[4]} len={f[5]}")
iv.sort(); cov=[]; last=0
for x,y in iv:
    if x>last+1: cov.append((last+1,x-1))
    last=max(last,y)
if last<len(A): cov.append((last+1,len(A)))
print(f"\n  regiones de A sin correspondencia en B:")
for x,y in cov:
    if y-x >= 30: print(f"    {x}-{y}  ({y-x+1} pb)")

print("\n"+"="*66)
print("V2. Los ST108 ecuatorianos: donde estan r3-T13 y r3-T50")
print("="*66)
apt=defaultdict(list)
for r in tsv(f"{B}/resultados/apt_v2/apt_hits_v2.tsv")[1:]: apt[r[0]].append(r)
cr=tsv(f"{B}/resultados/apt_v2/colocalizacion.tsv"); ci={n:i for i,n in enumerate(cr[0])}
for a in ("GCA_029955865.1","GCA_029956025.1"):
    print(f"\n  --- {a} ---")
    for r in apt.get(a,[]):
        print(f"    {r[1]:10s} contig={r[2]:26s} len={r[5]:>7s} id={r[6]} cov={r[7]}")
    for r in cr[1:]:
        if r[ci['genoma']]==a and r[ci['gen']]=="blaOXA-72":
            print(f"    blaOXA-72 -> contig {r[ci['gen_contig']]}  estado={r[ci['estado']]}")

print("\n"+"="*66)
print("V3. Contigs frente a genomas en los co-localizados")
print("="*66)
S1=tsv(f"{B}/resultados/TablaS1_900_genomas.tsv"); h=S1[0]
clin={r[0]:dict(zip(h,r+[""]*len(h))) for r in S1[1:]}
hits=defaultdict(set)
for r in tsv(f"{B}/tmp_rescan/hits.tsv"):
    g,c=r[0].split("|",1); hits[g.split("|")[0]].add(c)
import re
ACC=re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m=ACC.search(x); return m.group(1) if m else x.strip()
aptc=defaultdict(lambda: defaultdict(set))
for r in tsv(f"{B}/resultados/apt_v2/apt_hits_v2.tsv")[1:]: aptc[norm(r[0])][r[2]].add(r[1])
ng=nc=0; dup=[]
por_pais=Counter(); ctg_pais=Counter()
for g,cs in hits.items():
    g=norm(g)
    if g not in clin: continue
    val=[c for c in cs if aptc[g].get(c)]
    if val:
        ng+=1; nc+=len(val); por_pais[clin[g]["pais"]]+=1; ctg_pais[clin[g]["pais"]]+=len(val)
        if len(val)>1: dup.append((g,clin[g]["pais"],len(val)))
print(f"  genomas co-localizados: {ng}   contigs: {nc}")
print(f"  genomas por pais : {dict(por_pais)}")
print(f"  contigs por pais : {dict(ctg_pais)}")
print(f"\n  genomas con mas de un contig ({len(dup)}):")
for g,p,n in sorted(dup, key=lambda x:x[1]): print(f"    {g:20s} {p:10s} {n} contigs")

print("\n"+"="*66)
print("V4. Evaluabilidad real de la deteccion de genes de resistencia")
print("="*66)
print(f"  portadores de blaOXA-72 por AMRFinderPlus estandar : 64")
print(f"  portadores tras el rastreo dirigido                : 73")
print(f"  falsos negativos del procedimiento estandar        : 9  ({100*9/73:.1f} % de los portadores)")
print(f"\n  La Tabla 4 no puede afirmar 100 % de evaluabilidad para genes de resistencia:")
print(f"  la deteccion produjo un resultado en los 900 genomas, pero ese resultado")
print(f"  fue incorrecto en 9 casos. Evaluabilidad y exactitud son cosas distintas.")

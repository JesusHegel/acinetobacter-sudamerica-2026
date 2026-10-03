#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Punto 21: BioProjects multipais, independencia de las ausencias replicadas."""
import re, os, glob
from collections import Counter, defaultdict
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
sinbp = [g for g in clin if g not in bp]
print(f"clinicos: {len(clin)}   sin BioProject asignado: {len(sinbp)} {sinbp[:5]}")

byp = defaultdict(Counter)
for g in clin:
    if g in bp: byp[bp[g]][clin[g]["pais"]] += 1
print(f"BioProjects unicos en el conjunto clinico: {len(byp)}   (V37 reporto 104)")
suma = sum(len(v) for v in byp.values())
print(f"suma de pares proyecto-pais: {suma}   (Tabla 1 suma 114)")

multi = {p:c for p,c in byp.items() if len(c)>1}
print(f"\n=== A. BioProjects MULTIPAIS: {len(multi)} ===")
for p,c in sorted(multi.items(), key=lambda kv:-sum(kv[1].values())):
    print(f"  {p:16s} n={sum(c.values()):4d}  {dict(c)}")

print("\n=== B. contraste DENTRO de proyectos multipais (proyecto constante) ===")
for p,c in sorted(multi.items(), key=lambda kv:-sum(kv[1].values())):
    if sum(c.values())<8: continue
    print(f"\n  {p}  (n={sum(c.values())})")
    print(f"    {'pais':12s} {'n':>4s} {'OXA-72':>8s} {'OXA-23':>8s} {'ST2':>7s} {'sin carb':>9s}")
    for pais in sorted(c):
        gs=[g for g in clin if bp.get(g)==p and clin[g]["pais"]==pais]
        n=len(gs)
        f=lambda k: 100*sum(1 for g in gs if k in clin[g]["genes"])/n
        st2=100*sum(1 for g in gs if clin[g]["st"] in ("2","2724"))/n
        sc=100*sum(1 for g in gs if clin[g]["carb"].lower() in ("no","-",""))/n
        print(f"    {pais:12s} {n:4d} {f('blaOXA-72'):7.1f}% {f('blaOXA-23'):7.1f}% "
              f"{st2:6.1f}% {sc:8.1f}%")

print("\n=== C. proyectos que sostienen las AUSENCIAS REPLICADAS ===")
for var,paises in (("blaOXA-72",["Argentina","Chile","Paraguay"]),
                   ("ST2",["Chile","Colombia","Bolivia"])):
    print(f"\n  --- ausencia de {var} ---")
    for pais in paises:
        gs=[g for g in clin if clin[g]["pais"]==pais]
        pr=defaultdict(list)
        for g in gs: pr[bp.get(g,"?")].append(g)
        for p,l in sorted(pr.items(), key=lambda kv:-len(kv[1])):
            if len(l)<8: continue
            pos=sum(1 for g in l if (var=="ST2" and clin[g]["st"] in ("2","2724"))
                                 or (var!="ST2" and var in clin[g]["genes"]))
            tag=" [MULTIPAIS]" if p in multi else ""
            print(f"    {pais:10s} {p:16s} n={len(l):3d}  positivos={pos:3d}{tag}")

print("\n=== D. PRJNA1404816: longitudes de contig (1 contig, 0 rep) ===")
cl=defaultdict(list)
for r in tsv(f"{B}/datos/contig_len.tsv"):
    if len(r)>=3: cl[norm(r[0])].append(int(r[2]))
for g in sorted(g for g in clin if bp.get(g)=="PRJNA1404816"):
    L=sorted(cl.get(g,[]),reverse=True)
    print(f"  {g:20s} contigs={len(L):3d}  mayores={L[:4]}  reps=[{clin[g]['rep']}]")

print("\n=== E. contigs de los 58 ensamblados propios ===")
fa=glob.glob(f"{B}/datos/ensamblados_63/*.fna")
print(f"  archivos .fna encontrados: {len(fa)}")
for f_ in sorted(fa)[:3]:
    n=sum(1 for l in open(f_) if l.startswith(">"))
    print(f"    {os.path.basename(f_):40s} {n} contigs")

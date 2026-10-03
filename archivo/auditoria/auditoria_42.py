#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cierre de la seccion 4.2: inventario, KL/OCL, rep, Tabla 2, eventos."""
import re, os
from collections import Counter, defaultdict
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]
rows = tsv(f"{B}/resultados/tabla_clinica_900.tsv")
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)", rows[0][0]): rows = rows[1:]
C=["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
clin={norm(r[0]):dict(zip(C,(r+[""]*10)[:10])) for r in rows}
bp={norm(r[0]):r[1] for r in tsv(f"{B}/datos/acc_bioproject_944.tsv")
    if len(r)>=2 and r[1].startswith("PRJ")}
N=len(clin)

print("=== A. inventario de carbapenemasas sobre los 900 (V4) ===")
inv=Counter()
for d in clin.values():
    for g in d["genes"].split(","):
        g=g.strip()
        if g and g not in ("-","NA"): inv[g]+=1
for g,n in inv.most_common(): print(f"  {g:16s} {n}")
sc=sum(1 for d in clin.values() if d["carb"].strip().lower() in ("no","-",""))
print(f"  sin carbapenemasa (col.5): {sc} ({100*sc/N:.1f} %)")

print("\n=== B. KL y OCL mas frecuentes ===")
for col in ("kl","ocl"):
    c=Counter(d[col] for d in clin.values())
    print(f"  {col.upper()} top5:", [(k,v,f"{100*v/N:.1f}%") for k,v in c.most_common(5)])

print("\n=== C. tipificacion de replicones ===")
con=sum(1 for d in clin.values() if d["rep"] not in ("","-","NA","ninguno"))
print(f"  con al menos un rep: {con} ({100*con/N:.1f} %)   (manuscrito: 88,6 %)")
rc=Counter(t.strip() for d in clin.values() for t in d["rep"].split(",")
           if t.strip() and t.strip() not in ("-","ninguno","NA"))
print("  tipos mas frecuentes:", rc.most_common(6))

print("\n=== D. eventos desduplicados (pais+ST+BioProject+anio) ===")
ev=defaultdict(list)
for g,d in clin.items(): ev[(d["pais"],d["st"],bp.get(g,"?"),d["anio"])].append(g)
print(f"  eventos totales: {len(ev)}   (manuscrito: 360)")
sinst=sum(1 for k in ev if not k[1].isdigit())
print(f"  eventos cuya clave tiene ST ausente: {sinst}  (punto 35)")
pe=Counter(k[0] for k in ev)
print("\n  === Tabla 1: genomas y eventos por pais ===")
print(f"  {'pais':12s} {'genomas':>8s} {'eventos':>8s} {'redund.':>8s} {'BioProj':>8s}")
for p,n in Counter(d["pais"] for d in clin.values()).most_common():
    npj=len({bp.get(g,'?') for g,d in clin.items() if d['pais']==p})
    print(f"  {p:12s} {n:8d} {pe[p]:8d} {100*(1-pe[p]/n):7.0f}% {npj:8d}")

print("\n  === Tabla 2: carbapenemasas por pais sobre eventos ===")
print(f"  {'pais':12s} {'ev':>4s} {'OXA-23':>8s} {'OXA-72':>8s} {'OXA-58':>8s} {'NDM':>8s}")
for p,_ in Counter(d["pais"] for d in clin.values()).most_common():
    ks=[k for k in ev if k[0]==p]; n=len(ks)
    def f(t):
        return 100*sum(1 for k in ks if any(t in clin[g]["genes"] for g in ev[k]))/n
    print(f"  {p:12s} {n:4d} {f('blaOXA-23'):7.1f}% {f('blaOXA-72'):7.1f}% "
          f"{f('blaOXA-58'):7.1f}% {f('blaNDM'):7.1f}%")

print("\n=== E. blaOXA-72 sin desduplicar ===")
for p in ("Peru","Perú","Ecuador"):
    gs=[g for g,d in clin.items() if d["pais"]==p]
    if not gs: continue
    k=sum(1 for g in gs if "blaOXA-72" in clin[g]["genes"])
    print(f"  {p}: {k}/{len(gs)} = {100*k/len(gs):.1f} %")

print("\n=== F. determinantes infrecuentes: linajes y anios ===")
raros=[g for g,d in clin.items() if any(x in d["genes"] for x in
       ("blaKPC","blaIMP","blaGES","blaOXA-241","blaOXA-231","blaOXA-143"))]
print(f"  genomas: {len(raros)}  linajes: {sorted({clin[g]['st'] for g in raros})}")
print(f"  anios: {sorted({clin[g]['anio'] for g in raros})}")
for g in raros: print(f"    {g:20s} {clin[g]['pais']:10s} ST{clin[g]['st']:6s} {clin[g]['anio']:5s} {clin[g]['genes']}")

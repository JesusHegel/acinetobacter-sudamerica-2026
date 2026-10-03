#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Puntos 49 (r3-T18) y 50 (regla de asignacion gen->evento)."""
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

print("=== PUNTO 49: los genomas con r3-T18 ===")
t18=[g for g,d in clin.items() if "r3-T18" in d["rep"].split(",")
     or "r3-T18" in [x.strip() for x in d["rep"].split(",")]]
print(f"  total: {len(t18)}   (4.6 solo explica 23 dentro de ST2 + 1 colombiano)")
print(f"  por pais: {dict(Counter(clin[g]['pais'] for g in t18))}")
print(f"  por ST:   {dict(Counter('ST'+clin[g]['st'] for g in t18))}")
print(f"  con blaOXA-72: {sum(1 for g in t18 if 'blaOXA-72' in clin[g]['genes'])}")
print("\n  --- los que NO son ST2/ST2724 ---")
for g in sorted(t18, key=lambda x:(clin[x]['pais'],clin[x]['st'])):
    if clin[g]['st'] not in ("2","2724"):
        d=clin[g]
        print(f"    {g:20s} {d['pais']:10s} ST{d['st']:6s} {d['anio']:5s} {d['kl']:6s} "
              f"{bp.get(g,'?'):14s} genes=[{d['genes']}]")

print("\n=== PUNTO 50: regla de asignacion gen -> evento ===")
ev=defaultdict(list)
for g,d in clin.items(): ev[(d["pais"],d["st"],bp.get(g,"?"),d["anio"])].append(g)
multi=[k for k,v in ev.items() if len(v)>1]
het=[k for k in multi if len({clin[g]["genes"] for g in ev[k]})>1]
print(f"  eventos: {len(ev)}   con >1 genoma: {len(multi)}   HETEROGENEOS: {len(het)}")
print(f"  (en {len(het)} eventos la regla elegida cambia el resultado)\n")
GEN=["blaOXA-23","blaOXA-72","blaOXA-58","blaNDM"]
print(f"  {'pais':11s} {'ev':>4s} " + " ".join(f"{g.replace('bla',''):>18s}" for g in GEN))
print(f"  {'':11s} {'':>4s} " + " ".join(f"{'cualq/repr/mayor':>18s}" for _ in GEN))
for p,_ in Counter(d["pais"] for d in clin.values()).most_common():
    ks=[k for k in ev if k[0]==p]; n=len(ks); out=[]
    for t in GEN:
        a=sum(1 for k in ks if any(t in clin[g]["genes"] for g in ev[k]))
        r=sum(1 for k in ks if t in clin[sorted(ev[k])[0]]["genes"])
        m=sum(1 for k in ks if sum(1 for g in ev[k] if t in clin[g]["genes"])*2>len(ev[k]))
        out.append(f"{100*a/n:5.1f}/{100*r/n:5.1f}/{100*m/n:5.1f}")
    print(f"  {p:11s} {n:4d} " + " ".join(f"{x:>18s}" for x in out))

print("\n=== eventos con clave de ST ausente (punto 35) ===")
sa=[k for k in ev if not k[1].isdigit()]
print(f"  {len(sa)} eventos, {sum(len(ev[k]) for k in sa)} genomas")
for k in sorted(sa)[:8]: print(f"    {k} -> {len(ev[k])} genomas")

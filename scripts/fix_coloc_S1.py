#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recalcula oxa72_contexto exigiendo rep en EL MISMO contig que el gen."""
import os, re, csv
from collections import defaultdict
B = os.path.expanduser("~/abaumannii")
def tsv(p):
    return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()

# contigs con blaOXA-72, del rastreo dirigido
hits = defaultdict(set)
for r in tsv(f"{B}/tmp_rescan/hits.tsv"):
    g, ctg = r[0].split("|", 1)
    hits[norm(g)].add(ctg)

# reps por contig
apt = defaultdict(lambda: defaultdict(set))
for r in tsv(f"{B}/resultados/apt_v2/apt_hits_v2.tsv")[1:]:
    apt[norm(r[0])][r[2]].add(r[1])

S1 = tsv(f"{B}/resultados/TablaS1_900_genomas.tsv")
hdr = S1[0]
rows = [dict(zip(hdr, r + [""]*len(hdr))) for r in S1[1:]]

port = coloc = noev = 0
cambios = []
for d in rows:
    a = d["accesion"]; prev = d["oxa72_contexto"]
    if a not in hits:
        d["oxa72_contexto"] = ""
        if prev: cambios.append((a, prev, "(vacio)"))
        continue
    port += 1
    reps = set()
    for c in hits[a]:                 # SOLO contigs que contienen el gen
        reps |= apt[a].get(c, set())
    nuevo = ",".join(sorted(reps)) if reps else "no_evaluable"
    if reps: coloc += 1
    else: noev += 1
    if nuevo != prev: cambios.append((a, prev, nuevo))
    d["oxa72_contexto"] = nuevo

print(f"=== recuento con criterio de mismo contig ===")
print(f"  portadores      : {port}")
print(f"  co-localizados  : {coloc}")
print(f"  no evaluables   : {noev}")
print(f"\n=== cambios respecto de la version anterior: {len(cambios)} ===")
for a, v, n in cambios[:20]:
    print(f"  {a:20s} [{v}] -> [{n}]")

with open(f"{B}/resultados/TablaS1_900_genomas.tsv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t")
    w.writeheader()
    for d in rows: w.writerow(d)
print(f"\n=== tabla reescrita ===")

# desglose por region
from collections import Counter
AND = {"Peru","Perú","Ecuador"}
cc2 = [d for d in rows if d["oxa72_contexto"] not in ("", "no_evaluable")
       and d["pais"] in AND and d["st_pasteur"] in ("2","2724")]
print(f"\n  CC2 andino co-localizado: {len(cc2)}  (manuscrito: 24)")
print(f"  reps del CC2 andino: {dict(Counter(d['oxa72_contexto'] for d in cc2))}")
br = [d for d in rows if d["oxa72_contexto"] not in ("", "no_evaluable") and d["pais"]=="Brasil"]
print(f"  Brasil co-localizado: {len(br)}  (manuscrito: 35)")
print(f"\n  por pais: {dict(Counter(d['pais'] for d in rows if d['oxa72_contexto'] not in ('','no_evaluable')))}")

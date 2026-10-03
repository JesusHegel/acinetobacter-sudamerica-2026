#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A1: tabla de metadatos para el deposito de los 58 ensamblados."""
import os, re, glob, csv
from collections import Counter
B = os.path.expanduser("~/abaumannii")
def tsv(p):
    return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]

M = tsv(f"{B}/metadatos/maestra_curada.tsv")
h = {n: i for i, n in enumerate(M[0])}
mae = {}
for r in M[1:]:
    if len(r) > h["Run"] and r[h["Run"]]:
        for run in r[h["Run"]].split(","):
            mae[run.strip()] = r

C = ["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
rows = tsv(f"{B}/resultados/tabla_clinica_900.tsv")
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)", rows[0][0]): rows = rows[1:]
clin = {r[0]: dict(zip(C, (r+[""]*10)[:10])) for r in rows}
props = sorted(g for g in clin if not g.startswith("GC"))

def g(r, col):
    i = h.get(col)
    if i is None or i >= len(r): return ""
    v = r[i].strip().strip('"')
    return "" if v.upper() in ("NA","NAN","-","") else v

out = []
falt = Counter()
for acc in props:
    r = mae.get(acc)
    if not r:
        falt["sin_fila_en_maestra"] += 1
        out.append({"run": acc}); continue
    d = {
        "run": acc,
        "biosample": g(r,"BioSample"),
        "bioproject_origen": g(r,"BioProject"),
        "organism": "Acinetobacter baumannii",
        "isolate": g(r,"Isolate") or g(r,"Strain"),
        "collection_date": g(r,"Collection date") or clin[acc]["anio"],
        "geo_loc_name": g(r,"Location") or clin[acc]["pais"],
        "isolation_source": g(r,"Isolation source"),
        "host": g(r,"Host"),
        "collected_by": g(r,"Collected by"),
        "platform": g(r,"Platform"),
        "st_pasteur": clin[acc]["st"],
        "kl": clin[acc]["kl"], "ocl": clin[acc]["ocl"],
        "carbapenemasas": clin[acc]["genes"],
    }
    for k in ("biosample","isolation_source","host","collected_by","collection_date"):
        if not d[k]: falt[k] += 1
    out.append(d)

cols = ["run","biosample","bioproject_origen","organism","isolate","collection_date",
        "geo_loc_name","isolation_source","host","collected_by","platform",
        "st_pasteur","kl","ocl","carbapenemasas"]
p = f"{B}/deposito/metadatos_58_ensamblados.tsv"
with open(p,"w",newline="",encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for d in out: w.writerow(d)

print(f"=== tabla escrita: {p}  ({len(out)} filas) ===\n")
print("=== campos faltantes ===")
print(f"  {dict(falt) if falt else 'ninguno: la tabla esta completa'}\n")
print("=== BioProjects de origen ===")
for k,v in Counter(d.get("bioproject_origen","?") for d in out).most_common():
    print(f"  {k:16s} {v:3d} genomas")
print("\n=== fuentes de aislamiento ===")
for k,v in Counter(d.get("isolation_source","(vacio)") for d in out).most_common(8):
    print(f"  {k:24s} {v}")
print("\n=== colectores ===")
for k,v in Counter(d.get("collected_by","(vacio)") for d in out).most_common():
    print(f"  {k:28s} {v}")
print("\n=== plataformas ===")
print(f"  {dict(Counter(d.get('platform','(vacio)') for d in out))}")
print("\n=== primeras 3 filas ===")
for d in out[:3]:
    print("  " + " | ".join(f"{k}={d.get(k,'')}" for k in cols[:10]))

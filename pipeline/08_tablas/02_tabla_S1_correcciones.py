#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Correcciones a la Tabla S1: oxa72 con rastreo dirigido, contigs propios, ciudad, ST Oxford."""
import os, re, csv, glob
from collections import Counter, defaultdict
B = os.path.expanduser("~/abaumannii")
def tsv(p):
    return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()

S1 = tsv(f"{B}/resultados/TablaS1_900_genomas.tsv")
hdr = S1[0]; rows = [dict(zip(hdr, r + [""]*len(hdr))) for r in S1[1:]]
print(f"filas cargadas: {len(rows)}")

# 1. oxa72 desde el rastreo dirigido
hits = defaultdict(list)
for r in tsv(f"{B}/tmp_rescan/hits.tsv"):
    g_, ctg = r[0].split("|", 1)
    hits[norm(g_)].append(ctg)
apt = defaultdict(lambda: defaultdict(set))
for r in tsv(f"{B}/resultados/apt_v2/apt_hits_v2.tsv")[1:]:
    apt[norm(r[0])][r[2]].add(r[1])
port = coloc = 0
for d in rows:
    a = d["accesion"]
    if a not in hits:
        d["oxa72_contexto"] = ""; continue
    port += 1
    reps = sorted({t for c in set(hits[a]) for t in apt[a].get(c, set())})
    if reps:
        d["oxa72_contexto"] = ",".join(reps); coloc += 1
    else:
        d["oxa72_contexto"] = "no_evaluable"
print(f"\n1. blaOXA-72: portadores={port} (esperado 73)  co-localizados={coloc} (esperado 63)")

# 2. contigs y longitud de los 58 propios
n = 0
for d in rows:
    if d["tipo_registro"] != "ensamblado_propio": continue
    f = f"{B}/datos/ensamblados_63/{d['accesion']}.fna"
    if not os.path.exists(f): continue
    L = []; cur = 0
    for l in open(f):
        if l.startswith(">"):
            if cur: L.append(cur)
            cur = 0
        else: cur += len(l.strip())
    if cur: L.append(cur)
    d["n_contigs"] = len(L); d["longitud_total"] = sum(L); n += 1
print(f"2. contigs calculados para {n} ensamblados propios")

# 3. ciudad desde la columna curada
M = tsv(f"{B}/metadatos/maestra_curada.tsv")
h = {x: i for i, x in enumerate(M[0])}
ciu = {}
for r in M[1:]:
    for col in ("Assembly", "Run"):
        i = h.get(col)
        if i is not None and i < len(r) and r[i]:
            for v in r[i].split(","):
                j = h.get("ciudad")
                val = r[j].strip() if j is not None and j < len(r) else ""
                ciu[v.strip()] = "" if val.upper() in ("NA","") else val
nc = 0
for d in rows:
    d["ciudad"] = ciu.get(d["accesion"], "")
    if d["ciudad"]: nc += 1
print(f"3. ciudad real disponible en {nc}/{len(rows)} ({100*nc/len(rows):.1f} %)")

# 4. ST Oxford
so = sum(1 for d in rows if d["st_oxford"] not in ("", "-", "NA"))
print(f"4. ST Oxford asignado: {so}/{len(rows)} ({100*so/len(rows):.1f} %)")
print(f"   valores vacios o '-': {dict(Counter(d['st_oxford'] for d in rows if d['st_oxford'] in ('','-','NA')))}")

p = f"{B}/resultados/TablaS1_900_genomas.tsv"
with open(p, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t")
    w.writeheader()
    for d in rows: w.writerow(d)
print(f"\n=== tabla reescrita ===")
print("\ncontroles finales:")
print(f"  portadores blaOXA-72     : {port}   (manuscrito: 73)")
print(f"  co-localizados           : {coloc}  (manuscrito: 63)")
print(f"  con n_contigs            : {sum(1 for d in rows if str(d['n_contigs']).strip() not in ('','NA'))}/900")
print(f"  con ciudad               : {nc}/900")

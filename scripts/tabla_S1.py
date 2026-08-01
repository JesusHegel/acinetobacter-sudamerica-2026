#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tabla Suplementaria S1: los 900 genomas clinicos con metadatos y tipificacion."""
import os, re, csv
from collections import Counter
B = os.path.expanduser("~/abaumannii")
def tsv(p):
    return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]

C = ["accesion","pais","st_pasteur","anio","carbapenemasa","genes_carb",
     "st_oxford","kl","ocl","rep_apt"]
rows = tsv(f"{B}/resultados/tabla_clinica_900.tsv")
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)", rows[0][0]): rows = rows[1:]
clin = {r[0]: dict(zip(C,(r+[""]*10)[:10])) for r in rows}
print(f"genomas clinicos: {len(clin)}")

bp = {r[0]: r[1] for r in tsv(f"{B}/datos/acc_bioproject_944.tsv")
      if len(r) >= 2 and r[1].startswith("PRJ")}
print(f"con BioProject: {sum(1 for g in clin if g in bp)}")

M = tsv(f"{B}/metadatos/maestra_curada.tsv")
h = {n: i for i, n in enumerate(M[0])}
mae = {}
for r in M[1:]:
    for col in ("Assembly", "Run"):
        i = h.get(col)
        if i is not None and i < len(r) and r[i]:
            for v in r[i].split(","):
                v = v.strip()
                if v: mae[v] = r
def g(r, col):
    i = h.get(col)
    if r is None or i is None or i >= len(r): return ""
    v = r[i].strip().strip('"')
    return "" if v.upper() in ("NA","NAN","-","") else v

ncont = Counter()
clen = {}
for r in tsv(f"{B}/datos/contig_len.tsv"):
    if len(r) >= 3:
        ncont[r[0]] += 1
        clen[r[0]] = clen.get(r[0], 0) + int(r[2])

col = {n: i for i, n in enumerate(tsv(f"{B}/resultados/apt_v2/colocalizacion.tsv")[0])}
cr = tsv(f"{B}/resultados/apt_v2/colocalizacion.tsv")[1:]
o72 = {}
for r in cr:
    if r[col["gen"]] == "blaOXA-72":
        acc = r[col["genoma"]]
        est = r[col["estado"]]
        reps = r[col["reps_mismo_contig"]]
        if est == "COLOCALIZADO" and reps not in ("-","","NA"):
            o72[acc] = reps
        else:
            o72.setdefault(acc, "no_evaluable")

OUT = []
falt = Counter()
for acc in sorted(clin):
    d = clin[acc]
    m = mae.get(acc)
    if m is None: falt["sin_metadatos"] += 1
    OUT.append({
        "accesion": acc,
        "tipo_registro": "ensamblado_publico" if acc.startswith("GC") else "ensamblado_propio",
        "bioproject": bp.get(acc, ""),
        "biosample": g(m, "BioSample"),
        "pais": d["pais"],
        "ciudad": g(m, "Location"),
        "anio": d["anio"],
        "fecha_colecta": g(m, "Collection date"),
        "fuente_aislamiento": g(m, "Isolation source"),
        "tipo_aislamiento": g(m, "Isolation type"),
        "hospedero": g(m, "Host"),
        "colector": g(m, "Collected by"),
        "plataforma": g(m, "Platform"),
        "n_contigs": ncont.get(acc, ""),
        "longitud_total": clen.get(acc, ""),
        "st_pasteur": d["st_pasteur"],
        "st_oxford": d["st_oxford"],
        "kl": d["kl"],
        "ocl": d["ocl"],
        "carbapenemasa_adquirida": d["carbapenemasa"],
        "genes_carbapenemasa": d["genes_carb"],
        "tipos_rep_apt": d["rep_apt"],
        "oxa72_contexto": o72.get(acc, ""),
    })

cols = list(OUT[0].keys())
p = f"{B}/resultados/TablaS1_900_genomas.tsv"
with open(p, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for d in OUT: w.writerow(d)
print(f"\n=== escrita: {p} ===")
print(f"  filas: {len(OUT)}   columnas: {len(cols)}")

print("\n=== completitud por columna ===")
for c in cols:
    n = sum(1 for d in OUT if str(d[c]).strip() not in ("", "NA", "-"))
    marca = "" if n == len(OUT) else f"   <-- faltan {len(OUT)-n}"
    print(f"  {c:26s} {n:4d}/{len(OUT)} ({100*n/len(OUT):5.1f} %){marca}")

print("\n=== controles de coherencia ===")
print(f"  ensamblados propios: {sum(1 for d in OUT if d['tipo_registro']=='ensamblado_propio')}  (esperado 58)")
print(f"  BioProjects unicos:  {len({d['bioproject'] for d in OUT if d['bioproject']})}  (esperado 104)")
print(f"  portadores blaOXA-72 con contexto: {sum(1 for d in OUT if d['oxa72_contexto'])}")
print(f"    de ellos co-localizados: {sum(1 for d in OUT if d['oxa72_contexto'] and d['oxa72_contexto']!='no_evaluable')}  (esperado 63)")
print(f"  sin ST Pasteur: {sum(1 for d in OUT if d['st_pasteur'] in ('-',''))}  (esperado 25)")
print(f"  paises: {len({d['pais'] for d in OUT})}  (esperado 9)")

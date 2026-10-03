#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B2: replicon del segundo ST108 ecuatoriano. B3: desglose de los 25 sin ST."""
import os, re
from collections import Counter, defaultdict
B = os.path.expanduser("~/abaumannii")
def tsv(p):
    return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]

S1 = tsv(f"{B}/resultados/TablaS1_900_genomas.tsv")
h = {n: i for i, n in enumerate(S1[0])}
rows = [dict(zip(S1[0], r + [""]*len(S1[0]))) for r in S1[1:]]

print("=== B2. genomas ecuatorianos ===")
for d in rows:
    if d["pais"] == "Ecuador":
        print(f"  {d['accesion']:20s} ST{d['st_pasteur']:6s} {d['kl']:7s} {d['anio']:5s} "
              f"{d['bioproject']:14s} reps=[{d['tipos_rep_apt']}]")
        print(f"      genes=[{d['genes_carbapenemasa']}]  oxa72_ctx=[{d['oxa72_contexto']}]")

print("\n=== B2b. los ST108: detalle de contigs ===")
apt = defaultdict(list)
for r in tsv(f"{B}/resultados/apt_v2/apt_hits_v2.tsv")[1:]:
    apt[r[0]].append(r)
cr = tsv(f"{B}/resultados/apt_v2/colocalizacion.tsv")
ci = {n: i for i, n in enumerate(cr[0])}
for d in rows:
    if d["pais"] == "Ecuador" and d["st_pasteur"] == "108":
        a = d["accesion"]
        print(f"\n  --- {a} ---")
        for r in apt.get(a, []):
            print(f"    rep {r[1]:10s} contig={r[2]:24s} len={r[5]:>8s} id={r[6]} cov={r[7]}")
        for r in cr[1:]:
            if r[ci["genoma"]] == a and r[ci["gen"]] == "blaOXA-72":
                print(f"    blaOXA-72: contig={r[ci['gen_contig']]} estado={r[ci['estado']]}")
                print(f"      mismo_contig=[{r[ci['reps_mismo_contig']]}]  otros=[{r[ci['reps_otros_contigs']]}]")

print("\n=== B3. los 25 genomas sin ST Pasteur ===")
sin = [d for d in rows if d["st_pasteur"] in ("-", "")]
print(f"  total: {len(sin)}")
print(f"  por pais: {dict(Counter(d['pais'] for d in sin))}")
print(f"  por tipo: {dict(Counter(d['tipo_registro'] for d in sin))}")
print(f"  por BioProject: {dict(Counter(d['bioproject'] for d in sin))}")

ML = None
for c in ("resultados/tipificacion/mlst_pasteur.tsv", "resultados/mlst_pasteur.tsv"):
    if os.path.exists(f"{B}/{c}"): ML = f"{B}/{c}"; break
print(f"\n  fichero MLST: {ML if ML else 'NO ENCONTRADO'}")
if ML:
    ACC = re.compile(r"(GC[AF]_\d+\.\d+|[EDS]RR\d+)")
    dat = {}
    for l in open(ML):
        f = l.rstrip("\n").split("\t")
        if len(f) < 3: continue
        m = ACC.search(f[0])
        if m: dat[m.group(1)] = f
    nuevo = truncado = ausente = 0
    print(f"\n  {'genoma':20s} {'pais':10s} perfil alelico")
    for d in sin:
        f = dat.get(d["accesion"])
        if not f: ausente += 1; continue
        alelos = " ".join(f[3:10])
        tiene_nuevo = "~" in alelos or "?" in alelos
        tiene_falta = "-" in alelos
        cat = "ALELO_NUEVO" if tiene_nuevo else ("LOCUS_AUSENTE" if tiene_falta else "COMBINACION_NUEVA")
        if cat == "ALELO_NUEVO": nuevo += 1
        elif cat == "LOCUS_AUSENTE": truncado += 1
        print(f"  {d['accesion']:20s} {d['pais']:10s} {cat:18s} {alelos}")
    comb = len(sin) - nuevo - truncado - ausente
    print(f"\n  === desglose ===")
    print(f"    combinacion alelica no descrita en PubMLST : {comb}")
    print(f"    alelo nuevo en al menos un locus           : {nuevo}")
    print(f"    locus no detectado (fragmentacion)         : {truncado}")
    print(f"    sin fila en el fichero MLST                : {ausente}")

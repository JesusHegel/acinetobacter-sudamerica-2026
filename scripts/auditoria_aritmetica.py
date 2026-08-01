#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Puntos 16-20 y 47: auditoria de cifras del manuscrito."""
import re, os, glob
from collections import Counter
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]

rows = tsv(f"{B}/resultados/tabla_clinica_900.tsv")
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)", rows[0][0]): rows = rows[1:]
C = ["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
clin = {norm(r[0]): dict(zip(C,(r+[""]*10)[:10])) for r in rows}
N = len(clin); print(f"conjunto clinico: {N}")

print("\n=== PUNTO 17: asignacion de ST Pasteur ===")
sin = [g for g,d in clin.items() if not d["st"].isdigit()]
con = N - len(sin)
print(f"  con ST: {con} ({100*con/N:.1f} %)   sin ST: {len(sin)}   (manuscrito: 97,2 % y 30)")
print(f"  valores literales:", dict(Counter(clin[g]['st'] for g in sin)))
print(f"  por pais:", dict(Counter(clin[g]['pais'] for g in sin)))
print(f"  propios (no GCA):", sum(1 for g in sin if not g.startswith("GC")))

print("\n=== PUNTO 18: cinco linajes mas frecuentes ===")
cnt = Counter(d["st"] for d in clin.values() if d["st"].isdigit())
top = cnt.most_common(5); s = sum(n for _,n in top)
for st,n in top: print(f"  ST{st:6s} n={n}")
print(f"  suma={s}   sobre 900 = {100*s/N:.1f} %   sobre asignados = {100*s/con:.1f} %   (manuscrito: 77 %)")

print("\n=== PUNTO 16: ST2 y ST79, tabla clinica vs matriz cgMLST ===")
for st,d in (("2","cgmlst_st2_eval"),("79","cgmlst_st79_eval")):
    cs = {g for g,v in clin.items() if v["st"]==st}
    mp = f"{B}/resultados/{d}/distance_matrix_symmetric.tsv"
    if not os.path.exists(mp): print(f"  ST{st}: falta {mp}"); continue
    with open(mp) as f: names=f.readline().rstrip("\n").split("\t")[1:]
    mat = {norm(n) for n in names}
    print(f"  ST{st}: clinicos={len(cs)}  matriz={len(mat)}")
    print(f"    en matriz y NO clinicos: {len(mat-cs)}  {sorted(mat-cs)[:6]}")
    print(f"    clinicos y NO en matriz: {len(cs-mat)}  {sorted(cs-mat)[:6]}")
    if cs-mat:
        print("    paises de los ausentes:", dict(Counter(clin[g]['pais'] for g in cs-mat)))
        print("    propios (no GCA):", sum(1 for g in cs-mat if not g.startswith("GC")))

print("\n=== PUNTO 20: esquema APT ===")
apt = f"{B}/AcinetobacterPlasmidTyping-main/2026_May_rep_DNA-seqs_V3.fasta"
h = [l[1:].strip() for l in open(apt) if l.startswith(">")]
def ty(x):
    m = re.match(r"([A-Za-z0-9]+)[-_](T\d+)", x)
    return f"{m.group(1)}-{m.group(2)}" if m else x.split("_")[0]
t = {ty(x) for x in h}
print(f"  secuencias de referencia: {len(h)}   tipos distintos: {len(t)}   (manuscrito: 257 / doc: 262)")
print(f"  prefijos:", dict(Counter(x.split('-')[0] for x in t)))

print("\n=== KL y OCL distintos (verificacion de V37) ===")
for col,exp in (("kl",65),("ocl",12)):
    v = {d[col] for d in clin.values() if d[col] and d[col] not in ("-","NA","")}
    print(f"  {col.upper()}: {len(v)} distintos (V37: {exp})   sin asignar: "
          f"{sum(1 for d in clin.values() if d[col] in ('','-','NA'))}")

print("\n=== PUNTO 47: los dos ST2742 KL14 sospechosos ===")
for a in ("GCA_058672075.1","GCA_058672145.1","GCA_058667195.1"):
    f = glob.glob(f"{B}/datos/genomas_945/ncbi_dataset/data/{a}/*.fna")
    if not f: print(f"  {a}: sin FASTA"); continue
    L=[];cur=0
    for l in open(f[0]):
        if l.startswith(">"):
            if cur: L.append(cur)
            cur=0
        else: cur+=len(l.strip())
    if cur: L.append(cur)
    L.sort(reverse=True)
    print(f"  {a}: contigs={len(L)} total={sum(L)} mayor={L[0]} N50_aprox={L[len(L)//4]}")

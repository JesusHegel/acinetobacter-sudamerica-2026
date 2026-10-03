#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A1: inventario de los ensamblados propios para deposito en NCBI."""
import os, re, glob
from collections import Counter
B = os.path.expanduser("~/abaumannii")
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]

rows = tsv(f"{B}/resultados/tabla_clinica_900.tsv")
hdr = rows[0] if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)", rows[0][0]) else None
if hdr: rows = rows[1:]
C = ["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
clin = {r[0]: dict(zip(C,(r+[""]*10)[:10])) for r in rows}
nogca = sorted(g for g in clin if not g.startswith("GC"))
print(f"=== A. los {len(nogca)} clinicos sin accesion GCA ===")
print("  prefijos:", dict(Counter(re.match(r'^([A-Za-z]+)', g).group(1) for g in nogca)))
print("  por pais:", dict(Counter(clin[g]['pais'] for g in nogca)))
raros = [g for g in nogca if not re.match(r'^(ERR|SRR|DRR)\d+$', g)]
print(f"  con formato inesperado: {raros if raros else 'ninguno'}")

fa = sorted(glob.glob(f"{B}/datos/ensamblados_63/*.fna"))
ids_fa = {os.path.basename(f)[:-4] for f in fa}
print(f"\n=== B. archivos .fna: {len(fa)} ===")
print(f"  en la tabla clinica y con .fna : {len(set(nogca) & ids_fa)}")
print(f"  con .fna pero NO clinicos      : {sorted(ids_fa - set(nogca))}")
print(f"  clinicos SIN .fna              : {sorted(set(nogca) - ids_fa)}")

print("\n=== C. control de calidad para NCBI ===")
print(f"  {'genoma':16s} {'contigs':>8s} {'total_pb':>10s} {'min':>7s} {'<200':>5s} {'N':>7s} {'%N':>6s}")
prob = []
for f in fa:
    acc = os.path.basename(f)[:-4]
    if acc not in clin: continue
    L=[]; cur=0; ns=0
    for l in open(f):
        if l.startswith(">"):
            if cur: L.append(cur)
            cur=0
        else:
            s=l.strip(); cur+=len(s); ns+=s.upper().count("N")
    if cur: L.append(cur)
    tot=sum(L); short=sum(1 for x in L if x<200); pn=100*ns/tot if tot else 0
    if short or pn>1 or tot<3.5e6 or tot>4.5e6: prob.append(acc)
    if len(prob)<=8 or short or pn>1:
        print(f"  {acc:16s} {len(L):8d} {tot:10d} {min(L):7d} {short:5d} {ns:7d} {pn:6.2f}")
print(f"\n  genomas con alguna incidencia: {len(prob)}  {prob[:10] if prob else ''}")

print("\n=== D. campos de metadatos para el envio ===")
meta = ["organism","strain/isolate","collection_date","geo_loc_name",
        "isolation_source","host","host_disease","collected_by"]
print("  BioSample (paquete Pathogen.cl.1.0) requiere:")
for m in meta: print(f"    - {m}")
print("\n  disponibles en tabla_clinica_900.tsv:")
print(f"    - pais y anio: SI")
print(f"    - fuente de aislamiento, hospedero, colector: NO (columnas ausentes)")
print(f"\n  anios de los ensamblados propios: "
      f"{sorted({clin[g]['anio'] for g in nogca if g in ids_fa})}")

print("\n=== E. ficheros de metadatos crudos disponibles ===")
for d in ("metadatos","datos"):
    p=f"{B}/{d}"
    if os.path.isdir(p):
        for f in sorted(os.listdir(p))[:12]:
            fp=os.path.join(p,f)
            if os.path.isfile(fp):
                print(f"  {d}/{f}  ({os.path.getsize(fp)//1024} KB)")

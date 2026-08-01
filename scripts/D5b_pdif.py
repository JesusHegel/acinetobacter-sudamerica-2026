#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sitios pdif: busca las dos orientaciones documentadas, con tolerancia."""
import os, glob
B = os.path.expanduser("~/abaumannii")
def fasta(p):
    d={};k=None;buf=[]
    for l in open(p):
        if l.startswith(">"):
            if k: d[k]="".join(buf)
            k=l[1:].split()[0]; buf=[]
        else: buf.append(l.strip())
    if k: d[k]="".join(buf)
    return d
def fpath(a):
    g=glob.glob(f"{B}/datos/genomas_945/ncbi_dataset/data/{a}/*.fna")
    if g: return g[0]
    p=f"{B}/datos/ensamblados_63/{a}.fna"
    return p if os.path.exists(p) else None
def mm(a,b): return sum(1 for x,y in zip(a,b) if x!=y)

# Orientaciones observadas en V28/V30
PAT=[("C>D","ATTTCGTATAA","TTATGTTAAAT"),
     ("D>C","ATTTAACATAA","TTATACGAAAT")]

def buscar(seq, tol):
    n=len(seq); ext=seq+seq[:40]; out=[]
    for i in range(n):
        v=ext[i:i+28]
        if len(v)<28: break
        for ori,izq,der in PAT:
            d=mm(v[:11],izq)+mm(v[17:],der)
            if d<=tol: out.append((i,ori,v[11:17],d,v))
    return out

ELEM=[("ANDINO r3-T18","GCA_051942065.1","JBFECW010000069.1"),
      ("ANDINO r3-T18 (Ecuador)","GCA_029955785.1","JASBHG010000073.1"),
      ("CHILE r3-T18 vacio","GCA_024139015.1","CP076810.1"),
      ("BRASIL r3-T1","GCA_009829545.1",None),
      ("BRASIL r3-T1+T14","GCA_015537245.1",None)]
cr=[l.rstrip("\n").split("\t") for l in open(f"{B}/resultados/apt_v2/colocalizacion.tsv") if l.strip()]
ci={n:i for i,n in enumerate(cr[0])}
def ctg_oxa(acc):
    for r in cr[1:]:
        if r[ci["genoma"]]==acc and r[ci["gen"]]=="blaOXA-72" and r[ci["estado"]]=="COLOCALIZADO":
            return r[ci["gen_contig"]]

for tol in (0,1,2,3):
    print(f"\n{'='*64}\nTOLERANCIA: {tol} desajuste(s) en total sobre los 22 pb de brazos\n{'='*64}")
    for tag,acc,ctg in ELEM:
        c = ctg or ctg_oxa(acc)
        p=fpath(acc); F=fasta(p) if p else {}
        if not c or c not in F:
            print(f"\n{tag}: no hallado"); continue
        seq=F[c]; hits=buscar(seq,tol)
        print(f"\n{tag}  ({len(seq)} pb): {len(hits)} sitio(s)")
        for pos,ori,esp,d,v in sorted(hits):
            print(f"  pos {pos:6d}  {ori}  espaciador {esp}  desaj {d}   {v[:11]} {v[11:17]} {v[17:]}")

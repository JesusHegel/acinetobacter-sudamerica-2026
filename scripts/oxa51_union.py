#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Union de la delecion: repeticiones directas, y barrido de hits a nivel de familia."""
import os, glob
from collections import Counter
B = os.path.expanduser("~/abaumannii")
OUT = f"{B}/resultados/verif/amr_relax"
def fasta(p):
    d={};k=None;buf=[]
    for l in open(p):
        if l.startswith(">"):
            if k: d[k]="".join(buf)
            k=l[1:].split()[0]; buf=[]
        else: buf.append(l.strip())
    if k: d[k]="".join(buf)
    return d
def rc(s): return s[::-1].translate(str.maketrans("ACGTacgt","TGCAtgca"))
def lcp(a,b):
    i=0
    while i<min(len(a),len(b)) and a[i]==b[i]: i+=1
    return i
def lcs(a,b):
    i=0
    while i<min(len(a),len(b)) and a[-1-i]==b[-1-i]: i+=1
    return i

T=[("GCA_058667195.1","KL14"),("GCA_058672075.1","KL14"),("GCA_058672145.1","KL14"),
   ("GCA_058672675.1","ctrl"),("GCA_058673275.1","ctrl")]
nt={}; gen={}
print("=== A. secuencias nucleotidicas ===")
for acc,tag in T:
    t=f"{OUT}/{acc}.tsv"
    if not os.path.exists(t): print(f"  {acc}: falta {t}"); continue
    g=glob.glob(f"{B}/datos/genomas_945/ncbi_dataset/data/{acc}/*.fna")[0]
    F=fasta(g); gen[acc]=F
    for l in open(t).readlines()[1:]:
        r=l.rstrip("\n").split("\t")
        if "OXA-51" in r[6]:
            s=F[r[1]][int(r[2])-1:int(r[3])]
            if r[4]=="-": s=rc(s)
            nt[acc]=s
            print(f"  {acc} [{tag:4s}] {r[1]} {r[2]}-{r[3]} {r[4]} len={len(s)}")

k14=[a for a,t in T if t=="KL14" and a in nt]
ctl=[a for a,t in T if t=="ctrl" and a in nt]
print("\n=== B. union de la delecion ===")
if k14 and ctl:
    A,Cc=nt[k14[0]],nt[ctl[0]]
    print(f"  KL14 identicos entre si (nt): {all(nt[a]==A for a in k14)}")
    print(f"  controles identicos entre si (nt): {all(nt[a]==Cc for a in ctl)}")
    p,s=lcp(A,Cc),lcs(A,Cc)
    print(f"  len KL14={len(A)}  len ctrl={len(Cc)}  delecion={len(Cc)-len(A)} pb")
    print(f"  prefijo comun={p} pb   sufijo comun={s} pb   p+s={p+s}")
    rep=p+s-len(A)
    if rep>0:
        print(f"  -> REPETICION DIRECTA de {rep} pb en la union: {A[p-rep:p]}")
        print(f"     (firma de deslizamiento de hebra; delecion biologica plausible)")
    else:
        print("  -> sin repeticion directa detectable en la union")
    d=Cc[p:len(Cc)-s]
    print(f"  segmento ausente ({len(d)} pb): {d}")
    print(f"  contexto ctrl:  ...{Cc[max(0,p-18):p]} | {d[:18]}...{d[-18:]} | {Cc[len(Cc)-s:len(Cc)-s+18]}...")
    print(f"  contexto KL14:  ...{A[max(0,p-18):p]} | {A[p:p+18]}...")
    print("\n=== C. el segmento ausente, presente en otra parte del genoma? ===")
    for acc in k14:
        hit=[(c,seq.find(d)) for c,seq in gen[acc].items() if d in seq or rc(d) in seq]
        print(f"  {acc}: {hit if hit else 'no presente (descarta colapso por duplicacion)'}")

print("\n=== D. los 9 ST2742: BioProject y evento desduplicado ===")
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]
bp={r[0].split()[0]:r[1] for r in tsv(f"{B}/datos/acc_bioproject_944.tsv")
    if len(r)>=2 and r[1].startswith("PRJ")}
rows=[r for r in tsv(f"{B}/resultados/tabla_clinica_900.tsv") if len(r)>2 and r[2] in ("2742","188")]
ev=Counter()
for r in rows:
    p=bp.get(r[0],"?"); ev[(r[1],r[2],p,r[3])]+=1
    print(f"  {r[0]:20s} ST{r[2]:5s} {r[7]:6s} {p}")
print(f"\n  eventos distinguibles (pais+ST+BioProject+anio): {len(ev)}")
for k,v in ev.items(): print(f"    {k} -> {v} genomas")

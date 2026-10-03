#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B1: comparacion del elemento andino con plasmidos circulares cerrados de r3-T18."""
import os, re, glob, subprocess
from collections import defaultdict
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

SEQ=[("ANDINO_PE_8111","GCA_051942065.1","JBFECW010000069.1"),
     ("ANDINO_PE_7851","GCA_053320315.1","JBROGJ010000074.1"),
     ("ANDINO_EC_7795","GCA_029955785.1","JASBHG010000073.1"),
     ("CERRADO_CL_CP076810","GCA_024139015.1","CP076810.1")]
os.makedirs(f"{B}/tmp_cerrados", exist_ok=True)
fa=f"{B}/tmp_cerrados/comparacion.fa"
with open(fa,"w") as o:
    n=0
    for tag,acc,ctg in SEQ:
        p=fpath(acc)
        if not p: print(f"  [!] sin FASTA: {acc}"); continue
        F=fasta(p)
        if ctg not in F: print(f"  [!] contig ausente: {acc}/{ctg}"); continue
        o.write(f">{tag}|{len(F[ctg])}\n{F[ctg]}\n"); n+=1
print(f"=== A. secuencias extraidas: {n} ===")
S=fasta(fa)
for k in S: print(f"  {k}")

print("\n=== B. cobertura real (intervalos fusionados) ===")
r=subprocess.run(["blastn","-query",fa,"-subject",fa,"-evalue","1e-5",
                  "-outfmt","6 qseqid sseqid pident length qlen qstart qend"],
                 capture_output=True,text=True)
pair=defaultdict(list); ident=defaultdict(float); ql={}
for l in r.stdout.strip().split("\n"):
    if not l: continue
    f=l.split("\t")
    if f[0]==f[1]: continue
    ql[f[0]]=int(f[4])
    pair[(f[0],f[1])].append((min(int(f[5]),int(f[6])),max(int(f[5]),int(f[6]))))
    ident[(f[0],f[1])]=max(ident[(f[0],f[1])],float(f[2]))
def merge(iv):
    iv=sorted(iv); out=[]
    for x,y in iv:
        if out and x<=out[-1][1]+1: out[-1][1]=max(out[-1][1],y)
        else: out.append([x,y])
    return sum(y-x+1 for x,y in out)
print(f"  {'query':22s} {'subject':22s} {'cob%':>6s} {'id%':>7s} {'no_alin':>8s}")
for k in sorted(pair):
    c=merge(pair[k]); q=ql[k[0]]
    print(f"  {k[0].split('|')[0]:22s} {k[1].split('|')[0]:22s} {100*c/q:6.1f} {ident[k]:7.2f} {q-c:8d}")

print("\n=== C. el contig cerrado, ¿es circular en el deposito? ===")
for tag,acc,ctg in SEQ:
    if not tag.startswith("CERRADO"): continue
    p=fpath(acc)
    if p:
        for h in fasta(p):
            if h==ctg: print(f"  identificador: {h}")

print("\n=== D. blaOXA-72 en cada secuencia ===")
AMR=glob.glob(f"{os.path.expanduser('~')}/miniforge3/envs/*/share/amrfinderplus/data/*/AMR_CDS.fa")
if AMR:
    Ac=fasta(AMR[0])
    sel={k:v for k,v in Ac.items() if re.search(r"blaOXA-72\b",k)}
    q=f"{B}/tmp_cerrados/oxa72.fa"
    with open(q,"w") as o:
        for k,v in list(sel.items())[:1]: o.write(f">blaOXA-72\n{v}\n")
    r2=subprocess.run(["blastn","-query",q,"-subject",fa,"-evalue","1e-5",
                       "-outfmt","6 sseqid pident length qlen sstart send"],
                      capture_output=True,text=True)
    got=set()
    for l in r2.stdout.strip().split("\n"):
        if not l: continue
        f=l.split("\t"); got.add(f[0])
        print(f"  {f[0].split('|')[0]:22s} id={f[1]:>7s} aln={f[2]}/{f[3]} pos={f[4]}-{f[5]}")
    for k in S:
        if k not in got: print(f"  {k.split('|')[0]:22s} SIN blaOXA-72")
else:
    print("  [!] catalogo AMR_CDS no encontrado")

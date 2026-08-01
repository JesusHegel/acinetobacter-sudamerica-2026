#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Resuelve: dos contigs de 7851 pb, uno con blaOXA-72 y otro sin el."""
import re, os, glob, subprocess
from collections import defaultdict
B = os.path.expanduser("~/abaumannii")
FA = f"{B}/tmp_plas/armazon_r3t18.fa"
AMR = glob.glob(f"{os.path.expanduser('~')}/miniforge3/envs/abaumannii/share/amrfinderplus/data/*/AMR_CDS.fa")

def fasta(p):
    d={};k=None;buf=[]
    for l in open(p):
        if l.startswith(">"):
            if k: d[k]="".join(buf)
            k=l[1:].rstrip(); buf=[]
        else: buf.append(l.strip())
    if k: d[k]="".join(buf)
    return d
S = fasta(FA)
def tag(k): return k.split("|")[0]

print("=== A. comparacion directa de las dos secuencias de 7851 pb ===")
a=[v for k,v in S.items() if tag(k)=="ANDES_con_gen"]
b=[v for k,v in S.items() if tag(k)=="PERU_vacio_7851"]
if a and b:
    A,Bs=a[0],b[0]
    print(f"  identicas: {A==Bs}   len {len(A)} vs {len(Bs)}")
    if A!=Bs:
        d=[i for i in range(min(len(A),len(Bs))) if A[i]!=Bs[i]]
        print(f"  posiciones distintas (mismo offset): {len(d)}  primeras: {d[:10]}")

print("\n=== B. cobertura REAL (intervalos fusionados, sin doble conteo) ===")
r=subprocess.run(["blastn","-query",FA,"-subject",FA,"-evalue","1e-5",
                  "-outfmt","6 qseqid sseqid pident length qlen qstart qend"],
                 capture_output=True,text=True)
pair=defaultdict(list); ql={}
for l in r.stdout.strip().split("\n"):
    if not l: continue
    f=l.split("\t")
    if f[0]==f[1]: continue
    ql[f[0]]=int(f[4]); pair[(f[0],f[1])].append((min(int(f[5]),int(f[6])),max(int(f[5]),int(f[6]))))
def merge(iv):
    iv=sorted(iv); out=[]
    for x,y in iv:
        if out and x<=out[-1][1]+1: out[-1][1]=max(out[-1][1],y)
        else: out.append([x,y])
    return sum(y-x+1 for x,y in out)
KEY=["ANDES_con_gen","PERU_vacio_7851","PERU_vacio_7928","PERU_vacio_8149",
     "ECUADOR_con_gen","CHILE_8356","CHILE_CP076810"]
print(f"  {'query':18s} {'subject':18s} {'cob_real%':>10s} {'pb_no_alineados':>16s}")
for k,v in sorted(pair.items()):
    if tag(k[0]) in KEY and tag(k[1]) in KEY:
        c=merge(v); q=ql[k[0]]
        if 100*c/q>=40:
            print(f"  {tag(k[0]):18s} {tag(k[1]):18s} {100*c/q:9.1f}% {q-c:16d}")

print("\n=== C. familia OXA-24/40 contra cada contig ===")
if AMR:
    Ac=fasta(AMR[0])
    sel={k:v for k,v in Ac.items() if re.search(r"OXA-72|OXA-24|OXA-40",k,re.I)}
    q=f"{B}/tmp_plas/oxa24_refs.fa"
    with open(q,"w") as o:
        for k,v in sel.items(): o.write(f">{k.split()[0]}\n{v}\n")
    print(f"  referencias usadas: {len(sel)}")
    r2=subprocess.run(["blastn","-query",q,"-subject",FA,"-evalue","1e-5",
                       "-outfmt","6 qseqid sseqid pident length qlen sstart send"],
                      capture_output=True,text=True)
    h=defaultdict(list)
    for l in r2.stdout.strip().split("\n"):
        if not l: continue
        f=l.split("\t"); h[f[1]].append((f[0],float(f[2]),int(f[3]),int(f[4]),int(f[5]),int(f[6])))
    for k in S:
        hh=h.get(k,[])
        if hh:
            bst=max(hh,key=lambda x:(x[1],x[2]))
            print(f"  {tag(k):18s} n={len(hh):3d}  id={bst[1]:6.2f} aln={bst[2]}/{bst[3]} pos={bst[4]}-{bst[5]}")
        else:
            print(f"  {tag(k):18s} SIN hit OXA-24/40")

print("\n=== D. region entre los dos sitios pdif ===")
P1=re.compile(r"ATTTAACATAA.{6}TTATACGAAAT"); P2=re.compile(r"ATTTCGTATAA.{6}TTATGTTAAAT")
for k,s in S.items():
    pos=sorted([(m.start(),"D>C") for m in P1.finditer(s)]+[(m.start(),"C>D") for m in P2.finditer(s)])
    if len(pos)>=2:
        i,j=pos[0][0],pos[1][0]
        print(f"  {tag(k):18s} sitios en {[p for p,_ in pos]}  region interna={j-i-28} pb")
    else:
        print(f"  {tag(k):18s} sitios={len(pos)} {pos if pos else ''}")

print("\n=== E. AMRFinderPlus relajado sobre los 3 genomas 'vacios' ===")
for acc in ("GCA_053320395.1","GCA_053373735.1","GCA_037545805.1"):
    g=glob.glob(f"{B}/datos/genomas_945/ncbi_dataset/data/{acc}/*.fna")
    if not g: print(f"  {acc}: sin FASTA"); continue
    o=f"/tmp/amr_{acc}.tsv"
    if not os.path.exists(o):
        subprocess.run(["amrfinder","-n",g[0],"--organism","Acinetobacter_baumannii","--plus",
                        "--ident_min","0.5","--coverage_min","0.3","-o",o],
                       stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    ox=[l.split("\t") for l in open(o).readlines()[1:] if "OXA" in l]
    print(f"  {acc}: {len(ox)} hits OXA")
    for r3 in ox: print(f"     {r3[5]:14s} {r3[1]:24s} cov={r3[15]:>7s} id={r3[16]:>7s}")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Puntos 4,5,7,45: caracterizacion del gen OXA-51-like de los ST2742 KL14."""
import os, glob, subprocess
B = os.path.expanduser("~/abaumannii")
OUT = f"{B}/resultados/verif/amr_relax"; os.makedirs(OUT, exist_ok=True)
T = [("GCA_058667195.1","KL14"),("GCA_058672075.1","KL14"),("GCA_058672145.1","KL14"),
     ("GCA_058672675.1","KL9-ctrl"),("GCA_058673275.1","KL22-ctrl")]
C={}
for i,b1 in enumerate("TCAG"):
    for j,b2 in enumerate("TCAG"):
        for k,b3 in enumerate("TCAG"):
            C[b1+b2+b3]="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[i*16+j*4+k]
def tr(s): return "".join(C.get(s[i:i+3].upper(),"X") for i in range(0,len(s)-2,3))
def rc(s): return s[::-1].translate(str.maketrans("ACGTacgt","TGCAtgca"))
def fasta(p):
    d={};k=None;buf=[]
    for l in open(p):
        if l.startswith(">"):
            if k: d[k]="".join(buf)
            k=l[1:].split()[0]; buf=[]
        else: buf.append(l.strip())
    if k: d[k]="".join(buf)
    return d
def lcp(a,b):
    i=0
    while i<min(len(a),len(b)) and a[i]==b[i]: i+=1
    return i
def lcs(a,b):
    i=0
    while i<min(len(a),len(b)) and a[-1-i]==b[-1-i]: i+=1
    return i

prot={}
for acc,tag in T:
    g=glob.glob(f"{B}/datos/genomas_945/ncbi_dataset/data/{acc}/*.fna")
    if not g: print(f"{acc}: FASTA no encontrado"); continue
    o=f"{OUT}/{acc}.tsv"
    if not os.path.exists(o) or os.path.getsize(o)==0:
        subprocess.run(["amrfinder","-n",g[0],"--organism","Acinetobacter_baumannii","--plus",
                        "--ident_min","0.5","--coverage_min","0.3","-o",o],
                       stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    print(f"\n=== {acc}  [{tag}] ===")
    rows=[l.rstrip("\n").split("\t") for l in open(o)][1:]
    for r in rows:
        if "OXA" not in r[6] and "OXA" not in r[5]: continue
        print(f"  symbol={r[5]:12s} scope={r[7]:6s} method={r[12]:8s} "
              f"tgt={r[13]:>4s} ref={r[14]:>4s} cov={r[15]:>7s} id={r[16]:>7s} aln={r[17]:>4s}")
        print(f"     closest: {r[18]} / {r[19]}")
        print(f"     HMM:     {r[20]} / {r[21]}")
        if "OXA-51" in r[6] or r[5]=="blaOXA-407":
            ctg,ini,fin,st=r[1],int(r[2]),int(r[3]),r[4]
            seq=fasta(g[0]).get(ctg,"")[ini-1:fin]
            if st=="-": seq=rc(seq)
            p=tr(seq); prot[acc]=p
            print(f"     nt={len(seq)} aa={len(p)} stop_interno={'SI' if '*' in p[:-1] else 'no'}")
            print(f"     {p}")

print("\n=== COMPARACION KL14 vs control ===")
kl14=[a for a,t in T if t=="KL14" and a in prot]
ctrl=[a for a,t in T if "ctrl" in t and a in prot]
if len(kl14)>1:
    print("  KL14 identicos entre si:", all(prot[kl14[0]]==prot[a] for a in kl14[1:]))
if kl14 and ctrl:
    A,Bp=prot[kl14[0]],prot[ctrl[0]]
    p,s=lcp(A,Bp),lcs(A,Bp)
    print(f"  KL14 len={len(A)}  control len={len(Bp)}  dif={len(Bp)-len(A)}")
    print(f"  prefijo comun={p}  sufijo comun={s}  p+s={p+s}")
    if p+s>=len(A):
        print(f"  -> DELECION LIMPIA de {len(Bp)-len(A)} aa tras el residuo {p}")
        print(f"     ausente: {Bp[p:len(Bp)-s]}")
    else:
        print("  -> no es una delecion simple; hay sustituciones")
    for m in ("STFK","KTG","SAV"):
        print(f"     motivo {m}: KL14 pos={A.find(m)}  control pos={Bp.find(m)}")

print("\n=== A. genes core en la columna de carbapenemasas (punto 45) ===")
for acc,tag in T[:1]:
    for l in open(f"{OUT}/{acc}.tsv").readlines()[1:]:
        r=l.rstrip("\n").split("\t")
        if len(r)>11 and "CARBAPENEM" in r[11]:
            print(f"  {r[5]:14s} scope={r[7]:6s} subclass={r[11]}")

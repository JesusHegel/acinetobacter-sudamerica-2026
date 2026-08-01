#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura 4: matriz de cobertura reciproca entre grupos de elementos."""
import os, glob, subprocess
from collections import defaultdict
B = os.path.expanduser("~/abaumannii")
def fasta(p, full=False):
    d={};k=None;buf=[]
    for l in open(p):
        if l.startswith(">"):
            if k: d[k]="".join(buf)
            k=(l[1:].rstrip() if full else l[1:].split()[0]); buf=[]
        else: buf.append(l.strip())
    if k: d[k]="".join(buf)
    return d
def fpath(a):
    g=glob.glob(f"{B}/datos/genomas_945/ncbi_dataset/data/{a}/*.fna")
    if g: return g[0]
    p=f"{B}/datos/ensamblados_63/{a}.fna"
    return p if os.path.exists(p) else None
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]

cr=tsv(f"{B}/resultados/apt_v2/colocalizacion.tsv"); ci={n:i for i,n in enumerate(cr[0])}
S1=tsv(f"{B}/resultados/TablaS1_900_genomas.tsv"); h=S1[0]
clin={r[0]:dict(zip(h,r+[""]*len(h))) for r in S1[1:]}
port=defaultdict(list)
for r in cr[1:]:
    if r[ci["gen"]]=="blaOXA-72" and r[ci["estado"]]=="COLOCALIZADO":
        a=r[ci["genoma"]]; reps=r[ci["reps_mismo_contig"]]
        try: ln=int(r[ci["contig_len"]])
        except ValueError: ln=0
        if a in clin: port[reps].append((ln, a, r[ci["gen_contig"]], clin[a]["pais"]))

print("=== representantes disponibles por combinacion de replicones ===")
for reps in sorted(port, key=lambda k:-len(port[k])):
    port[reps].sort(reverse=True)
    ln,a,c,p = port[reps][0]
    print(f"  {reps:16s} n={len(port[reps]):2d}  mayor: {a} ({p}, {ln} pb)")

SEQ=[]
for reps,tag in [("r3-T18","ANDINO r3-T18"),("r3-T1,r3-T14","BRASIL r3-T1+T14"),
                 ("r3-T1","BRASIL r3-T1"),("r3-T14","BRASIL r3-T14"),
                 ("r3-T13","r3-T13"),("r3-T13,r3-T50","ECUADOR r3-T50"),
                 ("r3-T8,r3-T97","COLOMBIA r3-T8+T97")]:
    if reps in port:
        ln,a,c,p = port[reps][0]
        SEQ.append((f"{tag}", a, c))
SEQ.append(("CHILE r3-T18 vacio","GCA_024139015.1","CP076810.1"))

os.makedirs(f"{B}/tmp_fig4", exist_ok=True)
fa=f"{B}/tmp_fig4/grupos.fa"
with open(fa,"w") as o:
    n=0
    for tag,acc,ctg in SEQ:
        p=fpath(acc); F=fasta(p) if p else {}
        if ctg in F:
            o.write(f">{tag.replace(' ','_')}#{len(F[ctg])}\n{F[ctg]}\n"); n+=1
        else: print(f"  [!] no hallado: {tag} ({acc}/{ctg})")
print(f"\nsecuencias extraidas: {n}")

S=fasta(fa,full=True)
tags=[k.split("#")[0].replace("_"," ") for k in S]
L={k.split("#")[0].replace("_"," "):int(k.split("#")[1]) for k in S}
r=subprocess.run(["blastn","-query",fa,"-subject",fa,"-evalue","1e-5",
                  "-outfmt","6 qseqid sseqid pident length qlen qstart qend"],
                 capture_output=True,text=True)
iv=defaultdict(list)
for l in r.stdout.strip().split("\n"):
    if not l: continue
    f=l.split("\t"); q=f[0].split("#")[0].replace("_"," "); s=f[1].split("#")[0].replace("_"," ")
    iv[(q,s)].append((min(int(f[5]),int(f[6])),max(int(f[5]),int(f[6]))))
def merge(v):
    v=sorted(v); out=[]
    for x,y in v:
        if out and x<=out[-1][1]+1: out[-1][1]=max(out[-1][1],y)
        else: out.append([x,y])
    return sum(y-x+1 for x,y in out)
M={(q,s): (100.0 if q==s else (100*merge(iv[(q,s)])/L[q] if (q,s) in iv else 0.0))
   for q in tags for s in tags}

print("\n=== matriz de cobertura (fila = consulta) ===")
print(" "*24 + " ".join(f"{t[:9]:>10s}" for t in tags))
for q in tags:
    print(f"{q:24s}" + " ".join(f"{M[(q,s)]:10.1f}" for s in tags))

N=len(tags); C=54; X0,Y0=210,165; W=X0+N*C+150; H=Y0+N*C+80
def col(v):
    return ("#08306b" if v>=95 else "#2171b5" if v>=80 else "#6baed6"
            if v>=50 else "#bdd7e7" if v>=25 else "#eff3ff" if v>0 else "#f7f7f7")
s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
   f'<rect width="{W}" height="{H}" fill="white"/>']
for j,t in enumerate(tags):
    x=X0+j*C+C/2
    s.append(f'<text x="{x}" y="{Y0-8}" font-family="Helvetica,Arial" font-size="9.5" '
             f'transform="rotate(-55 {x} {Y0-8})">{t}</text>')
for i,q in enumerate(tags):
    y=Y0+i*C
    s.append(f'<text x="{X0-8}" y="{y+C/2+4}" text-anchor="end" font-family="Helvetica,Arial" font-size="9.5">{q}</text>')
    s.append(f'<text x="{X0+N*C+8}" y="{y+C/2+4}" font-family="Helvetica,Arial" font-size="8.5" fill="#666">{L[q]} pb</text>')
    for j,t in enumerate(tags):
        v=M[(q,t)]; x=X0+j*C
        s.append(f'<rect x="{x}" y="{y}" width="{C}" height="{C}" fill="{col(v)}" stroke="white"/>')
        s.append(f'<text x="{x+C/2}" y="{y+C/2+4}" text-anchor="middle" font-family="Helvetica,Arial" '
                 f'font-size="9.5" fill="{"white" if v>=50 else "#333"}">{v:.0f}</text>')
s.append(f'<text x="{X0}" y="{Y0+N*C+26}" font-family="Helvetica,Arial" font-size="10">'
         f'Cobertura de la secuencia de fila alineada sobre la de columna (%)</text>')
s.append('</svg>')
open(f"{B}/resultados/figuras/Fig4_cobertura.svg","w").write("\n".join(s))
print("\nescrita: resultados/figuras/Fig4_cobertura.svg")

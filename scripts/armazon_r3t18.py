#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Ultimo: armazon r3-T18 vacio vs portador del modulo oxa24."""
import re, os, glob, subprocess
from collections import defaultdict
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]
def fpath(a):
    g = glob.glob(f"{B}/datos/genomas_945/ncbi_dataset/data/{a}/*.fna")
    if g: return g[0]
    p = f"{B}/datos/ensamblados_63/{a}.fna"
    return p if os.path.exists(p) else None
def fasta(p):
    d={};k=None;buf=[]
    for l in open(p):
        if l.startswith(">"):
            if k: d[k]="".join(buf)
            k=l[1:].split()[0]; buf=[]
        else: buf.append(l.strip())
    if k: d[k]="".join(buf)
    return d

rows = tsv(f"{B}/resultados/tabla_clinica_900.tsv")
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)", rows[0][0]): rows = rows[1:]
C=["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
clin={norm(r[0]):dict(zip(C,(r+[""]*10)[:10])) for r in rows}
cr=tsv(f"{B}/resultados/apt_v2/colocalizacion.tsv"); ci={n:i for i,n in enumerate(cr[0])}
o72=defaultdict(list)
for r in cr[1:]:
    if r[ci["gen"]]=="blaOXA-72": o72[norm(r[ci["genoma"]])].append(r)

VAC=["GCA_037545805.1","GCA_053320395.1","GCA_053373735.1","GCA_053373755.1",
     "GCA_053373775.1","GCA_006493735.1","GCA_053320255.1","GCA_053373795.1"]
print("=== A. los 8 peruanos con r3-T18 vacio: portan blaOXA-72? ===")
for g in VAC:
    d=clin.get(g,{})
    if g in o72:
        for r in o72[g]:
            print(f"  {g:20s} ST{d.get('st','?'):6s} SI -> estado={r[ci['estado']]:16s} "
                  f"contig={r[ci['gen_contig']]} reps_otros=[{r[ci['reps_otros_contigs']]}]")
    else:
        print(f"  {g:20s} ST{d.get('st','?'):6s} NO porta blaOXA-72   genes=[{d.get('genes','?')}]")

SEQ=[("ANDES_con_gen","GCA_053320315.1","JBROGJ010000074.1"),
     ("ANDES_con_gen2","GCA_051942065.1","JBFECW010000069.1"),
     ("ECUADOR_con_gen","GCA_029955785.1","JASBHG010000073.1"),
     ("PERU_vacio_7851","GCA_053320395.1","JBROGN010000054.1"),
     ("PERU_vacio_7928","GCA_053373735.1","JBRZLA010000060.1"),
     ("PERU_vacio_8149","GCA_037545805.1","JARBXK010000089.1"),
     ("CHILE_CP076810","GCA_024139015.1","CP076810.1"),
     ("CHILE_8356","GCA_018156055.1","JAGTAK010000059.1"),
     ("CHILE_8329","GCA_024151615.1","JAHICI010000063.1"),
     ("COLOMBIA_5372","ERR17057656","contig00027"),
     ("VENEZUELA_5261","GCA_016485725.1","DADAKP010000051.1")]
os.makedirs(f"{B}/tmp_plas", exist_ok=True)
fa=f"{B}/tmp_plas/armazon_r3t18.fa"
with open(fa,"w") as out:
    n=0
    for tag,acc,ctg in SEQ:
        p=fpath(acc)
        if not p: print(f"  [!] sin FASTA: {acc}"); continue
        F=fasta(p)
        if ctg not in F: print(f"  [!] contig ausente: {acc}/{ctg}"); continue
        out.write(f">{tag}|{acc}|{len(F[ctg])}\n{F[ctg]}\n"); n+=1
print(f"\n=== B. secuencias extraidas: {n} -> {fa} ===")

print("\n=== C. BLASTN todos contra todos ===")
r=subprocess.run(["blastn","-query",fa,"-subject",fa,"-evalue","1e-5",
                  "-outfmt","6 qseqid sseqid pident length qlen slen"],
                 capture_output=True,text=True)
cov=defaultdict(int); ident=defaultdict(float); ql={}
for l in r.stdout.strip().split("\n"):
    if not l: continue
    f_=l.split("\t")
    if f_[0]==f_[1]: continue
    k=(f_[0],f_[1]); cov[k]+=int(f_[3]); ql[f_[0]]=int(f_[4])
    ident[k]=max(ident[k],float(f_[2]))
tags=[s.split("|")[0] for s in ql] or [t for t,_,_ in SEQ]
print(f"  {'query':18s} {'subject':18s} {'cov%':>6s} {'id%':>7s}")
for k in sorted(cov, key=lambda x:(-cov[x]/max(ql.get(x[0],1),1))):
    c=100*min(cov[k],ql[k[0]])/ql[k[0]]
    if c>=10:
        print(f"  {k[0].split('|')[0]:18s} {k[1].split('|')[0]:18s} {c:6.1f} {ident[k]:7.2f}")
if r.stderr.strip(): print("  stderr:", r.stderr.strip()[:300])

print("\n=== D. blaOXA-72 y sitios pdif en cada secuencia ===")
S=fasta(fa)
AMR=glob.glob(f"{os.path.expanduser('~')}/miniforge3/envs/abaumannii/share/amrfinderplus/data/*/AMR_CDS*")
print(f"  catalogo AMR_CDS: {AMR[0] if AMR else 'NO ENCONTRADO'}")
for k,s in S.items():
    xd=len(re.findall(r"ATTTAACATAA......TTATACGAAAT",s))+len(re.findall(r"ATTTCGTATAA......TTATGTTAAAT",s))
    print(f"  {k.split('|')[0]:18s} len={len(s):6d}  sitios_pdif_exactos={xd}")

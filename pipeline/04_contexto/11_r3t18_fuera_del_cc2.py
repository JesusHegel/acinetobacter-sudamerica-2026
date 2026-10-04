#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Punto 51: r3-T18 fuera del CC2 andino; co-localizacion con blaOXA-58 en Chile."""
import re, os, glob, subprocess
from collections import Counter, defaultdict
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]
def fpath(a):
    g = glob.glob(f"{B}/datos/genomas_945/ncbi_dataset/data/{a}/*.fna")
    return g[0] if g else (f"{B}/datos/ensamblados_63/{a}.fna"
                           if os.path.exists(f"{B}/datos/ensamblados_63/{a}.fna") else None)
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
t18=[g for g,d in clin.items() if "r3-T18" in [x.strip() for x in d["rep"].split(",")]]

hr=tsv(f"{B}/resultados/apt_v2/apt_hits_v2.tsv"); hi={n:i for i,n in enumerate(hr[0])}
t18ctg=defaultdict(list)
for r in hr[1:]:
    if r[hi["tipo_rep"]]=="r3-T18":
        t18ctg[norm(r[hi["genoma"]])].append((r[hi["contig"]], r[hi["contig_len"]]))

cr=tsv(f"{B}/resultados/apt_v2/colocalizacion.tsv"); ci={n:i for i,n in enumerate(cr[0])}
coloc=defaultdict(list)
for r in cr[1:]: coloc[norm(r[ci["genoma"]])].append(r)

print("=== A. gen co-localizado EN EL CONTIG de r3-T18 ===")
print(f"  {'genoma':20s} {'pais':10s} {'ST':7s} {'contig':22s} {'len':>8s}  gen")
res=Counter(); chile_ctg=[]
for g in sorted(t18, key=lambda x:(clin[x]['pais'], clin[x]['st'], x)):
    d=clin[g]
    for ctg,cl in (t18ctg.get(g) or [("no_hallado","?")]):
        gs=[r[ci["gen"]] for r in coloc.get(g,[])
            if r[ci["gen_contig"]]==ctg
            and "r3-T18" in [x.strip() for x in r[ci["reps_mismo_contig"]].split(",")]]
        tag=",".join(sorted(set(gs))) if gs else "-"
        res[(d['pais'],tag)]+=1
        print(f"  {g:20s} {d['pais']:10s} ST{d['st']:5s} {ctg:22s} {cl:>8s}  {tag}")
        if d['pais']=="Chile" and gs: chile_ctg.append((g,ctg,cl))

print("\n=== B. resumen pais x gen ===")
for (p,t),v in sorted(res.items()): print(f"  {p:12s} {t:20s} {v}")

print(f"\n=== C. contigs chilenos con gen co-localizado en r3-T18: {len(chile_ctg)} ===")
if chile_ctg:
    ref=("GCA_051942065.1","JBFECW010000069.1")
    tmp=f"{B}/tmp_plas"; os.makedirs(tmp, exist_ok=True)
    fa=f"{tmp}/r3t18_chile_vs_andes.fa"
    with open(fa,"w") as out:
        n=0
        for acc,ctg,cl in chile_ctg:
            p=fpath(acc)
            if p and ctg in fasta(p):
                out.write(f">CHILE_{acc}_{ctg}\n{fasta(p)[ctg]}\n"); n+=1
        p=fpath(ref[0]); F=fasta(p) if p else {}
        if ref[1] in F: out.write(f">ANDES_{ref[0]}_{ref[1]}\n{F[ref[1]]}\n"); n+=1
    print(f"  secuencias escritas: {n} -> {fa}")
    print("\n  --- BLASTN todos contra todos (cobertura e identidad) ---")
    r=subprocess.run(["blastn","-query",fa,"-subject",fa,"-evalue","1e-5",
                      "-outfmt","6 qseqid sseqid pident length qlen qcovs"],
                     capture_output=True, text=True)
    seen=set()
    for l in r.stdout.strip().split("\n"):
        if not l: continue
        f_=l.split("\t")
        if f_[0]==f_[1] or (f_[0],f_[1]) in seen: continue
        seen.add((f_[0],f_[1]))
        print(f"    {f_[0][:34]:34s} vs {f_[1][:34]:34s} id={f_[2]:>7s} cov={f_[5]:>4s}% qlen={f_[4]}")
    if r.stderr.strip(): print("  stderr:", r.stderr.strip()[:200])
else:
    print("  ninguno: en Chile r3-T18 NO co-localiza con ningun gen de resistencia")

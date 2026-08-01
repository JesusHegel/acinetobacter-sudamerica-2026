#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura 6: modulo oxa24 en el elemento andino, el brasileno y el armazon vacio."""
import os, re, glob, subprocess
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

PAT=[("C>D","ATTTCGTATAA","TTATGTTAAAT"),("D>C","ATTTAACATAA","TTATACGAAAT")]
TOL=1
def pdif(seq):
    n=len(seq); ext=seq+seq[:40]; out=[]
    for i in range(n):
        v=ext[i:i+28]
        if len(v)<28: break
        for ori,izq,der in PAT:
            if mm(v[:11],izq)+mm(v[17:],der)<=TOL:
                out.append((i,ori,v[11:17])); break
    return sorted(out)

cr=[l.rstrip("\n").split("\t") for l in open(f"{B}/resultados/apt_v2/colocalizacion.tsv") if l.strip()]
ci={n:i for i,n in enumerate(cr[0])}
def ctg_oxa(acc):
    for r in cr[1:]:
        if r[ci["genoma"]]==acc and r[ci["gen"]]=="blaOXA-72" and r[ci["estado"]]=="COLOCALIZADO":
            return r[ci["gen_contig"]]

AMR=glob.glob(f"{os.path.expanduser('~')}/miniforge3/envs/*/share/amrfinderplus/data/*/AMR_CDS.fa")
os.makedirs(f"{B}/tmp_fig6", exist_ok=True)
q=f"{B}/tmp_fig6/genes.fa"
if AMR:
    Ac=fasta(AMR[0])
    with open(q,"w") as o:
        for pat,nm in [(r"blaOXA-72\b","blaOXA-72"),(r"mph\(E\)|\bmphE\b","mph(E)"),
                       (r"msr\(E\)|\bmsrE\b","msr(E)")]:
            for k,v in Ac.items():
                if re.search(pat,k): o.write(f">{nm}\n{v}\n"); break

ELEM=[("Elemento andino  r3-T18  (Peru, ST2)","GCA_051942065.1","JBFECW010000069.1"),
      ("Armazon chileno  r3-T18  (ST109, sin modulo)","GCA_024139015.1","CP076810.1"),
      ("Elemento brasileno  r3-T1  (ST79)","GCA_009829545.1",None)]
DATA=[]
for tag,acc,ctg in ELEM:
    c = ctg or ctg_oxa(acc)
    p=fpath(acc); F=fasta(p) if p else {}
    if c not in F: print(f"  [!] no hallado: {tag}"); continue
    seq=F[c]; L=len(seq); sit=pdif(seq); gen=[]
    if AMR:
        s=f"{B}/tmp_fig6/{acc}.fa"; open(s,"w").write(f">s\n{seq}\n")
        r=subprocess.run(["blastn","-query",q,"-subject",s,"-evalue","1e-20",
                          "-perc_identity","90","-outfmt","6 qseqid pident sstart send"],
                         capture_output=True,text=True)
        seen=set()
        for l in r.stdout.strip().split("\n"):
            if not l: continue
            f=l.split("\t")
            if f[0] in seen: continue
            seen.add(f[0]); gen.append((f[0],int(f[2]),int(f[3])))
    DATA.append((tag,L,sit,gen))
    print(f"{tag}\n  {L} pb | {len(sit)} pdif | genes: {[g[0] for g in gen]}")
    for pos,ori,esp in sit: print(f"    pdif {ori} pos {pos:6d}  espaciador {esp}")

MXL=max(d[1] for d in DATA)
W=980; H=90+len(DATA)*135
X0=60; BW=700
GC={"blaOXA-72":"#b2182b","mph(E)":"#4393c3","msr(E)":"#2166ac"}
s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
   f'<rect width="{W}" height="{H}" fill="white"/>']
y=58
for tag,L,sit,gen in DATA:
    bw=BW*L/MXL
    s.append(f'<text x="{X0}" y="{y-18}" font-family="Helvetica,Arial" font-size="11.5" font-weight="600">{tag}</text>')
    s.append(f'<text x="{X0+bw+10:.0f}" y="{y+17}" font-family="Helvetica,Arial" font-size="9" fill="#666">{L} pb</text>')
    s.append(f'<rect x="{X0}" y="{y}" width="{bw:.1f}" height="24" fill="#f2f2f2" stroke="#999" stroke-width="1"/>')
    for nm,a,b in gen:
        x1=X0+bw*min(a,b)/L; x2=X0+bw*max(a,b)/L
        s.append(f'<rect x="{x1:.1f}" y="{y+3}" width="{max(x2-x1,3):.1f}" height="18" fill="{GC.get(nm,"#888")}" rx="2"/>')
        s.append(f'<text x="{(x1+x2)/2:.1f}" y="{y-4}" text-anchor="middle" font-family="Helvetica,Arial" '
                 f'font-size="8.5" fill="{GC.get(nm,"#888")}" font-weight="600">{nm}</text>')
    for pos,ori,esp in sit:
        x=X0+bw*pos/L
        s.append(f'<polygon points="{x-5:.1f},{y+24} {x+5:.1f},{y+24} {x:.1f},{y+36}" fill="#e08214"/>')
        s.append(f'<text x="{x:.1f}" y="{y+48}" text-anchor="middle" font-family="Helvetica,Arial" font-size="8">{esp}</text>')
        s.append(f'<text x="{x:.1f}" y="{y+58}" text-anchor="middle" font-family="Helvetica,Arial" font-size="7.5" fill="#777">{ori}</text>')
    y+=135
s.append(f'<polygon points="{X0},{H-26} {X0+10},{H-26} {X0+5},{H-16}" fill="#e08214"/>')
s.append(f'<text x="{X0+16}" y="{H-17}" font-family="Helvetica,Arial" font-size="9.5">'
         f'sitio pdif (XerC/XerD) con su espaciador central y orientacion; hasta 1 desajuste en los brazos</text>')
s.append(f'<text x="{X0+560}" y="{H-17}" font-family="Helvetica,Arial" font-size="9.5" fill="#666">'
         f'escala comun; posicion 0 = inicio del contig</text>')
s.append('</svg>')
open(f"{B}/resultados/figuras/Fig6_modulo.svg","w").write("\n".join(s))
print(f"\nescrita: resultados/figuras/Fig6_modulo.svg")

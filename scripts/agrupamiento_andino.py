#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Agrupamiento de enlace simple del conjunto andino. Puntos 1, 2, 3, 42."""
import re, os
from collections import Counter, defaultdict

B = os.path.expanduser("~/abaumannii")
MATRIX = f"{B}/resultados/cgmlst_st2_eval/distance_matrix_symmetric.tsv"
COLOC  = f"{B}/resultados/apt_v2/colocalizacion.tsv"
CLIN   = f"{B}/resultados/tabla_clinica_900.tsv"
BIOP   = f"{B}/datos/acc_bioproject_944.tsv"
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()
def tsv(p):
    out=[]
    with open(p) as f:
        for l in f:
            l=l.rstrip("\n")
            if l.strip(): out.append(l.split("\t"))
    return out

# --- tabla clinica ---
C = ["accesion","pais","st","anio","carb","genes","st_ox","kl","ocl","rep"]
rows = tsv(CLIN)
if not re.match(r"^(GC[AF]_|ERR|SRR|DRR)", rows[0][0]): rows = rows[1:]
clin = {}
for r in rows:
    r = r + [""]*(10-len(r)); clin[norm(r[0])] = dict(zip(C, r[:10]))
print(f"tabla clinica: {len(clin)} genomas\n")

# --- bioproject ---
bp = {}
for r in tsv(BIOP):
    if len(r) >= 2 and r[1].startswith("PRJ"): bp[norm(r[0])] = r[1]

# --- colocalizacion ---
cr = tsv(COLOC); ch = cr[0]; ci = {n:i for i,n in enumerate(ch)}
o72 = [r for r in cr[1:] if r[ci["gen"]] == "blaOXA-72"]
print("=== A. blaOXA-72 en colocalizacion.tsv ===")
print(f"filas totales: {len(o72)}   genomas unicos: {len(set(norm(r[ci['genoma']]) for r in o72))}")
print("estado:", dict(Counter(r[ci["estado"]] for r in o72)))
COL = [r for r in o72 if r[ci["reps_mismo_contig"]] not in ("-","","NA")]
gcol = sorted(set(norm(r[ci["genoma"]]) for r in COL))
print(f"con rep en el mismo contig -> {len(COL)} filas, {len(gcol)} genomas unicos  (manuscrito: 54)")
print("por pais:", dict(Counter(clin.get(g,{}).get("pais","?") for g in gcol)))

print("\n=== B. contig_len de los co-localizados (punto 42) ===")
na = Counter()
for r in COL:
    g = norm(r[ci["genoma"]])
    if r[ci["contig_len"]] in ("NA","","-"): na[clin.get(g,{}).get("pais","?")] += 1
print("filas con contig_len NA por pais:", dict(na) if na else "ninguna")
print("IDs no-GCA entre los co-localizados:", [g for g in gcol if not g.startswith("GC")] or "ninguno")

# --- subconjunto andino CC2 ---
AND = {"Peru","Perú","Ecuador"}
cc2 = sorted(g for g in gcol if clin.get(g,{}).get("pais") in AND
             and clin[g]["st"] in ("2","2724"))
print(f"\n=== C. subconjunto andino CC2 con blaOXA-72 co-localizado: {len(cc2)} genomas ===")
for g in cc2:
    d = clin[g]; print(f"  {g:22s} {d['pais']:8s} ST{d['st']:5s} {d['anio']:5s} {d['kl']:6s} {bp.get(g,'?'):14s} {d['rep']}")

# --- matriz ---
with open(MATRIX) as f:
    names = f.readline().rstrip("\n").split("\t")[1:]
    D = {}
    for l in f:
        p = l.rstrip("\n").split("\t"); D[p[0]] = p[1:]
idx = {n:i for i,n in enumerate(names)}
m2n = {norm(n): n for n in names}
def dist(a,b): return int(float(D[a][idx[b]]))
print(f"\nmatriz ST2: {len(names)} genomas")

andes_st2 = sorted(g for g,d in clin.items() if d["pais"] in AND and d["st"]=="2" and norm(g) in m2n)
print(f"TODOS los ST2 andinos en la matriz: {len(andes_st2)}   <- denominador de la Tabla 4")
inm  = [g for g in cc2 if g in m2n]
fuera= [g for g in cc2 if g not in m2n]
print(f"del subconjunto CC2 portador: {len(inm)} en la matriz, {len(fuera)} fuera -> {fuera}")

N = [m2n[g] for g in inm]
acc = {m2n[g]: g for g in inm}

# --- V13: distancia minima de cada peruano al bloque ecuatoriano ---
ec = [n for n in N if clin[acc[n]]["pais"]=="Ecuador"]
pe = [n for n in N if clin[acc[n]]["pais"] in ("Peru","Perú")]
print(f"\n=== D. distancia minima al bloque ecuatoriano (n_ec={len(ec)}, n_pe={len(pe)}) ===")
for n in sorted(pe, key=lambda x: min(dist(x,e) for e in ec)):
    print(f"  {acc[n]:22s} min={min(dist(n,e) for e in ec):4d}  {clin[acc[n]]['anio']}  {clin[acc[n]]['kl']}")
if len(ec)>1:
    print(f"  [Ecuador entre si: max={max(dist(a,b) for a in ec for b in ec if a!=b)}]")

# --- enlace simple ---
def clusters(nodes, t):
    par = {n:n for n in nodes}
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i,a in enumerate(nodes):
        for b in nodes[i+1:]:
            if dist(a,b) <= t: par[f(a)] = f(b)
    g = defaultdict(list)
    for n in nodes: g[f(n)].append(n)
    return sorted(g.values(), key=len, reverse=True)

print("\n=== E. numero de grupos por umbral ===")
print("umbral  grupos")
for t in range(0,26):
    print(f"{t:4d}   {len(clusters(N,t)):4d}")

for t in (5,8,10):
    cl = clusters(N,t)
    print(f"\n=== F. composicion a umbral {t}: {len(cl)} grupos ===")
    for i,gr in enumerate(cl,1):
        print(f"  grupo {i} (n={len(gr)}):")
        for n in sorted(gr, key=lambda x: acc[x]):
            d = clin[acc[n]]
            print(f"     {acc[n]:22s} {d['pais']:8s} {d['anio']:5s} {d['kl']:6s} {bp.get(acc[n],'?')}")

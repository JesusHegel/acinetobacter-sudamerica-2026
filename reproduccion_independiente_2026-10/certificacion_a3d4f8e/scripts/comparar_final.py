#!/usr/bin/env python3
"""Compara las salidas de la reejecucion con las publicadas en el repositorio."""
import os, re
B = os.environ["ABAUMANNII_BASE"]
R, M = f"{B}/repo/resultados", f"{B}/resultados"
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x or ""); return m.group(1) if m else (x or "").strip()
def filas(p):
    return [l.rstrip("\r\n").split("\t")
            for l in open(p, encoding="utf-8", errors="replace") if l.strip()]
print("=== tablas, comparacion fila a fila ===")
for t in ("TablaS1_900_genomas.tsv", "TablaS2_proyectos.tsv", "TablaS3_peru.tsv",
          "TablaS4_evaluabilidad.tsv", "TablaS5_ensamblados_propios.tsv",
          "tabla_clinica_900.tsv", "anexo_isaba1_oxa23.tsv"):
    a, b = f"{M}/{t}", f"{R}/{t}"
    if not (os.path.exists(a) and os.path.exists(b)):
        print(f"  {t:32s} no comparable"); continue
    A, C = filas(a), filas(b)
    if len(A) != len(C):
        print(f"  {t:32s} FILAS {len(A)} vs {len(C)}"); continue
    d = [i for i, (x, y) in enumerate(zip(A, C)) if x != y]
    est = "IDENTICA" if not d else f"difieren {len(d)} filas"
    print(f"  {t:32s} {len(A)-1:4d} filas   {est}")
print()
print("=== matriz de distancias de los ST2 ===")
def matriz(p):
    with open(p) as f:
        nom = [norm(x) for x in f.readline().rstrip("\r\n").split("\t")[1:]]
        d = {}
        for l in f:
            c = l.rstrip("\r\n").split("\t")
            for b_, v in zip(nom, c[1:]):
                d[frozenset((norm(c[0]), b_))] = int(float(v))
    return nom, d
na, da = matriz(f"{M}/cgmlst_st2_eval/distance_matrix_symmetric.tsv")
nb, db = matriz(f"{R}/cgmlst_st2_eval/distance_matrix_symmetric.tsv")
print(f"  genomas {len(na)} vs {len(nb)}   mismo conjunto: {set(na)==set(nb)}")
print(f"  pares {len(da)}   IDENTICOS: {da==db}")
if da != db:
    print(f"  pares que difieren: {sum(1 for k in da if da.get(k)!=db.get(k))}")
print()
print("=== arboles: comparacion por particiones ===")
from Bio import Phylo
def part(p):
    arb = Phylo.read(p, "newick")
    h = {norm(t.name) for t in arb.get_terminals()}
    return h, {frozenset(norm(t.name) for t in c.get_terminals())
               for c in arb.get_nonterminals()
               if 1 < len(c.get_terminals()) < len(h)}
for nom, rel in (("ST2 (213)", "grapetree/st2_msa.nwk"),
                 ("900 genomas", "grapetree_900/mst_900.nwk")):
    ha, pa = part(f"{M}/{rel}")
    hb, pb = part(f"{R}/{rel}")
    comp = len(pa & pb)
    tot = len(pa | pb)
    print(f"  {nom:12s} hojas {len(ha)}/{len(hb)} iguales: {ha==hb}   "
          f"particiones {len(pa)}/{len(pb)}  compartidas {comp}  "
          f"({100*comp/max(tot,1):.1f} %)")

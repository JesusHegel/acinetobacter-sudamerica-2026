#!/usr/bin/env python3
"""Por que difieren los perfiles de los 900 y el arbol de los ST2."""
import re, os
from collections import Counter
from Bio import Phylo
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x or ""); return m.group(1) if m else (x or "").strip()
def perfiles(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\r\n").split("\t")[1:]
        return hdr, {norm(l.split("\t")[0]): l.rstrip("\r\n").split("\t")[1:]
                     for l in f}
print("=== A. en que difieren los perfiles de los 900 ===")
ha, pa = perfiles(f"{B}/resultados/cgmlst_900_eval/cgMLST95.tsv")
hb, pb = perfiles(f"{B}/repo/resultados/cgmlst_900_eval/cgMLST95.tsv")
print(f"  mismos loci en el mismo orden: {ha == hb}")
patron = Counter(); porgen = {}; celdas = 0
for g in pa:
    n = 0
    for i, (x, y) in enumerate(zip(pa[g], pb[g])):
        if x == y: continue
        n += 1; celdas += 1
        cx = "0" if x in ("0", "", "-") else "n"
        cy = "0" if y in ("0", "", "-") else "n"
        patron[f"mio={cx} publicado={cy}"] += 1
    if n: porgen[g] = n
print(f"  genomas que difieren: {len(porgen)} de {len(pa)}")
print(f"  celdas que difieren: {celdas} de {len(pa)*len(ha)}"
      f"  ({100*celdas/(len(pa)*len(ha)):.3f} %)")
print(f"  loci distintos por genoma: min={min(porgen.values())} "
      f"max={max(porgen.values())} mediana="
      f"{sorted(porgen.values())[len(porgen)//2]}")
print(f"  patron de la diferencia: {dict(patron)}")
loci = Counter()
for g, _ in porgen.items():
    for i, (x, y) in enumerate(zip(pa[g], pb[g])):
        if x != y: loci[ha[i]] += 1
print(f"  loci implicados: {len(loci)}   los 5 mas frecuentes: {loci.most_common(5)}")
print("\n=== B. el arbol de los ST2 ===")
def part(p):
    arb = Phylo.read(p, "newick")
    h = {norm(t.name) for t in arb.get_terminals()}
    return h, {frozenset(norm(t.name) for t in c.get_terminals())
               for c in arb.get_nonterminals()
               if 1 < len(c.get_terminals()) < len(h)}
ja, qa = part(f"{B}/resultados/grapetree/st2_msa.nwk")
jb, qb = part(f"{B}/repo/resultados/grapetree/st2_msa.nwk")
print(f"  hojas {len(ja)} vs {len(jb)}   mismo conjunto: {ja == jb}")
print(f"  particiones {len(qa)} vs {len(qb)}   compartidas {len(qa & qb)}")
print(f"  coincidencia: {100*len(qa & qb)/max(len(qa | qb), 1):.1f} %")

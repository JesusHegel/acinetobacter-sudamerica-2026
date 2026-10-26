#!/usr/bin/env python3
"""Compara dos arboles newick: hojas, particiones y el clado andino."""
import re, os, sys
from Bio import Phylo
B = os.path.expanduser("~/abaumannii")
MIO = f"{B}/resultados/grapetree_900/mst_900.nwk"
PUB = f"{B}/repo/resultados/grapetree_900/mst_900.nwk"
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x or ""); return m.group(1) if m else (x or "").strip()
def carga(p):
    arb = Phylo.read(p, "newick")
    hojas = {norm(t.name) for t in arb.get_terminals()}
    part = set()
    for cl in arb.get_nonterminals():
        s = frozenset(norm(t.name) for t in cl.get_terminals())
        if 1 < len(s) < len(hojas): part.add(s)
    return hojas, part
ha, pa = carga(MIO)
hb, pb = carga(PUB)
print(f"hojas: mio={len(ha)}  publicado={len(hb)}   mismo conjunto: {ha == hb}")
if ha != hb:
    print(f"  solo en el mio: {sorted(ha - hb)[:5]}")
    print(f"  solo en el publicado: {sorted(hb - ha)[:5]}")
print(f"particiones: mio={len(pa)}  publicado={len(pb)}  compartidas={len(pa & pb)}")
print(f"  solo mias={len(pa - pb)}   solo suyas={len(pb - pa)}")
print(f"  coincidencia: {100*len(pa & pb)/max(len(pa | pb), 1):.1f} %")
AND = {"GCA_026311855.1","GCA_026311875.1","GCA_026340685.1","GCA_026340785.1",
       "GCA_026340935.1","GCA_029955785.1","GCA_030016875.1","GCA_037545705.1",
       "GCA_037545745.1","GCA_037545755.1","GCA_037545805.1","GCA_051941745.1",
       "GCA_051942065.1","GCA_053320195.1","GCA_053320295.1","GCA_053320315.1",
       "GCA_053320335.1","GCA_053320355.1","GCA_053320375.1","GCA_053320395.1",
       "GCA_053320415.1","GCA_053373735.1","GCA_053373755.1","GCA_053373775.1"}
def minimo(part, obj):
    c = [s for s in part if obj <= s]
    return min(c, key=len) if c else None
for et, part in (("mio", pa), ("publicado", pb)):
    m = minimo(part, AND)
    print(f"clado minimo que contiene los 24 portadores andinos [{et}]: "
          f"{len(m) if m else 'NO HALLADO'} genomas")
ma, mb = minimo(pa, AND), minimo(pb, AND)
print(f"ese clado es el mismo conjunto: {ma == mb}")

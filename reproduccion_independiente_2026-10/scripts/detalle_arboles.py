#!/usr/bin/env python3
"""Detalle de la comparacion: los extras del clado andino y las
particiones que difieren entre los dos arboles."""
import re, os, csv
from collections import Counter
from Bio import Phylo
B = os.path.expanduser("~/abaumannii")
MIO = f"{B}/resultados/grapetree_900/mst_900.nwk"
PUB = f"{B}/repo/resultados/grapetree_900/mst_900.nwk"
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x or ""); return m.group(1) if m else (x or "").strip()
S1 = {d["accesion"]: d for d in csv.DictReader(
      open(f"{B}/resultados/TablaS1_900_genomas.tsv"), delimiter="\t")}
def et(g):
    d = S1.get(g, {})
    return f"{g} {d.get('pais','?')} ST{d.get('st_pasteur','?')} {d.get('anio','?')}"
def carga(p):
    arb = Phylo.read(p, "newick")
    hojas = {norm(t.name) for t in arb.get_terminals()}
    part = set()
    for cl in arb.get_nonterminals():
        s = frozenset(norm(t.name) for t in cl.get_terminals())
        if 1 < len(s) < len(hojas): part.add(s)
    return part
pa, pb = carga(MIO), carga(PUB)
AND = {g for g, d in S1.items()
       if d["oxa72_contexto"] not in ("", "no_evaluable")
       and d["pais"] in ("Peru", "Ecuador") and d["st_pasteur"] in ("2", "2724")}
print(f"portadores andinos del CC2: {len(AND)}")
cl = min((s for s in pa if AND <= s), key=len)
print(f"clado minimo: {len(cl)} genomas; los {len(cl)-len(AND)} que no son portadores:")
for g in sorted(cl - AND):
    print(f"   {et(g)}  contexto={S1.get(g,{}).get('oxa72_contexto','') or 'sin gen'}")
print("\nparticiones solo en el mio:")
for s in sorted(pa - pb, key=len):
    print(f"   n={len(s):3d}  {dict(Counter(S1.get(g,{}).get('pais','?') for g in s))}")
print("particiones solo en el publicado:")
for s in sorted(pb - pa, key=len):
    print(f"   n={len(s):3d}  {dict(Counter(S1.get(g,{}).get('pais','?') for g in s))}")
tocan = [s for s in (pa - pb) | (pb - pa) if s & AND]
print(f"\nparticiones discrepantes que tocan algun portador andino: {len(tocan)}")

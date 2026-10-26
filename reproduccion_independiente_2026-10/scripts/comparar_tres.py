#!/usr/bin/env python3
"""Compara tres salidas con las publicadas, normalizando los nombres."""
import re, os, hashlib
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x or ""); return m.group(1) if m else (x or "").strip()
def matriz(p):
    with open(p) as f:
        nom = [norm(x) for x in f.readline().rstrip("\r\n").split("\t")[1:]]
        d = {}
        for l in f:
            c = l.rstrip("\r\n").split("\t")
            for b, v in zip(nom, c[1:]):
                d[frozenset((norm(c[0]), b))] = int(float(v))
    return nom, d
def perfiles(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\r\n").split("\t")[1:]
        d = {}
        for l in f:
            c = l.rstrip("\r\n").split("\t")
            d[norm(c[0])] = hashlib.md5("\t".join(c[1:]).encode()).hexdigest()
    return hdr, d
def arbol(p):
    return re.sub(r"(GC[AF]_\d+\.\d+)_[A-Za-z0-9_.]*genomic", r"\1",
                  open(p).read().strip())
print("=== 1. matriz de distancias de los ST2 ===")
na, da = matriz(f"{B}/resultados/cgmlst_st2_eval/distance_matrix_symmetric.tsv")
nb, db = matriz(f"{B}/repo/resultados/cgmlst_st2_eval/distance_matrix_symmetric.tsv")
print(f"  genomas {len(na)} vs {len(nb)}   mismo conjunto: {set(na) == set(nb)}")
print(f"  pares {len(da)}   IDENTICOS: {da == db}")
if da != db:
    print(f"  pares que difieren: {sum(1 for k in da if da.get(k) != db.get(k))}")
print("=== 2. arbol de los ST2 ===")
aa, ab = (arbol(f"{B}/resultados/grapetree/st2_msa.nwk"),
          arbol(f"{B}/repo/resultados/grapetree/st2_msa.nwk"))
print(f"  longitud {len(aa)} vs {len(ab)}   IDENTICOS: {aa == ab}")
print("=== 3. perfiles cgMLST de los 900 ===")
ha, pa = perfiles(f"{B}/resultados/cgmlst_900_eval/cgMLST95.tsv")
hb, pb = perfiles(f"{B}/repo/resultados/cgmlst_900_eval/cgMLST95.tsv")
print(f"  loci {len(ha)} vs {len(hb)}   genomas {len(pa)} vs {len(pb)}")
print(f"  mismo conjunto de genomas: {set(pa) == set(pb)}")
if set(pa) == set(pb) and len(ha) == len(hb):
    dif = [g for g in pa if pa[g] != pb[g]]
    print(f"  genomas con perfil distinto: {len(dif)}   {dif[:4]}")

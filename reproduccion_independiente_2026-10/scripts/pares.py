#!/usr/bin/env python3
"""Cuantos pares de genomas cambian de distancia por los loci discrepantes."""
import re, os
from collections import Counter
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x or ""); return m.group(1) if m else (x or "").strip()
def perfiles(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\r\n").split("\t")[1:]
        return hdr, {norm(l.split("\t")[0]): l.rstrip("\r\n").split("\t")[1:]
                     for l in f}
ha, pa = perfiles(f"{B}/resultados/cgmlst_900_eval/cgMLST95.tsv")
hb, pb = perfiles(f"{B}/repo/resultados/cgmlst_900_eval/cgMLST95.tsv")
gen = sorted(pa)
afect = [i for i in range(len(ha)) if any(pa[g][i] != pb[g][i] for g in gen)]
print(f"loci afectados: {len(afect)}")
def falta(v): return v in ("0", "", "-")
cambio = {}
for i in afect:
    A = {g: pa[g][i] for g in gen}
    C = {g: pb[g][i] for g in gen}
    mov = {g for g in gen if A[g] != C[g]}
    clases_a = len({v for v in A.values() if not falta(v)})
    clases_c = len({v for v in C.values() if not falta(v)})
    pares = {frozenset((g, h)) for g in mov for h in gen if h != g}
    n = 0
    for k in pares:
        g, h = tuple(k)
        da = 0 if (falta(A[g]) or falta(A[h])) else int(A[g] != A[h])
        dc = 0 if (falta(C[g]) or falta(C[h])) else int(C[g] != C[h])
        if da != dc:
            cambio[k] = cambio.get(k, 0) + (da - dc); n += 1
    print(f"  {ha[i]}: {len(mov):3d} genomas movidos  "
          f"alelos distintos mio={clases_a:3d} suyo={clases_c:3d}  "
          f"pares afectados={n}")
tot = len(gen) * (len(gen) - 1) // 2
print(f"\npares cuya distancia cambia: {len(cambio)} de {tot} "
      f"({100*len(cambio)/tot:.4f} %)")
if cambio:
    print(f"  magnitud: {dict(Counter(cambio.values()))}")
    print(f"  cambio maximo en valor absoluto: "
          f"{max(abs(x) for x in cambio.values())} loci sobre {len(ha)}")

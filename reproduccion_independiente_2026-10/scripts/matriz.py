#!/usr/bin/env python3
"""Genera cgMLST_profiles.tsv y distance_matrix_symmetric.tsv a partir
del cgMLST100.tsv que produce ExtractCgMLST."""
import os, sys, numpy as np
B = os.path.expanduser("~/abaumannii")
EV = sys.argv[1] if len(sys.argv) > 1 else f"{B}/resultados/cgmlst_st2_eval"
ENTRADA = f"{EV}/cgMLST100.tsv"
PERFILES = f"{EV}/cgMLST_profiles.tsv"
MATRIZ   = f"{EV}/distance_matrix_symmetric.tsv"
filas = [l.rstrip("\n").split("\t") for l in open(ENTRADA) if l.strip()]
hdr, datos = filas[0], filas[1:]
nombres = [r[0] for r in datos]
print(f"perfiles: {len(nombres)} genomas x {len(hdr)-1} loci")
with open(PERFILES, "w", newline="") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in datos:
        fh.write("\t".join(r) + "\n")
M = np.array([r[1:] for r in datos], dtype=object)
n = len(nombres)
D = np.zeros((n, n), dtype=np.int32)
for i in range(n):
    D[i, i+1:] = (M[i] != M[i+1:]).sum(axis=1)
D += D.T
with open(MATRIZ, "w", newline="") as fh:
    fh.write("FILE\t" + "\t".join(nombres) + "\n")
    for i, nom in enumerate(nombres):
        fh.write(nom + "\t" + "\t".join(map(str, D[i])) + "\n")
print(f"escrito: {PERFILES}")
print(f"escrito: {MATRIZ}   ({n} x {n})")
print(f"distancia minima no nula: {D[D>0].min()}   maxima: {D.max()}")

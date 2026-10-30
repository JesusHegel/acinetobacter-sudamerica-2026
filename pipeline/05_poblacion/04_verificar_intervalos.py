#!/usr/bin/env python3
"""Verifica dos afirmaciones del manuscrito: el segundo intervalo del
conjunto andino y las distancias de los tres portadores no evaluables."""
import csv, os, re
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
def norm(x):
    m = ACC.search(x); return m.group(1) if m else x.strip()
S1 = {d["accesion"]: d for d in csv.DictReader(
      open(f"{B}/resultados/TablaS1_900_genomas.tsv"), delimiter="\t")}
AND = {"Peru", "Perú", "Ecuador"}
with open(f"{B}/resultados/cgmlst_st2_eval/distance_matrix_symmetric.tsv") as f:
    nombres = f.readline().rstrip("\n").split("\t")[1:]
    D = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in f}
idx = {n: i for i, n in enumerate(nombres)}
def dist(a, b): return int(float(D[a][idx[b]]))
def pais(n): return S1.get(norm(n), {}).get("pais", "")
def ctx(n):  return S1.get(norm(n), {}).get("oxa72_contexto", "")
port  = [n for n in nombres if pais(n) in AND and ctx(n) not in ("", "no_evaluable")]
resto = [n for n in nombres if pais(n) not in AND]
intra = [dist(a, b) for i, a in enumerate(port) for b in port[i+1:]]
inter = [dist(a, b) for a in port for b in resto]
print(f"andinos portadores: {len(port)}   resto del continente: {len(resto)}")
print(f"SEGUNDO INTERVALO: intra maxima={max(intra)}  ->  inter minima={min(inter)}")
noev = [n for n in nombres if pais(n) in AND and ctx(n) == "no_evaluable"]
print(f"\nno evaluables andinos en la matriz: {len(noev)}")
for n in noev:
    d, q = min((dist(n, p), norm(p)) for p in port)
    print(f"  {norm(n)}  {S1[norm(n)]['anio']}  a {d} loci de {q}")

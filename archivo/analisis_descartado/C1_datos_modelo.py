#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""C1: tabla de entrada del modelo multinivel, con conteos enteros."""
import os, re, csv
from collections import defaultdict, Counter
B = os.path.expanduser("~/abaumannii")
def tsv(p):
    return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]

S1 = tsv(f"{B}/resultados/TablaS1_900_genomas.tsv")
hdr = S1[0]
rows = [dict(zip(hdr, r + [""]*len(hdr))) for r in S1[1:]]
print(f"genomas: {len(rows)}")

grp = defaultdict(list)
for d in rows:
    grp[(d["pais"], d["bioproject"])].append(d)
elig = {k: v for k, v in grp.items() if len(v) >= 8}
print(f"grupos pais-proyecto con n>=8: {len(elig)}   (esperado 22)")
print(f"  paises representados: {sorted({k[0] for k in elig})}")

VARS = {
    "sin_carbapenemasa": lambda d: d["carbapenemasa_adquirida"].strip().lower() in ("no","-",""),
    "ST2":               lambda d: d["st_pasteur"] in ("2","2724"),
    "blaOXA23":          lambda d: "blaOXA-23" in d["genes_carbapenemasa"],
    "blaOXA72":          lambda d: "blaOXA-72" in d["genes_carbapenemasa"],
}

out = []
for (pais, proj), gs in sorted(elig.items()):
    for v, f in VARS.items():
        k = sum(1 for d in gs if f(d))
        out.append({"pais": pais, "proyecto": proj, "variable": v,
                    "positivos": k, "n": len(gs), "pct": round(100*k/len(gs), 1)})

p = f"{B}/resultados/modelo_entrada.tsv"
with open(p, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["pais","proyecto","variable","positivos","n","pct"], delimiter="\t")
    w.writeheader()
    for d in out: w.writerow(d)
print(f"\n=== escrita: {p}  ({len(out)} filas = {len(elig)} proyectos x {len(VARS)} variables) ===")

print("\n=== proyectos por pais (determina que paises admiten varianza propia) ===")
npj = Counter(k[0] for k in elig)
for pais, n in npj.most_common():
    marca = "" if n >= 2 else "   <-- un solo proyecto: sin varianza estimable"
    print(f"  {pais:12s} {n} proyecto(s){marca}")

print("\n=== rangos observados por pais y variable ===")
for v in VARS:
    print(f"\n  --- {v} ---")
    for pais in sorted(npj, key=lambda x: -npj[x]):
        vals = sorted(d["pct"] for d in out if d["variable"] == v and d["pais"] == pais)
        if len(vals) >= 2:
            print(f"    {pais:12s} n_proy={len(vals)}  rango {vals[0]:5.1f} - {vals[-1]:5.1f}   {vals}")

print("\n=== comprobacion frente al manuscrito ===")
chk = [("Argentina","sin_carbapenemasa","0.0, 0.0, 90.9, 100.0"),
       ("Brasil","ST2","0.0 - 100.0"),
       ("Chile","blaOXA72","0.0 - 0.0"),
       ("Paraguay","blaOXA72","0.0 - 0.0")]
for pais, v, esp in chk:
    vals = sorted(d["pct"] for d in out if d["variable"]==v and d["pais"]==pais)
    print(f"  {pais:10s} {v:18s} obtenido: {vals}   manuscrito: {esp}")

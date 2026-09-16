#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F12: coordenadas de carbapenemasas adquiridas para el cruce con ISAba1.
Mismo procedimiento que F10 (blaOXA-51-like), aplicado a los determinantes
adquiridos. Salida: un fichero de coordenadas por familia de gen.
"""
import os, glob, csv, collections

B = os.path.expanduser("~/abaumannii")
S1 = {r['accesion'].strip() for r in
      csv.DictReader(open(f"{B}/resultados/TablaS1_900_genomas.tsv"), delimiter='\t')}
print(f"genomas en la Tabla S1: {len(S1)}")

# Familias a evaluar: etiqueta -> patron en la columna 'Element name'
FAM = {
    "oxa23": "OXA-23 family",
    "oxa24": "OXA-24 family",
    "oxa58": "beta-lactamase OXA-58",
    "oxa143": "OXA-143 family",
    "ndm":   "metallo-beta-lactamase NDM",
}

coords = collections.defaultdict(list)
for f in glob.glob(f"{B}/analisis_amr4/*.tsv"):
    g = os.path.basename(f)[:-4]
    if g not in S1:
        continue
    for i, l in enumerate(open(f)):
        if i == 0:
            continue
        c = l.rstrip("\n").split("\t")
        if len(c) <= 6:
            continue
        for tag, pat in FAM.items():
            if pat in c[6]:
                # corregir prefijo de contig en los ensamblados propios
                ctg = c[1]
                if ctg.startswith(g + "_"):
                    ctg = ctg[len(g) + 1:]
                coords[tag].append((g, ctg, c[2], c[3], c[4], c[5]))

# control: toda coordenada debe tener su longitud de contig
L = set()
for l in open(f"{B}/datos/contig_len_900.tsv"):
    a = l.rstrip("\n").split("\t")
    L.add(a[0] + "|" + a[1])

print()
for tag in FAM:
    rows = sorted(coords[tag])
    gs = {r[0] for r in rows}
    sin = sum(1 for r in rows if r[0] + "|" + r[1] not in L)
    out = f"{B}/datos/{tag}_coords_900.tsv"
    with open(out, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    aviso = f"  <<< {sin} sin longitud de contig" if sin else ""
    print(f"{tag:8s} {len(rows):5d} copias en {len(gs):4d} genomas{aviso}")
print(f"\nescrito en {B}/datos/*_coords_900.tsv")

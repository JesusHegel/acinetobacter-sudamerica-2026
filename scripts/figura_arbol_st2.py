#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Arbol de expansion minima de los 213 genomas ST2 (GrapeTree MSTreeV2),
coloreado por pais y por plasmido portador de blaOXA-72."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo
import csv, os

B = os.path.expanduser("~/abaumannii")
meta = {r['ID']: r for r in csv.DictReader(
        open(f"{B}/resultados/metadatos_grapetree_st2.tsv"), delimiter='\t')}

COL_PAIS = {"Peru":"#C0392B", "Ecuador":"#E67E22", "Brasil":"#2F5597",
            "Paraguay":"#7A9CC6", "Argentina":"#5B9BD5", "Venezuela":"#9DC3E6"}
COL_PLAS = {"r3-T18":"#C0392B", "no_evaluable":"#E8A33D", "-":"#CFCFCF"}
ETI_PLAS = {"r3-T18":"r3-T18 (portador)", "no_evaluable":"no evaluable",
            "-":"sin blaOXA-72"}

def dibuja(ax, campo, colores, titulo, etiquetas=None):
    arb = Phylo.read(f"{B}/resultados/grapetree/st2_msa.nwk", "newick")
    Phylo.draw(arb, axes=ax, do_show=False, label_func=lambda c: "")
    # posiciones de las hojas tras el dibujo
    ys = {}
    for i, t in enumerate(arb.get_terminals(order="postorder"), start=1):
        ys[t.name] = i
    xs = arb.depths()
    for t in arb.get_terminals():
        m = meta.get(t.name)
        c = colores.get(m[campo], "#BBBBBB") if m else "#BBBBBB"
        ax.plot(xs.get(t, 0), ys.get(t.name, 0), 'o', color=c, ms=4.5,
                mec='white', mew=0.4, zorder=5)
    ax.set_title(titulo, loc="left", fontsize=11, fontweight="bold", color="#1F3864")
    ax.set_xlabel("diferencias alelicas"); ax.set_ylabel("")
    ax.set_yticks([])
    for s in ("top","right","left"): ax.spines[s].set_visible(False)
    lab = etiquetas or {k:k for k in colores}
    ax.legend(handles=[Line2D([0],[0], marker='o', ls='', color=v, label=lab[k], ms=6)
                       for k,v in colores.items()],
              loc="lower right", frameon=False, fontsize=8)

fig, axes = plt.subplots(1, 2, figsize=(15, 10))
dibuja(axes[0], "pais", COL_PAIS, "Por pais de origen")
dibuja(axes[1], "plasmido_oxa72", COL_PLAS,
       "Por plasmido portador de blaOXA-72", ETI_PLAS)
plt.tight_layout()
plt.savefig(f"{B}/resultados/grapetree/arbol_st2.png", dpi=250)
plt.savefig(f"{B}/resultados/grapetree/arbol_st2.svg")
print("escrito: arbol_st2.png / .svg")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paneles A y B sobre la red de expansion minima de los 900 genomas.
A: color segun la carbapenemasa adquirida.
B: color segun el replicon en que reside blaOXA-72.
Ambos paneles comparten la misma topologia."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx, csv, os, collections
from Bio import Phylo

B = os.path.expanduser("~/abaumannii")
S1 = {r['accesion'].strip(): r for r in
      csv.DictReader(open(f"{B}/resultados/TablaS1_900_genomas.tsv"), delimiter='\t')}

def acc(n):
    return "_".join(n.split("_")[:2]) if n.startswith("GCA") else n.replace(".fna","")

def carb(n):
    g = S1.get(acc(n), {}).get('genes_carbapenemasa', '')
    ctx = S1.get(acc(n), {}).get('oxa72_contexto', '').strip()
    if 'blaOXA-72' in g or ctx: return "oxa72"
    if 'blaOXA-23' in g: return "oxa23"
    if 'blaOXA-58' in g: return "oxa58"
    if 'blaNDM' in g: return "ndm"
    if g and g != 'ninguna': return "otra"
    return "ninguna"

def plas(n):
    if carb(n) != "oxa72": return "no"
    ctx = S1.get(acc(n), {}).get('oxa72_contexto', '').strip()
    if not ctx or ctx == 'no_evaluable': return "no_eval"
    if 'r3-T18' in ctx: return "r3-T18"
    if 'r3-T1' in ctx or 'r3-T14' in ctx: return "r3-T1/T14"
    return "otro"

arb = Phylo.read(f"{B}/resultados/grapetree_900/mst_900.nwk", "newick")
G = nx.Graph()
def rec(cl, p=None):
    n = cl.name or f"i{id(cl)}"
    G.add_node(n, hoja=bool(cl.name))
    if p: G.add_edge(p, n)
    for h in cl.clades: rec(h, n)
rec(arb.root)

hojas = [n for n in G if G.nodes[n]["hoja"]]
otros = [n for n in G if not G.nodes[n]["hoja"]]
pos = nx.spring_layout(G, seed=11, k=0.09, iterations=200)

CA = {"oxa23":"#2F5597","oxa72":"#C0392B","oxa58":"#27AE60",
      "ndm":"#8E44AD","otra":"#E8A33D","ninguna":"#DCDCDC"}
EA = {"oxa23":"blaOXA-23","oxa72":"blaOXA-72","oxa58":"blaOXA-58",
      "ndm":"blaNDM","otra":"otra carbapenemasa","ninguna":"sin carbapenemasa adquirida"}
CB = {"r3-T18":"#C0392B","r3-T1/T14":"#2F5597","otro":"#E8A33D",
      "no_eval":"#8E44AD","no":"#DCDCDC"}
EB = {"r3-T18":"r3-T18 (andino)","r3-T1/T14":"r3-T1 / r3-T14 (brasileno)",
      "otro":"otro replicon","no_eval":"no evaluable","no":"sin blaOXA-72"}

fig, axes = plt.subplots(1, 2, figsize=(20, 10))
for ax, fn, COL, ETI, tit, letra in [
        (axes[0], carb, CA, EA, "carbapenemasa adquirida", "A"),
        (axes[1], plas, CB, EB, "plasmido portador de blaOXA-72", "B")]:
    cnt = collections.Counter(fn(n) for n in hojas)
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#CCCCCC", width=0.5)
    fondo  = [n for n in hojas if fn(n) in ("ninguna", "no")]
    frente = [n for n in hojas if fn(n) not in ("ninguna", "no")]
    nx.draw_networkx_nodes(G, pos, nodelist=fondo, node_color="#DCDCDC",
                           node_size=26, edgecolors="white", linewidths=0.3, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=frente,
                           node_color=[COL[fn(n)] for n in frente],
                           node_size=55, edgecolors="white", linewidths=0.4, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=otros, node_color="#FFF",
                           node_size=5, edgecolors="#CCCCCC", linewidths=0.3, ax=ax)
    ax.set_title(f"{letra}. Color segun {tit}", fontsize=13, fontweight="bold",
                 color="#1F3864", loc="left")
    ax.legend(handles=[Line2D([0],[0], marker='o', ls='', color=COL[k], ms=8,
                              label=f"{ETI[k]} (n={cnt.get(k,0)})")
                       for k in COL if cnt.get(k,0) > 0],
              loc="lower right", frameon=False, fontsize=9)
    ax.axis("off")

fig.suptitle("Red de expansion minima · 900 genomas clinicos de A. baumannii · cgMLST, 2133 loci",
             fontsize=14, fontweight="bold", color="#1F3864", x=0.02, ha="left")
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(f"{B}/resultados/grapetree_900/paneles_A_B.png", dpi=250)
plt.savefig(f"{B}/resultados/grapetree_900/paneles_A_B.svg")
print("escrito: paneles_A_B.png / .svg")

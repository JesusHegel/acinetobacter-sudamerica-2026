#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Red de expansion minima de los 213 genomas ST2 (GrapeTree MSTreeV2),
dibujada con networkx. Tres versiones: por pais, por plasmido y por anio."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx, csv, os, re

B = os.path.expanduser("~/abaumannii")
meta = {r['ID']: r for r in csv.DictReader(
        open(f"{B}/resultados/metadatos_grapetree_st2.tsv"), delimiter='\t')}

# --- leer el newick como grafo (nodos + longitudes de rama) ---
from Bio import Phylo
arb = Phylo.read(f"{B}/resultados/grapetree/st2_msa.nwk", "newick")
G = nx.Graph()
def rec(cl, padre=None):
    nombre = cl.name or f"i{id(cl)}"
    G.add_node(nombre, hoja=bool(cl.name))
    if padre is not None:
        G.add_edge(padre, nombre, w=cl.branch_length or 0)
    for h in cl.clades:
        rec(h, nombre)
rec(arb.root)
print(f"nodos: {G.number_of_nodes()}  aristas: {G.number_of_edges()}")

PAL = {
 "pais": ({"Peru":"#C0392B","Ecuador":"#E67E22","Brasil":"#2F5597",
           "Paraguay":"#7A9CC6","Argentina":"#5B9BD5","Venezuela":"#9DC3E6"},
          "País de origen"),
 "plasmido_oxa72": ({"r3-T18":"#C0392B","no_evaluable":"#E8A33D","-":"#D5D5D5"},
          "Plásmido portador de blaOXA-72"),
}
ETI = {"r3-T18":"r3-T18","no_evaluable":"no evaluable","-":"sin blaOXA-72"}

pos = nx.spring_layout(G, weight=None, seed=7, iterations=250)
for campo, (cols, tit) in PAL.items():
    fig, ax = plt.subplots(figsize=(11, 9))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#B8B8B8", width=0.7)
    hojas = [n for n in G if G.nodes[n]["hoja"]]
    cs = [cols.get(meta[n][campo], "#BBBBBB") if n in meta else "#EEEEEE" for n in hojas]
    nx.draw_networkx_nodes(G, pos, nodelist=hojas, node_color=cs, node_size=90,
                           edgecolors="white", linewidths=0.6, ax=ax)
    otros = [n for n in G if not G.nodes[n]["hoja"]]
    nx.draw_networkx_nodes(G, pos, nodelist=otros, node_color="#FFFFFF",
                           node_size=12, edgecolors="#B8B8B8", linewidths=0.5, ax=ax)
    ax.set_title(f"Red de expansión mínima · 213 genomas ST2 · {tit}",
                 fontsize=12, fontweight="bold", color="#1F3864", loc="left")
    ax.legend(handles=[Line2D([0],[0],marker='o',ls='',color=v,
                              label=ETI.get(k,k), ms=8) for k,v in cols.items()],
              loc="lower right", frameon=False, fontsize=9)
    ax.axis("off"); plt.tight_layout()
    out = f"{B}/resultados/grapetree/red_st2_{campo}"
    plt.savefig(out+".png", dpi=250); plt.savefig(out+".svg"); plt.close()
    print("escrito:", os.path.basename(out)+".png/.svg")

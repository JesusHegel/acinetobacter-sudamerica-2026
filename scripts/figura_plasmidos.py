#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Mapa comparado del elemento andino r3-T18 y del armazon chileno cerrado."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord

C = {"rep": "#2F5597", "mob": "#7A9CC6", "res": "#C0392B",
     "hip": "#D9D9D9", "frag": "#F0F0F0", "pdif": "#E8A33D"}

andino = [
    (3, 92, +1, "frag", "frag N"),
    (312, 599, -1, "hip", ""),
    (596, 1546, -1, "rep", "rep r3-T18"),
    (2225, 3880, -1, "mob", "mobQ"),
    (4059, 4352, +1, "hip", ""),
    (4386, 5039, -1, "hip", ""),
    (5166, 5471, -1, "hip", ""),
    (5615, 5938, +1, "hip", ""),
    (5931, 6233, +1, "hip", ""),
    (6287, 6610, +1, "hip", ""),
    (7150, 7977, +1, "res", "blaOXA-72"),
    (8044, 8109, +1, "frag", "frag C"),
]
pdif_and = [137, 4368, 7049]

chileno = [
    (1, 951, +1, "rep", "rep r3-T18"),
    (948, 1235, +1, "hip", ""),
    (1435, 1893, -1, "hip", ""),
    (1893, 2627, -1, "hip", ""),
    (3165, 3488, -1, "hip", ""),
    (3542, 3844, -1, "hip", ""),
    (3837, 4160, -1, "hip", ""),
    (4304, 4609, +1, "hip", ""),
    (4736, 5389, +1, "hip", ""),
    (5423, 5716, -1, "hip", ""),
    (5895, 7547, +1, "mob", "mobQ"),
]
pdif_chi = [1381, 2697, 5378]

def rec(genes, L, pdifs):
    fs = [GraphicFeature(start=s, end=e, strand=d, color=C[k],
                         label=lab if lab else None, linewidth=0.6)
          for s, e, d, k, lab in genes]
    for p in pdifs:
        fs.append(GraphicFeature(start=p, end=p + 28, strand=0,
                                 color=C["pdif"], label=None, linewidth=0))
    return GraphicRecord(sequence_length=L, features=fs)

fig, axes = plt.subplots(2, 1, figsize=(13, 5.2))
rec(andino, 8111, pdif_and).plot(ax=axes[0], with_ruler=False)
axes[0].set_title("Elemento andino r3-T18  ·  8111 pb  ·  con módulo blaOXA-72",
                  loc="left", fontsize=11, fontweight="bold", color="#1F3864")
rec(chileno, 8229, pdif_chi).plot(ax=axes[1])
axes[1].set_title("Armazón chileno p3UC20804 (CP076810.1)  ·  8229 pb  ·  sin módulo",
                  loc="left", fontsize=11, fontweight="bold", color="#1F3864")

from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=C["rep"], label="rep"),
                    Patch(color=C["mob"], label="mobQ"),
                    Patch(color=C["res"], label="blaOXA-72"),
                    Patch(color=C["hip"], label="hipotética"),
                    Patch(color=C["pdif"], label="sitio pdif")],
           loc="lower center", ncol=5, frameon=False, fontsize=9)
plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.savefig("figura_plasmidos.png", dpi=300)
plt.savefig("figura_plasmidos.svg")
print("escrito: figura_plasmidos.png / .svg")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura 1: diagrama de flujo de seleccion del conjunto de datos."""
import os
B = os.path.expanduser("~/abaumannii")
os.makedirs(f"{B}/resultados/figuras", exist_ok=True)

W, H = 900, 1000
def box(x, y, w, h, txt, sub="", fill="#ffffff", stroke="#333"):
    s = f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>'
    ty = y + (h/2 if not sub else h/2 - 7)
    s += f'<text x="{x+w/2}" y="{ty}" text-anchor="middle" font-family="Helvetica,Arial" font-size="13" font-weight="600">{txt}</text>'
    if sub:
        s += f'<text x="{x+w/2}" y="{ty+15}" text-anchor="middle" font-family="Helvetica,Arial" font-size="10.5" fill="#444">{sub}</text>'
    return s
def arrow(x, y1, y2):
    return (f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2-7}" stroke="#333" stroke-width="1.2"/>'
            f'<polygon points="{x-4},{y2-7} {x+4},{y2-7} {x},{y2}" fill="#333"/>')
def excl(x, y, w, h, txt):
    s = f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" fill="#f7f0e8" stroke="#a8815a" stroke-width="1"/>'
    for i, l in enumerate(txt.split("\n")):
        s += (f'<text x="{x+10}" y="{y+16+i*13}" font-family="Helvetica,Arial" '
              f'font-size="10" fill="#5c4426">{l}</text>')
    return s

CX, BW, BX = 300, 300, 150
p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
     f'<rect width="{W}" height="{H}" fill="white"/>']

y = 20
p.append(box(BX, y, BW, 46, "1119 registros recuperados",
             "NCBI Pathogen Detection, 28/07/2026", fill="#eef2f7"))
p.append(excl(500, y+2, 340, 30, "174 sin ensamblado depositado\n(evaluados para ensamblado propio)"))
p.append(arrow(CX, y+46, y+86)); y += 86

p.append(box(BX, y, BW, 40, "945 ensamblados descargados"))
p.append(excl(500, y-4, 340, 44,
              "20 de otras especies del complejo ACB (ANI)\n"
              "  12 A. seifertii · 5 A. pittii · 3 A. nosocomialis\n"
              "39 con contaminacion > 5 % (CheckM2)"))
p.append(arrow(CX, y+40, y+95)); y += 95

p.append(box(BX, y, BW, 40, "886 ensamblados publicos"))
p.append(arrow(CX, y+40, y+80)); y += 80

p.append(box(BX, y, BW, 52, "+ 58 ensamblados propios",
             "Shovill 1.4.2 / SPAdes 3.15.5", fill="#eaf3ea", stroke="#4a7c4a"))
p.append(excl(500, y-2, 340, 56,
              "63 ensamblados desde lecturas crudas:\n"
              "  Colombia 54 · Peru 4 · Argentina 3 · Bolivia 1 · Chile 1\n"
              "menos 5 A. nosocomialis (mismo programa colombiano)\n"
              "5 BioProjects de origen"))
p.append(arrow(CX, y+52, y+100)); y += 100

p.append(box(BX, y, BW, 40, "944 genomas tipificados"))
p.append(excl(500, y-4, 340, 44,
              "44 aislados no clinicos:\n"
              "  21 ambientales · 23 de ambiente hospitalario\n"
              "  Brasil 26 · Argentina 18"))
p.append(arrow(CX, y+40, y+95)); y += 95

p.append(box(BX, y, BW, 56, "900 genomas clinicos",
             "9 paises · 104 BioProjects · 360 eventos", fill="#dde7f2", stroke="#35506e"))
y += 90

p.append(f'<text x="{CX}" y="{y}" text-anchor="middle" font-family="Helvetica,Arial" '
         f'font-size="11.5" font-weight="600">Subconjuntos de analisis</text>')
y += 16
for txt, sub in [("22 BioProjects con n &#8805; 8",
                  "176 observaciones sobre 8 variables (seccion 3.3)"),
                 ("213 genomas ST2 &#183; 149 ST79",
                  "cgMLST, 2390 loci (seccion 3.9)"),
                 ("73 portadores de bla-OXA-72",
                  "63 con replicon co-localizado (secciones 3.5 y 3.6)")]:
    p.append(box(90, y, 420, 38, txt, sub, fill="#fafafa"))
    y += 46

p.append('</svg>')
out = f"{B}/resultados/figuras/Fig1_flujo.svg"
open(out, "w").write("\n".join(p))
print(f"escrita: {out}")
print("\ncomprobacion de cifras:")
for a, b in [("1119-174","945"),("945-20-39","886"),("886+58","944"),("944-44","900")]:
    print(f"  {a} = {b}")

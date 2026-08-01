#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura 3: A) replicon portador de blaOXA-72 por region. B) r3-T18 con y sin gen."""
import os
from collections import Counter, defaultdict
B = os.path.expanduser("~/abaumannii")
def tsv(p):
    return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]
S1 = tsv(f"{B}/resultados/TablaS1_900_genomas.tsv")
h = S1[0]; rows = [dict(zip(h, r + [""]*len(h))) for r in S1[1:]]
AND = {"Peru","Perú","Ecuador"}

col = [d for d in rows if d["oxa72_contexto"] not in ("", "no_evaluable")]
def reg(d):
    if d["pais"] in AND and d["st_pasteur"] in ("2","2724"): return "CC2 andino"
    if d["pais"] in AND: return "Andes, otros linajes"
    return d["pais"]
A = defaultdict(Counter)
for d in col: A[reg(d)][d["oxa72_contexto"]] += 1
print("=== PANEL A: replicon portador por region ===")
for r in A:
    print(f"  {r:22s} n={sum(A[r].values()):3d}  {dict(A[r])}")

t18 = [d for d in rows if "r3-T18" in [x.strip() for x in d["tipos_rep_apt"].split(",")]]
Bp = defaultdict(lambda: [0,0])
for d in t18:
    k = (d["pais"], "ST"+d["st_pasteur"] if d["st_pasteur"] not in ("-","") else "ST no asignado")
    Bp[k][0 if "r3-T18" in d["oxa72_contexto"] else 1] += 1
print(f"\n=== PANEL B: los {len(t18)} genomas con r3-T18 ===")
for k in sorted(Bp, key=lambda x: (-sum(Bp[x]), x)):
    con, sin = Bp[k]
    print(f"  {k[0]:11s} {k[1]:16s} con gen={con:2d}  sin gen={sin:2d}")

# --- SVG ---
COL = {"r3-T18":"#2166ac","r3-T1":"#b2182b","r3-T14":"#d6604d","r3-T1,r3-T14":"#8c2d04",
       "r3-T13":"#f4a582","r3-T50":"#5aae61","r3-T8,r3-T97":"#9970ab"}
W,H = 880, 700
s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
     f'<rect width="{W}" height="{H}" fill="white"/>',
     f'<text x="20" y="26" font-family="Helvetica,Arial" font-size="15" font-weight="700">A</text>',
     f'<text x="42" y="26" font-family="Helvetica,Arial" font-size="12.5">Tipo de replicon portador de bla-OXA-72 (n = {len(col)} genomas)</text>']
ordA = sorted(A, key=lambda r: -sum(A[r].values()))
mx = max(sum(A[r].values()) for r in A); X0, BW = 190, 560; y = 48
for r in ordA:
    tot = sum(A[r].values()); x = X0
    s.append(f'<text x="{X0-10}" y="{y+15}" text-anchor="end" font-family="Helvetica,Arial" font-size="11">{r}</text>')
    for rep, n in A[r].most_common():
        w = BW * n / mx
        s.append(f'<rect x="{x}" y="{y}" width="{w:.1f}" height="22" fill="{COL.get(rep,"#999")}" stroke="white"/>')
        if w > 34:
            s.append(f'<text x="{x+w/2:.1f}" y="{y+15}" text-anchor="middle" font-family="Helvetica,Arial" font-size="9.5" fill="white" font-weight="600">{n}</text>')
        x += w
    s.append(f'<text x="{x+8:.1f}" y="{y+15}" font-family="Helvetica,Arial" font-size="10" fill="#555">n={tot}</text>')
    y += 30
y += 6
lx = X0
for rep in ["r3-T18","r3-T1,r3-T14","r3-T1","r3-T14","r3-T13","r3-T50","r3-T8,r3-T97"]:
    s.append(f'<rect x="{lx}" y="{y}" width="11" height="11" fill="{COL[rep]}"/>')
    s.append(f'<text x="{lx+15}" y="{y+10}" font-family="Helvetica,Arial" font-size="9.5">{rep}</text>')
    lx += 24 + len(rep)*5.6
y += 44

s.append(f'<text x="20" y="{y}" font-family="Helvetica,Arial" font-size="15" font-weight="700">B</text>')
s.append(f'<text x="42" y="{y}" font-family="Helvetica,Arial" font-size="12.5">Distribucion del replicon r3-T18 (n = {len(t18)} genomas)</text>')
y += 22
mxB = max(sum(v) for v in Bp.values())
for k in sorted(Bp, key=lambda x: (-sum(Bp[x]), x)):
    con, sin = Bp[k]
    s.append(f'<text x="{X0-10}" y="{y+14}" text-anchor="end" font-family="Helvetica,Arial" font-size="10.5">{k[0]} &#183; {k[1]}</text>')
    x = X0
    for n, c, lb in ((con,"#2166ac","con"),(sin,"#c9ced4","sin")):
        if n:
            w = BW * n / mxB
            s.append(f'<rect x="{x}" y="{y}" width="{w:.1f}" height="20" fill="{c}" stroke="white"/>')
            if w > 20:
                s.append(f'<text x="{x+w/2:.1f}" y="{y+14}" text-anchor="middle" font-family="Helvetica,Arial" font-size="9.5" fill="{"white" if lb=="con" else "#333"}" font-weight="600">{n}</text>')
            x += w
    y += 26
y += 8
s.append(f'<rect x="{X0}" y="{y}" width="11" height="11" fill="#2166ac"/>')
s.append(f'<text x="{X0+15}" y="{y+10}" font-family="Helvetica,Arial" font-size="9.5">con bla-OXA-72 co-localizado</text>')
s.append(f'<rect x="{X0+190}" y="{y}" width="11" height="11" fill="#c9ced4"/>')
s.append(f'<text x="{X0+205}" y="{y+10}" font-family="Helvetica,Arial" font-size="9.5">sin gen de resistencia asociado</text>')
s.append('</svg>')
out = f"{B}/resultados/figuras/Fig3_replicones.svg"
open(out,"w").write("\n".join(s))
print(f"\nescrita: {out}")

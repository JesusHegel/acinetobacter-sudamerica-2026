#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""C5: exclusividad mutua andino/brasileno sobre EVENTOS, no genomas."""
import os
from math import comb
def tsv(p):
    return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]
B = os.path.expanduser("~/abaumannii")
S1 = tsv(f"{B}/resultados/TablaS1_900_genomas.tsv")
h = S1[0]; rows = [dict(zip(h, r + [""]*len(h))) for r in S1[1:]]

AND = {"Peru","Perú","Ecuador"}
col = [d for d in rows if d["oxa72_contexto"] not in ("", "no_evaluable")]
and_t18 = sum(1 for d in col if d["pais"] in AND and "r3-T18" in d["oxa72_contexto"])
and_otr = sum(1 for d in col if d["pais"] in AND and "r3-T18" not in d["oxa72_contexto"])
bra_t18 = sum(1 for d in col if d["pais"] == "Brasil" and "r3-T18" in d["oxa72_contexto"])
bra_otr = sum(1 for d in col if d["pais"] == "Brasil" and "r3-T18" not in d["oxa72_contexto"])

def fisher(a, b, c, d):
    n = a+b+c+d
    obs = comb(a+b, a) * comb(c+d, c) / comb(n, a+c)
    p = 0.0
    for i in range(max(0, a+c-(c+d)), min(a+b, a+c)+1):
        pr = comb(a+b, i) * comb(c+d, a+c-i) / comb(n, a+c)
        if pr <= obs * (1 + 1e-9): p += pr
    return p

print("=== A. sobre GENOMAS (version inflada, NO usar como principal) ===")
print(f"          r3-T18   otro")
print(f"  Andes   {and_t18:6d} {and_otr:6d}")
print(f"  Brasil  {bra_t18:6d} {bra_otr:6d}")
print(f"  p = {fisher(and_t18, and_otr, bra_t18, bra_otr):.3g}")

print("\n=== B. sobre EVENTOS (version defendible) ===")
EV_AND = 9   # enlace completo, umbral 8, seccion 3.9
import collections
bra_ev = len({(d["pais"], d["st_pasteur"], d["bioproject"], d["anio"])
              for d in col if d["pais"] == "Brasil"})
print(f"  eventos andinos (cgMLST, enlace completo): {EV_AND}, todos r3-T18")
print(f"  eventos brasilenos (pais+ST+proyecto+anio): {bra_ev}, ninguno r3-T18")
print(f"          r3-T18   otro")
print(f"  Andes   {EV_AND:6d} {0:6d}")
print(f"  Brasil  {0:6d} {bra_ev:6d}")
print(f"  p = {fisher(EV_AND, 0, 0, bra_ev):.3g}")

print("\n=== C. salvedad obligatoria ===")
print("  Region y linaje estan confundidos: los portadores andinos son CC2 y los")
print("  brasilenos pertenecen a CC79, CC15, ST730 e IC-6. El test no distingue")
print("  'depende de la region' de 'depende del linaje'.")
lin = collections.Counter(d["st_pasteur"] for d in col if d["pais"] == "Brasil")
print(f"\n  linajes brasilenos portadores: {dict(lin)}")
lin_a = collections.Counter(d["st_pasteur"] for d in col if d["pais"] in AND)
print(f"  linajes andinos portadores   : {dict(lin_a)}")

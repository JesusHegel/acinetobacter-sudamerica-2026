#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bootstrap del cociente de varianzas intranacional/internacional (2000 replicas).
Resultado: estimador inestable, no reportado en el manuscrito (seccion 2.5.5)."""
import random, statistics as st
import os
random.seed(42)
rows=[]
for L in open(os.path.expanduser('~/abaumannii/datos/proyectos_var.tsv')):
    f=L.rstrip('\n').split('\t')
    if len(f)<7: continue
    rows.append((f[0],f[1],int(f[2]),int(f[3]),int(f[4]),int(f[5]),int(f[6])))
VARS=[('sinC',3),('ST2',4),('OXA72',5),('OXA23',6)]
def ratio(sample, idx):
    byc={}
    for r in sample:
        p=r[idx]/r[2]
        byc.setdefault(r[0],[]).append(p)
    W=[st.variance(v) for v in byc.values() if len(v)>=2]
    Wm=sum(W)/len(W) if W else 0.0
    means=[sum(v)/len(v) for v in byc.values()]
    B=st.variance(means) if len(means)>=2 else 0.0
    return Wm/B if B>0 else float('nan')
print(f"{'variable':8} {'ratio':>7} {'IC95 inf':>9} {'IC95 sup':>9}")
for name,idx in VARS:
    obs=ratio(rows,idx)
    bs=[]
    for _ in range(2000):
        s=[random.choice(rows) for _ in rows]
        r=ratio(s,idx)
        if r==r: bs.append(r)
    bs.sort()
    lo=bs[int(.025*len(bs))]; hi=bs[int(.975*len(bs))]
    print(f"{name:8} {obs:7.2f} {lo:9.2f} {hi:9.2f}")

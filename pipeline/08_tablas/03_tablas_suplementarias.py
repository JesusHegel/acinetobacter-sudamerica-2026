#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tablas suplementarias S2, S3, S4, S5."""
import os, csv
from collections import Counter, defaultdict
B = os.path.expanduser("~/abaumannii")
def tsv(p): return [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]
def wr(nm, cols, rows):
    p=f"{B}/resultados/{nm}"
    with open(p,"w",newline="",encoding="utf-8") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(cols); w.writerows(rows)
    print(f"  {nm}: {len(rows)} filas x {len(cols)} columnas")

S1=tsv(f"{B}/resultados/TablaS1_900_genomas.tsv"); h=S1[0]
R=[dict(zip(h,r+[""]*len(h))) for r in S1[1:]]

# --- S2: los 22 proyectos con las 8 variables ---
grp=defaultdict(list)
for d in R: grp[(d["pais"],d["bioproject"])].append(d)
elig={k:v for k,v in grp.items() if len(v)>=8}
def dom(gs,c):
    C=Counter(d[c] for d in gs if d[c] not in ("-",""))
    if not C: return "", 0.0
    k,n=C.most_common(1)[0]; return k, round(100*n/len(gs),1)
rows=[]
for (pais,proj),gs in sorted(elig.items()):
    n=len(gs)
    p=lambda f: round(100*sum(1 for d in gs if f(d))/n,1)
    kl,klp=dom(gs,"kl"); st,stp=dom(gs,"st_pasteur")
    rep=Counter(t.strip() for d in gs for t in d["tipos_rep_apt"].split(",")
                if t.strip() not in ("-","","ninguno"))
    rt,rn = rep.most_common(1)[0] if rep else ("",0)
    ct=[int(d["n_contigs"]) for d in gs if str(d["n_contigs"]).isdigit()]
    anios=sorted({d["anio"] for d in gs if d["anio"] not in ("","NA")})
    rows.append([pais,proj,n,
        p(lambda d: d["carbapenemasa_adquirida"].lower() in ("no","-","")),
        p(lambda d: d["st_pasteur"] in ("2","2724")),
        p(lambda d: "blaOXA-23" in d["genes_carbapenemasa"]),
        p(lambda d: "blaOXA-72" in d["genes_carbapenemasa"]),
        kl, klp, "ST"+st if st else "", stp,
        len({d["st_pasteur"] for d in gs if d["st_pasteur"] not in ("-","")}),
        rt, round(100*rn/n,1),
        round(sum(ct)/len(ct),1) if ct else "",
        f"{anios[0]}-{anios[-1]}" if anios else ""])
wr("TablaS2_proyectos.tsv",
   ["pais","bioproject","n_genomas","pct_sin_carbapenemasa","pct_ST2","pct_blaOXA23",
    "pct_blaOXA72","kl_dominante","pct_kl_dominante","st_dominante","pct_st_dominante",
    "n_linajes","replicon_dominante","pct_replicon_dominante","media_contigs","periodo"],
   sorted(rows,key=lambda r:(r[0],-r[2])))

# --- S3: conjunto peruano ---
pe=[d for d in R if d["pais"]=="Peru"]
wr("TablaS3_peru.tsv",
   ["accesion","bioproject","anio","st_pasteur","kl","ocl","genes_carbapenemasa",
    "tipos_rep_apt","oxa72_contexto","n_contigs","fuente_aislamiento"],
   [[d["accesion"],d["bioproject"],d["anio"],"ST"+d["st_pasteur"] if d["st_pasteur"] not in ("-","") else "-",
     d["kl"],d["ocl"],d["genes_carbapenemasa"],d["tipos_rep_apt"],
     d["oxa72_contexto"] or "-",d["n_contigs"],d["fuente_aislamiento"] or "-"]
    for d in sorted(pe,key=lambda x:(x["anio"],x["accesion"]))])

# --- S4: evaluabilidad ---
n=len(R)
ev=[["MLST Pasteur",sum(1 for d in R if d["st_pasteur"] not in ("-","")),n,"Perfiles alelicos no descritos"],
    ["MLST Oxford",sum(1 for d in R if d["st_oxford"] not in ("-","")),n,"Perfiles no descritos y paralogos del esquema"],
    ["Tipo capsular KL",sum(1 for d in R if d["kl"] not in ("-","")),n,"-"],
    ["Locus OCL",sum(1 for d in R if d["ocl"] not in ("-","")),n,"-"],
    ["Tipificacion plasmidica APT",
     sum(1 for d in R if d["tipos_rep_apt"] not in ("-","","ninguno")),n,"Ausencia de genes rep tipificables"],
    ["Contexto de blaOXA-72",
     sum(1 for d in R if d["oxa72_contexto"] not in ("","no_evaluable")),
     sum(1 for d in R if d["oxa72_contexto"]!=""),"Gen y rep en contigs distintos"],
    ["Contexto de ISAba1 (copias)",44,988,"Fragmentacion del ensamblado"]]
wr("TablaS4_evaluabilidad.tsv",
   ["analisis","evaluables","total","pct","limitacion_dominante"],
   [[a,e,t,round(100*e/t,1),l] for a,e,t,l in ev])

# --- S5: los 58 ensamblados ---
M=tsv(f"{B}/deposito/metadatos_58_ensamblados.tsv"); mh=M[0]
wr("TablaS5_ensamblados_propios.tsv", mh, M[1:])

print("\n=== S4: evaluabilidad ===")
for a,e,t,l in ev: print(f"  {a:32s} {e:4d}/{t:4d}  {100*e/t:5.1f} %")
print(f"\n=== S2: {len(rows)} proyectos ===")
print(f"  paises: {dict(Counter(r[0] for r in rows))}")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Construye la tabla clinica del conjunto: una fila por genoma con los
metadatos curados y el resultado de los procedimientos de tipificacion.
Es la entrada de los pasos 04, 05, 06 y 08.

Entrada:  repo/metadatos/maestra_curada.tsv
          resultados/tipificacion/*.tsv
          resultados/apt_v2/apt_hits_v2.tsv
Salida:   resultados/tabla_clinica_900.tsv
"""
import csv, os, re, collections

B = os.path.expanduser("~/abaumannii")
T = f"{B}/resultados/tipificacion"
SALIDA = f"{B}/resultados/tabla_clinica_900.tsv"

# La familia blaOXA-51-like esta presente en todos los genomas de la especie
# y su expresion basal no confiere resistencia clinica, de modo que se excluye
# del recuento de carbapenemasas adquiridas (apartado 2.4 de metodos).
#
# La exclusion se aplica por la familia que declara el catalogo y no por el
# nombre de la variante: los nombres no permiten reconocer la pertenencia
# (blaOXA-100, blaOXA-407 o blaOXA-241 son todas de esta familia).
INTRINSECA = "OXA-51 family"

def acc_de(ruta):
    m = re.search(r"(GC[AF]_\d+\.\d+)", ruta)
    if m: return m.group(1)
    return os.path.basename(ruta).replace(".fna", "")

def lee_mlst(*ficheros):
    d = {}
    for f in ficheros:
      if not os.path.exists(f): continue
      for l in open(f):
        c = l.rstrip("\n").split("\t")
        if len(c) < 3: continue
        st = c[2].strip()
        d[acc_de(c[0])] = st if st else "-"
    return d

def lee_kaptive(*ficheros):
    d = {}
    for f in ficheros:
      if not os.path.exists(f): continue
      with open(f) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        col = next((c for c in (r.fieldnames or []) if "Best match locus" in c), None)
        if not col: return d
        for row in r:
            v = (row.get(col) or "").strip()
            d[acc_de(row.get("Assembly", ""))] = v if v else "NA"
    return d

meta = {}
with open(f"{B}/repo/metadatos/maestra_curada.tsv", encoding="utf-8") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        a = (row.get("Assembly") or "").strip() or (row.get("Run") or "").strip()
        if a: meta[a] = row

# Los 58 ensamblados propios se tipificaron en una ejecucion separada de
# los publicos, de modo que sus salidas estan en ficheros con sufijo _58.
K = f"{B}/resultados/kaptive"
pas = lee_mlst(f"{T}/mlst_pasteur.tsv", f"{T}/mlst_pasteur_58.tsv")
oxf = lee_mlst(f"{T}/mlst_oxford.tsv",  f"{T}/mlst_oxford_58.tsv")
kl  = lee_kaptive(f"{T}/kaptive_kl.tsv",  f"{K}/kl_886.tsv",  f"{K}/kl_58.tsv")
ocl = lee_kaptive(f"{T}/kaptive_ocl.tsv", f"{K}/ocl_886.tsv", f"{K}/ocl_58.tsv")

carb = collections.defaultdict(set)
for amr_f in (f"{T}/amrfinder_todos.tsv", f"{T}/amrfinder_58.tsv"):
  if not os.path.exists(amr_f): continue
  with open(amr_f) as fh:
    for l in fh:
        c = l.rstrip("\n").split("\t")
        if len(c) < 19: continue
        # Columna 13: subclase. Solo las carbapenemasas.
        if c[12].strip().upper() != "CARBAPENEM": continue
        # Columna 8: familia declarada por el catalogo. La exclusion de la
        # familia intrinseca se aplica por familia y no por nombre de
        # variante, dado que los nombres no permiten reconocerla.
        if INTRINSECA in c[7]: continue
        # Se descarta el nombre generico "blaOXA", que la herramienta emplea
        # cuando el alineamiento es parcial y no permite asignar ni familia
        # ni variante. Los nombres que si identifican familia (blaNDM,
        # blaKPC) o variante (blaOXA-23, blaOXA-253) se conservan.
        gen = c[6].strip()
        if gen == "blaOXA": continue
        carb[c[0].strip()].add(gen)

rep = collections.defaultdict(set)
apt = f"{B}/resultados/apt_v2/apt_hits_v2.tsv"
if os.path.exists(apt):
    for l in open(apt):
        c = l.rstrip("\n").split("\t")
        if len(c) < 2: continue
        t = c[1].strip()
        if re.match(r"^r[13P]-T\d+$", t):
            rep[acc_de(c[0])].add(t)

# El paso 03/01 normaliza el origen a tres categorias en espanol.
CLINICO = {"clinico"}

# El conjunto lo definen los genomas efectivamente tipificados, que son los
# que superaron la confirmacion de especie y el control de contaminacion.
tipificados = set(pas) | set(kl) | set(carb)
filas = []
for a, m in meta.items():
    if (m.get("origen") or "").strip().lower() not in CLINICO: continue
    # Solo los genomas que superaron el control de calidad del paso 03.
    if a not in tipificados: continue
    g = sorted(carb.get(a, []))
    filas.append([a,
        (m.get("pais") or "NA").strip(), pas.get(a, "-"),
        (m.get("anio") or "NA").strip(),
        "si" if g else "no", ",".join(g) if g else "ninguna",
        oxf.get(a, "-"), kl.get(a, "NA"), ocl.get(a, "NA"),
        ",".join(sorted(rep.get(a, []))) if rep.get(a) else "ninguno"])

filas.sort(key=lambda r: r[0])
with open(SALIDA, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accesion","pais","st_pasteur","anio","carbapenemasa",
                "genes_carb","st_oxford","kl","ocl","rep_apt"])
    w.writerows(filas)

print(f"escrito: {SALIDA}")
print(f"  genomas clinicos: {len(filas)}")
if len(filas) != 900: print("  AVISO: se esperaban 900 genomas")

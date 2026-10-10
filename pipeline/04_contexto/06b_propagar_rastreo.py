#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Propaga el resultado del rastreo dirigido de blaOXA-72 a la tabla clinica.

AMRFinderPlus etiqueta como 'blaOXA' generico los alineamientos parciales que
quedan partidos entre fragmentos de un contig. El rastreo dirigido por BLAST
(05_rastreo_dirigido.sh, -perc_identity 95) los recupera. 06_informe_rastreo.py
solo informa; este paso escribe el resultado, de modo que genes_carb coincida
con oxa72_contexto de la Tabla S1.

Entrada:  tmp_rescan/hits.tsv
Salida:   resultados/tabla_clinica_900.tsv     (reescrita en sitio)
          resultados/TablaS1_900_genomas.tsv   (reescrita en sitio)
"""
import re, os, csv, shutil
B = os.path.expanduser("~/abaumannii")
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
TC  = f"{B}/resultados/tabla_clinica_900.tsv"

det = set()
for l in open(f"{B}/tmp_rescan/hits.tsv"):
    if l.strip():
        m = ACC.search(l.split("\t")[0])
        if m: det.add(m.group(1))

shutil.copy(TC, TC + ".bak")
with open(TC, newline="", encoding="utf-8") as fh:
    rows = list(csv.reader(fh, delimiter="\t"))
hdr, datos = rows[0], rows[1:]
iG, iC = hdr.index("genes_carb"), hdr.index("carbapenemasa")

n_gen = n_carb = 0
for r in datos:
    m = ACC.search(r[0])
    if not m or m.group(1) not in det:       continue
    if "blaOXA-72" in r[iG]:                 continue
    r[iG] = "blaOXA-72" if r[iG] in ("", "ninguna") else r[iG] + ",blaOXA-72"
    n_gen += 1
    if r[iC] != "si":
        r[iC] = "si"; n_carb += 1

with open(TC, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(hdr); w.writerows(datos)

# --- misma propagacion sobre la Tabla S1 -------------------------------
S1 = f"{B}/resultados/TablaS1_900_genomas.tsv"
shutil.copy(S1, S1 + ".bak")
with open(S1, newline="", encoding="utf-8") as fh:
    f1 = list(csv.reader(fh, delimiter="\t"))
h1, d1 = f1[0], f1[1:]
jG, jC = h1.index("genes_carbapenemasa"), h1.index("carbapenemasa_adquirida")

m_gen = m_carb = 0
for r in d1:
    m = ACC.search(r[0])
    if not m or m.group(1) not in det:       continue
    if "blaOXA-72" in r[jG]:                 continue
    r[jG] = "blaOXA-72" if r[jG] in ("", "ninguna") else r[jG] + ",blaOXA-72"
    m_gen += 1
    if r[jC] != "si":
        r[jC] = "si"; m_carb += 1

with open(S1, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(h1); w.writerows(d1)

print(f"  detectados por BLAST: {len(det)}")
print(f"  [S1] genes corregidos: {m_gen}   no->si: {m_carb}")
print(f"  genes_carb corregidos:    {n_gen}")
print(f"  carbapenemasa no -> si:   {n_carb}")
print(f"  portadores blaOXA-72 ahora: {sum(1 for r in datos if 'blaOXA-72' in r[iG])}")

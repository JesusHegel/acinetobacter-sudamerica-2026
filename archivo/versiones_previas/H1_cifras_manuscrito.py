#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H1: extrae de los ficheros verificados todas las cifras que aparecen en el
manuscrito. Fuente unica de verdad para la redaccion: ninguna cifra se escribe
a mano, todas se leen de resultados/.
Salida: resultados/CIFRAS_MANUSCRITO.txt
"""
import os, csv, glob, collections
B = os.path.expanduser("~/abaumannii")
L = []
def p(s=""): L.append(s); print(s)

S1 = [{k:(v or '').strip() for k,v in r.items()}
      for r in csv.DictReader(open(f"{B}/resultados/TablaS1_900_genomas.tsv"), delimiter='\t')]
n = len(S1)
p("="*70); p("CIFRAS VERIFICADAS DEL MANUSCRITO"); p("="*70)

p("\n--- 3.1 COMPOSICION ---")
p(f"genomas clinicos: {n}")
for k,v in collections.Counter(d['pais'] for d in S1).most_common(): p(f"  {k}: {v}")
p(f"BioProjects distintos: {len({d['bioproject'] for d in S1})}")
c = collections.Counter((d['pais'],d['bioproject']) for d in S1)
p(f"pares proyecto-pais: {len(c)}")
p(f"estratos con n>=8: {sum(1 for v in c.values() if v>=8)}  (genomas: {sum(v for v in c.values() if v>=8)})")
p(f"BioProjects distintos entre esos estratos: {len({k[1] for k,v in c.items() if v>=8})}")
ev = collections.defaultdict(list)
for d in S1: ev[(d['pais'],d['st_pasteur'],d['bioproject'],d['anio'])].append(d)
p(f"eventos desduplicados: {len(ev)}  (redundancia {100*(1-len(ev)/n):.0f} %)")
p(f"propios: {sum(1 for d in S1 if d['tipo_registro']=='ensamblado_propio')}  "
  f"publicos: {sum(1 for d in S1 if d['tipo_registro']=='ensamblado_publico')}")

p("\n--- 3.2 DETERMINANTES ---")
g = collections.Counter()
for d in S1:
    for t in d['genes_carbapenemasa'].split(','):
        t = t.strip()
        if t and t != 'ninguna': g[t] += 1
for k,v in g.most_common(): p(f"  {k}: {v}")
sc = sum(1 for d in S1 if d['carbapenemasa_adquirida']=='no')
p(f"sin carbapenemasa: {sc} ({100*sc/n:.1f} %)")
st = collections.Counter(d['st_pasteur'] for d in S1)
asig = n - st.get('-',0)
p(f"MLST Pasteur: {asig}/{n} = {100*asig/n:.1f} %   sin asignar: {st.get('-',0)}")
p(f"  top5: {[(k,v) for k,v in st.most_common(6) if k!='-'][:5]}")
top5 = sum(v for k,v in st.most_common(6) if k!='-')
p(f"  suma top5: {top5} = {100*top5/n:.1f} % del total, {100*top5/asig:.1f} % de los asignados")
so = sum(1 for d in S1 if d['st_oxford'] not in ('-',''))
p(f"MLST Oxford: {so}/{n} = {100*so/n:.1f} %")
kl = collections.Counter(d['kl'] for d in S1)
p(f"KL distintos: {len([k for k in kl if k not in ('-','')])}   "
  f"OCL distintos: {len({d['ocl'] for d in S1 if d['ocl'] not in ('-','')})}")
for k,v in kl.most_common(3): p(f"  {k}: {v} ({100*v/n:.1f} %)")
rep = collections.Counter(); con = 0
for d in S1:
    ts = [t.strip() for t in d['tipos_rep_apt'].split(',') if t.strip() not in ('','-','ninguno')]
    if ts: con += 1
    for t in ts: rep[t] += 1
p(f"con rep tipificable: {con}/{n} = {100*con/n:.1f} %")
for k,v in rep.most_common(6): p(f"  {k}: {v}")

p("\n--- 3.5/3.6 blaOXA-72 ---")
port = [d for d in S1 if d['oxa72_contexto'] not in ('','-')]
coloc = [d for d in port if d['oxa72_contexto'] != 'no_evaluable']
p(f"portadores: {len(port)}   co-localizados: {len(coloc)} ({100*len(coloc)/len(port):.1f} %)")
p(f"por pais: {dict(collections.Counter(d['pais'] for d in port))}")
and_ = [d for d in coloc if d['oxa72_contexto']=='r3-T18']
bra = [d for d in coloc if d['pais']=='Brasil']
p(f"andinos r3-T18: {len(and_)}   brasilenos: {len(bra)}")
p(f"  ST de los andinos: {dict(collections.Counter('ST'+d['st_pasteur'] for d in and_))}")
p(f"  contexto brasileno: {dict(collections.Counter(d['oxa72_contexto'] for d in bra))}")
t18 = [d for d in S1 if 'r3-T18' in d['tipos_rep_apt']]
p(f"genomas con r3-T18: {len(t18)}  paises: {len({d['pais'] for d in t18})}  "
  f"linajes: {len({d['st_pasteur'] for d in t18 if d['st_pasteur'] not in ('-','')})}")
p(f"  de esos, con el gen: {sum(1 for d in t18 if d['oxa72_contexto'] not in ('','-'))}")

p("\n--- 3.3 MODELOS ---")
for l in open(f"{B}/resultados/modelo/varianzas.tsv"): p("  " + l.rstrip())

p("\n--- 3.10 EVALUABILIDAD ---")
for l in open(f"{B}/resultados/TablaS4_evaluabilidad.tsv"): p("  " + l.rstrip())
e = [l.rstrip("\n").split("\t")[6] for l in open(f"{B}/resultados/isaba1_oxa51_900.tsv") if l.strip()]
p(f"ISAba1: {collections.Counter(e)}")

p("\n--- VERIFICACIONES ADICIONALES (scripts F) ---")
for f in ("verif/fisher_homogeneo.tsv","verif/cgmlst_estructura.tsv","amr_v4_900.tsv"):
    q = f"{B}/resultados/{f}"
    if os.path.exists(q):
        p(f"\n  [{f}]")
        for i,l in enumerate(open(q)):
            if i < 8: p("    " + l.rstrip())

open(f"{B}/resultados/CIFRAS_MANUSCRITO.txt","w").write("\n".join(L))
p(f"\n\nescrito: resultados/CIFRAS_MANUSCRITO.txt")

import csv, os, collections
B = os.path.expanduser("~/abaumannii")
meta = {r['accesion'].strip(): r for r in
        csv.DictReader(open(f"{B}/resultados/TablaS1_900_genomas.tsv"), delimiter='\t')}
hits = collections.defaultdict(list)
for l in open(f"{B}/resultados/isaba1_oxa23_900.tsv"):
    c = l.rstrip("\n").split("\t")
    if c[6] == "presente":
        hits[c[0]].append(c[2])
filas = []
for g, ctgs in hits.items():
    m = meta.get(g, {})
    filas.append((m.get('pais',''), 'ST'+m.get('st_pasteur','?'), m.get('anio',''),
                  g, m.get('bioproject',''), m.get('kl',''), len(ctgs)))
filas.sort()
out = f"{B}/resultados/anexo_isaba1_oxa23.tsv"
with open(out, "w") as fh:
    fh.write("pais\tlinaje\tanio\taccesion\tbioproject\tkl\tcopias_con_ISAba1\n")
    for f in filas:
        fh.write("\t".join(str(x) for x in f) + "\n")
print(f"{'Pais':<11}{'Linaje':<9}{'Anio':<6}{'Accesion':<19}{'BioProject':<14}{'KL':<8}{'Cop':>4}")
print("-"*74)
for f in filas:
    print(f"{f[0]:<11}{f[1]:<9}{f[2]:<6}{f[3]:<19}{f[4]:<14}{f[5]:<8}{f[6]:>4}")
print(f"\n{len(filas)} genomas · {sum(f[6] for f in filas)} copias · {len({f[4] for f in filas})} BioProjects")
print(f"escrito: {out}")

#!/usr/bin/env bash
# =====================================================================
# apt_colocalizacion.sh
# Tipificacion plasmidica APT con registro del contig, y cruce con
# AMRFinderPlus para determinar co-localizacion gen-replicon.
#
# Sustituye al analisis previo, cuyo intermedio no conservaba el contig
# del genoma y por tanto no era verificable.
#
# Uso:  bash scripts/apt_colocalizacion.sh
# Entorno: conda activate abaumannii   (necesita blastn y makeblastdb)
# =====================================================================
set -euo pipefail

BASE=~/abaumannii
ESQUEMA=$BASE/AcinetobacterPlasmidTyping-main/2026_May_rep_DNA-seqs_V3.fasta
AMR=$BASE/resultados/tipificacion/amrfinder_todos.tsv
AMR58=$BASE/resultados/tipificacion/amrfinder_58.tsv
OUT=$BASE/resultados/apt_v2
TMP=$BASE/tmp_apt

# criterios del esquema APT
MIN_ID=95
MIN_COV=90

mkdir -p "$OUT" "$TMP"

echo "=== PASO 1/6  Localizando genomas ==="
find "$BASE/datos/genomas_945" -name "*_genomic.fna" > "$TMP/rutas.txt" 2>/dev/null || true
find "$BASE/datos/ensamblados_63" -name "*.fna" >> "$TMP/rutas.txt" 2>/dev/null || true
echo "    genomas encontrados: $(wc -l < "$TMP/rutas.txt")"

echo "=== PASO 2/6  Construyendo FASTA combinado con contigs prefijados ==="
# El prefijo garantiza identificadores unicos: los ensamblados propios de
# Shovill usan nombres tipo contig00001, repetidos entre genomas.
: > "$TMP/todos_contigs.fna"
while read -r ruta; do
    # accesion = nombre de carpeta (publicos) o del archivo (propios)
    if [[ "$ruta" == *ensamblados_63* ]]; then
        acc=$(basename "$ruta" .fna)
    else
        acc=$(basename "$(dirname "$ruta")")
    fi
    awk -v A="$acc" '/^>/{sub(/^>/,""); split($0,p," "); print ">"A"|"p[1]; next} {print}' \
        "$ruta" >> "$TMP/todos_contigs.fna"
done < "$TMP/rutas.txt"
echo "    contigs totales: $(grep -c '^>' "$TMP/todos_contigs.fna")"

echo "=== PASO 3/6  Base de datos BLAST ==="
makeblastdb -in "$TMP/todos_contigs.fna" -dbtype nucl -out "$TMP/db_contigs" > /dev/null
echo "    lista"

echo "=== PASO 4/6  BLAST del esquema APT (puede tardar) ==="
blastn -query "$ESQUEMA" -db "$TMP/db_contigs" \
       -outfmt '6 qseqid sseqid pident length qlen slen sstart send evalue bitscore' \
       -perc_identity "$MIN_ID" -max_target_seqs 100000 \
       -num_threads 6 -evalue 1e-20 \
       > "$TMP/blast_crudo.tsv"
echo "    alineamientos: $(wc -l < "$TMP/blast_crudo.tsv")"

echo "=== PASO 5/6  Filtrado y asignacion de tipo ==="
python3 - "$TMP/blast_crudo.tsv" "$OUT/apt_hits_v2.tsv" "$MIN_COV" <<'PYEOF'
import sys, re

entrada, salida, min_cov = sys.argv[1], sys.argv[2], float(sys.argv[3])

# el tipo es todo lo anterior a la accesion de referencia del esquema
PAT = re.compile(r'^(.+?)_([A-Z]{2}_?[A-Z]*\d{6,}\.\d+)_')

def tipo_de(qid):
    m = PAT.match(qid)
    t = m.group(1) if m else qid.split('_')[0]
    # normaliza r3_T3 -> r3-T3 (variantes con guion bajo en el esquema)
    return re.sub(r'^(r[0-9A-Za-z]+)_T', r'\1-T', t)

filas = []
with open(entrada) as f:
    for linea in f:
        c = linea.rstrip('\n').split('\t')
        qseqid, sseqid, pident, length, qlen, slen, sstart, send = c[:8]
        cov = 100.0 * int(length) / int(qlen)
        if cov < min_cov:
            continue
        genoma, contig = sseqid.split('|', 1)
        ini, fin = sorted((int(sstart), int(send)))
        filas.append((genoma, tipo_de(qseqid), contig, ini, fin,
                      int(slen), round(float(pident), 2), round(cov, 1)))

# conserva el mejor alineamiento por (genoma, tipo, contig)
mejor = {}
for f_ in filas:
    k = (f_[0], f_[1], f_[2])
    if k not in mejor or f_[6] > mejor[k][6]:
        mejor[k] = f_

with open(salida, 'w') as out:
    out.write("genoma\ttipo_rep\tcontig\trep_ini\trep_fin\tcontig_len\tpident\tcov\n")
    for f_ in sorted(mejor.values()):
        out.write("\t".join(map(str, f_)) + "\n")

print(f"    alineamientos que pasan el filtro: {len(filas)}")
print(f"    hits unicos (genoma, tipo, contig): {len(mejor)}")
PYEOF

echo "=== PASO 6/6  Cruce con AMRFinderPlus: co-localizacion ==="
python3 - "$OUT/apt_hits_v2.tsv" "$AMR" "$AMR58" "$OUT/colocalizacion.tsv" <<'PYEOF'
import sys
from collections import defaultdict

apt_f, amr_f, amr58_f, salida = sys.argv[1:5]

# columnas de amrfinder (sin encabezado): 1 genoma, 3 contig, 4 ini, 5 fin, 7 gen
GENES = ("blaOXA-72", "blaOXA-23", "blaOXA-58", "blaNDM", "blaOXA-143",
         "blaOXA-24", "blaKPC", "blaIMP", "blaGES", "blaOXA-40")

genes = []
for ruta in (amr_f, amr58_f):
    try:
        with open(ruta) as f:
            for linea in f:
                c = linea.rstrip('\n').split('\t')
                if len(c) < 7:
                    continue
                nombre = c[6]
                if not any(nombre.startswith(g) for g in GENES):
                    continue
                try:
                    ini, fin = sorted((int(c[3]), int(c[4])))
                except ValueError:
                    continue
                genes.append((c[0], c[2], ini, fin, nombre))
    except FileNotFoundError:
        print(f"    AVISO: no se encontro {ruta}")

# reps indexados por (genoma, contig)
reps = defaultdict(list)
contig_len = {}
with open(apt_f) as f:
    next(f)
    for linea in f:
        g, tipo, contig, ini, fin, clen, pid, cov = linea.rstrip('\n').split('\t')
        reps[(g, contig)].append((tipo, int(ini), int(fin)))
        contig_len[(g, contig)] = int(clen)

# reps por genoma, para saber si el genoma tiene algun rep en otro contig
reps_por_genoma = defaultdict(set)
for (g, _c), lst in reps.items():
    for t, _i, _f in lst:
        reps_por_genoma[g].add(t)

with open(salida, 'w') as out:
    out.write("genoma\tgen\tgen_contig\tgen_ini\tgen_fin\tcontig_len\t"
              "estado\treps_mismo_contig\tdist_min_pb\treps_otros_contigs\n")
    for g, contig, ini, fin, nombre in sorted(genes):
        clen = contig_len.get((g, contig), "NA")
        aqui = reps.get((g, contig), [])
        otros = sorted(reps_por_genoma.get(g, set()))
        if aqui:
            estado = "COLOCALIZADO"
            tipos = ",".join(sorted({t for t, _i, _f in aqui}))
            dmin = min(min(abs(ini - rf), abs(ri - fin)) if (rf < ini or ri > fin) else 0
                       for _t, ri, rf in aqui)
        else:
            estado = "NO_EVALUABLE" if otros else "SIN_REP_EN_GENOMA"
            tipos = "-"
            dmin = "NA"
        out.write(f"{g}\t{nombre}\t{contig}\t{ini}\t{fin}\t{clen}\t"
                  f"{estado}\t{tipos}\t{dmin}\t{','.join(otros) if otros else '-'}\n")

print(f"    genes de carbapenemasa procesados: {len(genes)}")
PYEOF

echo
echo "=== RESUMEN ==="
awk -F'\t' 'NR>1{n[$2"\t"$7]++} END{for(k in n) print n[k]"\t"k}' \
    "$OUT/colocalizacion.tsv" | sort -k2,2 -k1,1nr

echo
echo "Salidas:"
echo "  $OUT/apt_hits_v2.tsv      tipificacion con contig, coordenadas y longitud"
echo "  $OUT/colocalizacion.tsv   estado de co-localizacion por gen"
echo
echo "Para liberar espacio:  rm -rf $TMP"

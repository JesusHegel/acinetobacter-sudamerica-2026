#!/usr/bin/env bash
# =====================================================================
# Tipificacion de los 900 genomas: especie, MLST, capsula y resistoma.
#
# Entrada:  datos/genomas_900/        (pasos 01 y 02)
# Produce:  resultados/tipificacion/
# Duracion: 2 a 4 h, dominada por AMRFinderPlus
#
# Emplea dos entornos conda, indicados en cada bloque.
# =====================================================================
set -euo pipefail
BASE=~/abaumannii
IN=$BASE/datos/genomas_900
OUT=$BASE/resultados/tipificacion
REF=$BASE/repo/pipeline/01_datos/referencias_acb
mkdir -p "$OUT"

# ---------------------------------------------------------------------
# 1. Confirmacion de especie            entorno: qc
# ---------------------------------------------------------------------
# ANI frente a las cinco referencias del complejo A. calcoaceticus-baumannii.
# Umbral de 95 % para delimitacion de especie; se asigna la referencia con
# mayor identidad. Las referencias estan en 01_datos/referencias_acb/
echo "=== 1/4  Confirmacion de especie (skani) ==="
ls "$IN"/*.fna > "$OUT/lista_genomas.txt"
ls "$REF"/*.fna > "$OUT/lista_referencias.txt"
skani dist --ql "$OUT/lista_genomas.txt" --rl "$OUT/lista_referencias.txt" \
  -t 4 -o "$OUT/ani_skani.tsv"

# ---------------------------------------------------------------------
# 2. Tipificacion de secuencia multilocus          entorno: abaumannii
# ---------------------------------------------------------------------
# abaumannii_2 = esquema Pasteur (adoptado para el analisis)
# abaumannii   = esquema Oxford  (solo para comparar cobertura)
echo "=== 2/4  MLST ==="
mlst --scheme abaumannii_2 "$IN"/*.fna > "$OUT/mlst_pasteur.tsv"
mlst --scheme abaumannii   "$IN"/*.fna > "$OUT/mlst_oxford.tsv"

# ---------------------------------------------------------------------
# 3. Tipificacion capsular                         entorno: abaumannii
# ---------------------------------------------------------------------
# Las bases vienen incluidas en el propio paquete de Kaptive.
KDB=$(python3 -c "import kaptive,os;print(os.path.join(os.path.dirname(kaptive.__file__),'data'))")
echo "=== 3/4  Kaptive  (bases en $KDB) ==="
kaptive assembly "$KDB/Acinetobacter_baumannii_k_locus_primary_reference.gbk" \
  "$IN"/*.fna -o "$OUT/kaptive_kl.tsv"
kaptive assembly "$KDB/Acinetobacter_baumannii_OC_locus_primary_reference.gbk" \
  "$IN"/*.fna -o "$OUT/kaptive_ocl.tsv"

# ---------------------------------------------------------------------
# 4. Genes de resistencia                          entorno: abaumannii
# ---------------------------------------------------------------------
echo "=== 4/4  AMRFinderPlus ==="
: > "$OUT/amrfinder_todos.tsv"
for f in "$IN"/*.fna; do
  n=$(basename "$f" .fna)
  amrfinder -n "$f" --organism Acinetobacter_baumannii --plus --name "$n" \
    | tail -n +2 >> "$OUT/amrfinder_todos.tsv"
done

echo
echo "Salidas en $OUT"
ls -la "$OUT"

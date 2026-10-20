#!/usr/bin/env bash
# P0.3 — partitioned heritability, real ASD, one focal annotation + full baselineLD v2.2.
#
# Fixes vs first attempt:
#   - --print-coefficients added. The statistic of interest is the annotation's
#     coefficient z-score (tau*), NOT enrichment: enrichment is confounded by
#     annotation size and by overlap with baselineLD categories.
#   - the focal annotation is the LAST row of .results, labelled "L2_1" (baselineLD
#     occupies rows 1..97). The first attempt parsed row 1 and reported the baseline.
set -u
BASE=$(cd "$(dirname "$0")/.." && pwd)
DATA=${DATA:-$BASE/data}
REF=$DATA/ldsc_ref
PY2=${LDSC_PYTHON:-python2}
LDSC=ldsc/ldsc.py
ANN=$BASE/outputs/sldsc/annot
RES=$BASE/outputs/sldsc/results
mkdir -p $RES
SUMSTATS=$BASE/outputs/rg/ASDreal.sumstats.gz
[ -s "$SUMSTATS" ] || { echo "FATAL: $SUMSTATS missing"; exit 1; }

for S in $(ls $ANN/*.1.annot.gz | xargs -n1 basename | sed 's/\.1\.annot\.gz//'); do
  n=$(ls $ANN/${S}.*.l2.ldscore.gz 2>/dev/null | wc -l)
  [ "$n" -ne 22 ] && { echo "[skip] $S ($n/22)"; continue; }
  [ -s $RES/${S}.results ] && [ -s $RES/${S}.ok ] && { echo "[have] $S"; continue; }
  echo "[h2] $S"
  $PY2 $LDSC --h2 $SUMSTATS \
    --ref-ld-chr $REF/baselineLD.,$ANN/${S}. \
    --w-ld-chr $REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --overlap-annot --print-coefficients \
    --frqfile-chr $REF/1000G_Phase3_frq/1000G.EUR.QC. \
    --out $RES/${S} > /dev/null 2>&1
  [ -s $RES/${S}.results ] && touch $RES/${S}.ok || echo "   $S FAILED"
done
echo "SLDSC_H2_DONE"

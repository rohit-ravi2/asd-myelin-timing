#!/usr/bin/env bash
# MAGMA pipeline for ASD GWAS x lipid/myelin gene sets.
# Angle 1 of Test B. Runs both Grove2019 and Matoba2020.
# NOTE ON INDEPENDENCE: Matoba2020 = SPARK + iPSYCH + PGC; Grove2019 = iPSYCH + PGC.
# Matoba is essentially a SUPERSET of Grove, not an independent replication.
# Treat agreement between them as consistency, never as replication.
set -u
BASE=$(cd "$(dirname "$0")/.." && pwd)
DATA=${DATA:-$BASE/data}
COR=$DATA/asd-corpus
MAGMA=$COR/tools/magma/magma
GLOC=$COR/tools/magma_geneloc/NCBI37.3.gene.loc
REF=$COR/tools/magma_ref_eur/g1000_eur
SETS=$BASE/annot/asd_genesets.magma.set
WORK=$BASE/outputs/magma
mkdir -p $WORK
cd $WORK

# ---------------------------------------------------------- 1. prepare inputs
if [ ! -s grove.snploc ]; then
  echo "[prep] Grove2019"
  # REAL Grove2019 (iPSYCH-PGC_ASD_Nov2017): CHR SNP BP A1 A2 INFO OR SE P; N=46350
  zcat $DATA/asd-corpus/genetics/ASD_Grove2019_REAL.gz | \
    awk 'NR>1 && $2 ~ /^rs/ {print $2"\t"$1"\t"$3}' > grove.snploc
  zcat $DATA/asd-corpus/genetics/ASD_Grove2019_REAL.gz | \
    awk 'BEGIN{print "SNP\tP\tN"} NR>1 && $2 ~ /^rs/ && $9+0>0 {print $2"\t"$9"\t46350"}' > grove.pval
  wc -l grove.snploc grove.pval
fi

if [ ! -s matoba.snploc ]; then
  echo "[prep] Matoba2020"
  awk -F'\t' 'NR>1 {
      n=split($3,a,"_"); rs=a[n];
      if (rs ~ /^rs/) print rs"\t"$1"\t"$2
    }' $DATA/dvad-corpus/asd_scope/gwas/ASD_Matoba2020_SPARK_iPSYCH_PGC.tsv > matoba.snploc
  awk -F'\t' 'BEGIN{print "SNP\tP\tN"} NR>1 {
      n=split($3,a,"_"); rs=a[n];
      if (rs ~ /^rs/ && $8+0>0) print rs"\t"$8"\t"$10
    }' $DATA/dvad-corpus/asd_scope/gwas/ASD_Matoba2020_SPARK_iPSYCH_PGC.tsv > matoba.pval
  wc -l matoba.snploc matoba.pval
fi

# ---------------------------------------------------------- 2. annotate + gene + set
for G in grove matoba; do
  echo "=============== $G ==============="

  # window sensitivity: 0kb (gene body only) and 35/10kb (standard) -- ANGLE
  for W in "0,0:w0" "35,10:w35"; do
    WIN=${W%%:*}; TAG=${W##*:}
    if [ ! -s ${G}.${TAG}.genes.annot ]; then
      $MAGMA --annotate window=$WIN --snp-loc ${G}.snploc --gene-loc $GLOC \
             --out ${G}.${TAG} > /dev/null 2>&1
    fi
    if [ ! -s ${G}.${TAG}.genes.raw ]; then
      echo "[gene] $G $TAG"
      $MAGMA --bfile $REF --pval ${G}.pval ncol=N \
             --gene-annot ${G}.${TAG}.genes.annot --out ${G}.${TAG} 2>&1 | tail -3
    fi
    # competitive gene-set analysis (default; corrects for gene size/density)
    if [ ! -s ${G}.${TAG}.gsa.out ]; then
      echo "[set-competitive] $G $TAG"
      $MAGMA --gene-results ${G}.${TAG}.genes.raw --set-annot $SETS \
             --out ${G}.${TAG} 2>&1 | tail -2
    fi
    # self-contained as a second angle
    if [ ! -s ${G}.${TAG}.selfcont.gsa.out ]; then
      echo "[set-selfcontained] $G $TAG"
      $MAGMA --gene-results ${G}.${TAG}.genes.raw --set-annot $SETS \
             --model self-contained --out ${G}.${TAG}.selfcont 2>&1 | tail -2
    fi
  done
done
echo "MAGMA_DONE"

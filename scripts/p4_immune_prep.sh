#!/usr/bin/env bash
# MAGMA gene analysis for the ALLERGY (immune-liability) GWAS, so the identical
# gate-competition and OL-pseudotime machinery can be run on an immune trait.
# If immune liability localises to OL where ASD does not, the immune->OL genetic
# bridge is real and ASD simply enters elsewhere. If neither localises, the
# genetic form of the immune->OL bridge is absent in both.
set -u
BASE=$(cd "$(dirname "$0")/.." && pwd); DATA=${DATA:-$BASE/data}
COR=$DATA/asd-corpus; W=$BASE/outputs/magma
cd $W
if [ ! -s allergy.pval ]; then
  zcat $DATA/dvad-corpus/asd_myelin/ASD_Grove2019.gz | \
    awk 'NR>1 && $10 ~ /^rs/ {print $10"\t"$2"\t"$3}' > allergy.snploc
  zcat $DATA/dvad-corpus/asd_myelin/ASD_Grove2019.gz | \
    awk 'BEGIN{print "SNP\tP\tN"} NR>1 && $10 ~ /^rs/ && $8+0>0 {print $10"\t"$8"\t"$13}' > allergy.pval
fi
$COR/tools/magma/magma --annotate window=35,10 --snp-loc allergy.snploc \
   --gene-loc $COR/tools/magma_geneloc/NCBI37.3.gene.loc --out allergy.w35 >/dev/null 2>&1
$COR/tools/magma/magma --bfile $COR/tools/magma_ref_eur/g1000_eur --pval allergy.pval ncol=N \
   --gene-annot allergy.w35.genes.annot --out allergy.w35 2>&1 | tail -2
echo ALLERGY_MAGMA_DONE

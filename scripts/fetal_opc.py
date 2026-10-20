#!/usr/bin/env python3
"""Is there a PRENATAL window for the convergence genes in the OL lineage?
GSE155488 = human FETAL oligodendrocyte-lineage cells. If beta-catenin/TCF and
the OL-convergence genes are already expressed in fetal OPCs, a prenatal route
exists. If they only come on postnatally, the timing objection stands."""
import gzip
import numpy as np, pandas as pd
from pathlib import Path
G=Path("data/asd-corpus/celltype/GSE155488")
C=Path("data/asd-corpus"); P=Path(".")
df=pd.read_csv(G/"GSE155488_All_fetal_OL_lineage_subset.csv.gz",index_col=0)
print(f"fetal OL-lineage matrix: {df.shape[0]} genes x {df.shape[1]} cells")
expr=(df>0).mean(axis=1)          # fraction of fetal OL-lineage cells expressing
mean_e=df.mean(axis=1)
GENES={
 "beta-catenin/TCF (the robust finding)": ["CTNNB1","LDB1","TCF4","TCF7L2","TLE3"],
 "ENSHEATHMENT convergence (SFARI_SCORE1)": ["ANK2","CTNNB1","HRAS","MTOR","NF1","POMGNT1",
                                             "PTEN","SCN2A","SCN8A","SKI","TSC1"],
 "OPC identity markers": ["PDGFRA","CSPG4","SOX10","OLIG1","OLIG2","NKX2-2"],
 "mature myelin (should be LOW prenatally)": ["MBP","PLP1","MOG","MAG","MOBP","CNP","MYRF"],
 "negative controls (non-brain)": ["ALB","MYH7","KRT14","INS","CELA3A"],
}
print(f"\n{'gene':<12}{'% fetal OL cells expressing':>30}{'mean expr':>12}")
for label,gl in GENES.items():
    print(f"\n-- {label}")
    for g in gl:
        if g in expr.index:
            print(f"  {g:<12}{100*expr[g]:>28.1f}%{mean_e[g]:>12.3f}")
        else:
            print(f"  {g:<12}{'not in matrix':>29}")
    present=[g for g in gl if g in expr.index]
    if present:
        print(f"  {'SET MEAN':<12}{100*expr[present].mean():>28.1f}%")
# where do the convergence genes sit in the overall distribution?
conv=[g for g in GENES["ENSHEATHMENT convergence (SFARI_SCORE1)"] if g in expr.index]
bcat=[g for g in GENES["beta-catenin/TCF (the robust finding)"] if g in expr.index]
allq=expr.rank(pct=True)
print(f"\n=== expression percentile within the fetal OL transcriptome ===")
print(f"  beta-catenin/TCF genes : median {100*allq[bcat].median():.0f}th pct")
print(f"  convergence genes      : median {100*allq[conv].median():.0f}th pct")
print(f"  OPC identity markers   : median {100*allq[[g for g in GENES['OPC identity markers'] if g in expr.index]].median():.0f}th pct")
print(f"  mature myelin genes    : median {100*allq[[g for g in GENES['mature myelin (should be LOW prenatally)'] if g in expr.index]].median():.0f}th pct")
print(f"  genome median          : 50th pct by construction")

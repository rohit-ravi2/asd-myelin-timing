#!/usr/bin/env python3
"""Characterise the top ASD genes: protein class, adult cell-type expression,
oligodendrocyte-lineage enrichment, fetal OPC expression, framework assignment.
Then find the genes that are OL-enriched but NOT yet in the framework."""
import csv, gzip, glob, os
import numpy as np, pandas as pd
from pathlib import Path
C=Path("data/asd-corpus"); S=Path("data/dvad-corpus/asd_scope/rare_variant")
V=Path("data/dvad-corpus/asd_scope/velmeshev2019"); P=Path(".")
z=np.load(P/"outputs/testB/pseudobulk.npz",allow_pickle=True)
PB,ncell,keys,genes=z["PB"],z["ncell"],z["keys"],z["genes"]
info=pd.read_csv(P/"outputs/testB/individuals.csv"); info["individual"]=info.individual.astype(str)
dx=info.set_index("individual").diagnosis.to_dict()
ind=np.array([k.split("||")[0] for k in keys]); ct=np.array([k.split("||")[1] for k in keys])
ok=(ncell>=50)&np.array([dx.get(i)=="Control" for i in ind])     # CONTROLS only = normal expression
cpm=np.zeros_like(PB); s=PB.sum(axis=0); nz=s>0
cpm[:,nz]=PB[:,nz]/s[nz]*1e6
lg=np.log2(cpm+1)
CT=sorted(set(ct[ok]))
M=pd.DataFrame({c: lg[:,ok&(ct==c)].mean(axis=1) for c in CT}, index=genes)
M=M.groupby(level=0).max()
OLL=["OPC","Oligodendrocytes"]
M["OL_lineage"]=M[OLL].mean(axis=1)
M["all_mean"]=M[CT].mean(axis=1)
M["OL_enrich"]=M.OL_lineage-M.all_mean          # log2 difference
M["PV"]=M["IN-PV"] if "IN-PV" in M.columns else np.nan
M["PV_enrich"]=M.PV-M.all_mean
# fetal OPC
fet=pd.read_csv(C/"celltype/GSE155488/GSE155488_All_fetal_OL_lineage_subset.csv.gz",index_col=0)
fetfrac=(fet>0).mean(axis=1)
# protein class from Velmeshev marker annotation
pclass={}
for f in glob.glob(str(V/"markers"/"*.tsv.gz")):
    with gzip.open(f,"rt") as fh:
        hdr=fh.readline().rstrip("\n").split("\t")
        try: si,ci=hdr.index("symbol"),hdr.index("_hprdClass")
        except ValueError: continue
        for line in fh:
            p=line.rstrip("\n").split("\t")
            if len(p)>max(si,ci) and p[ci].strip() and p[si] not in pclass:
                pclass[p[si]]=p[ci].strip()[:34]
gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
g=lambda t:gmt.get(t,set())
FW={"OL":g("GOBP_ENSHEATHMENT_OF_NEURONS")|g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOCC_MYELIN_SHEATH"),
    "BCAT":g("GOCC_BETA_CATENIN_TCF_COMPLEX"),
    "mTOR/RAS":{"MTOR","PTEN","TSC1","TSC2","NF1","HRAS","AKT3","PIK3CA","BRAF","SYNGAP1"},
    "SYN":g("GOCC_SYNAPSE"),"CHROM":g("GOBP_CHROMATIN_REMODELING")|g("GOCC_HISTONE_DEACETYLASE_COMPLEX")}
sat={r["hugoGene"] for r in csv.DictReader(open(S/"Satterstrom2020_102genes.csv")) if r.get("hugoGene")}
fu={r["HGNC"] for r in csv.DictReader(open(S/"Fu2022_72_studywide_genes.csv")) if r.get("HGNC")}
GOLD=sorted(sat|fu)
rows=[]
for gene in GOLD:
    if gene not in M.index: continue
    fw=[k for k,v in FW.items() if gene in v]
    rows.append({"gene":gene,"protein_class":pclass.get(gene,""),
        "OL_enrich":round(float(M.OL_enrich[gene]),2),
        "PV_enrich":round(float(M.PV_enrich[gene]),2),
        "OPC":round(float(M.OPC[gene]),2),"Oligo":round(float(M.Oligodendrocytes[gene]),2),
        "all_mean":round(float(M.all_mean[gene]),2),
        "fetalOPC_pct":round(100*float(fetfrac[gene]),1) if gene in fetfrac.index else np.nan,
        "framework":"+".join(fw) if fw else "—"})
R=pd.DataFrame(rows).sort_values("OL_enrich",ascending=False)
print(f"gold-standard genes with expression data: {len(R)}\n")
print("=== TOP 22 by oligodendrocyte-lineage enrichment (adult human cortex, controls) ===")
print(R.head(22).to_string(index=False))
print("\n=== genes OL-enriched (>0.3 log2) but NOT in the OL/BCAT/mTOR framework ===")
cand=R[(R.OL_enrich>0.3)&(~R.framework.str.contains("OL|BCAT|mTOR"))]
print(cand.to_string(index=False))
print(f"\n=== framework genes' OL enrichment, for comparison ===")
print(R[R.framework.str.contains("OL|BCAT|mTOR")].to_string(index=False))
R.to_csv(P/"outputs/gene_characterisation.csv",index=False)
print(f"\nmedian OL_enrich across all gold genes: {R.OL_enrich.median():+.2f}")
print("wrote outputs/gene_characterisation.csv")

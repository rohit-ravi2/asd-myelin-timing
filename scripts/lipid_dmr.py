#!/usr/bin/env python3
"""Are LIPID TRANSPORT genes differentially methylated in ASD, and in which
direction? Uses GSE298243 PFC single-nucleus methylation: 2.5M 500bp bins with
cell type, cytosine context, gene symbol, diagnosis significance, and predicted
ASD vs control methylation (so direction is recoverable).

Background = bins of the SAME cell type and cytosine context (as the original
pipeline specifies), so enrichment is not confounded by coverage differences."""
import gzip, json
import numpy as np, pandas as pd
from scipy import stats
from pathlib import Path
C=Path("data/asd-corpus"); D=Path("data/dvad-corpus/asd_scope/celltype/GSE298243_Age_DMRs_PrespecifiedRegions.tsv.gz")
P=Path(".")
cols=["CellType","CytoContext","Symbol","ReportSignif_Diagnosis",
      "Pred_CTL_Young","Pred_ASD_Young","Pred_CTL_Mid","Pred_ASD_Mid","Pred_CTL_Old","Pred_ASD_Old"]
print("loading DMR table...",flush=True)
df=pd.read_csv(D,sep="\t",usecols=cols,low_memory=False)
df["Symbol"]=df.Symbol.astype(str).str.strip()
df=df[df.Symbol.ne("")&df.Symbol.ne("nan")]
df["sig"]=df.ReportSignif_Diagnosis.astype(str).str.strip().eq("*")
for a in ["Young","Mid","Old"]:
    df[f"d{a}"]=pd.to_numeric(df[f"Pred_ASD_{a}"],errors="coerce")-pd.to_numeric(df[f"Pred_CTL_{a}"],errors="coerce")
df["delta"]=df[["dYoung","dMid","dOld"]].mean(axis=1)
print(f"bins with a gene symbol: {len(df):,}; diagnosis-significant: {df.sig.sum():,}")
gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
g=lambda t:gmt.get(t,set())
SETS={
 "LIPID_TRANSPORT": g("GOBP_LIPID_TRANSPORT")|g("GOBP_LIPID_LOCALIZATION")|g("GOBP_INTRACELLULAR_LIPID_TRANSPORT"),
 "CHOL_STEROL_TRANSPORT": g("GOBP_CHOLESTEROL_TRANSPORT")|g("GOBP_STEROL_TRANSPORT"),
 "LIPOPROTEIN_PARTICLE": g("GOCC_HIGH_DENSITY_LIPOPROTEIN_PARTICLE")|g("GOCC_LOW_DENSITY_LIPOPROTEIN_PARTICLE")|g("GOCC_PLASMA_LIPOPROTEIN_PARTICLE"),
 "FATTY_ACID_TRANSPORT": g("GOBP_FATTY_ACID_TRANSPORT")|g("GOBP_LONG_CHAIN_FATTY_ACID_TRANSPORT"),
 "LIPID_BIOSYNTHESIS": g("GOBP_LIPID_BIOSYNTHETIC_PROCESS"),
 "ABC_SLC_lipid_curated": {"ABCA1","ABCA2","ABCA7","ABCG1","ABCG4","ABCD1","ABCD2","ABCD3","APOE","APOA1",
     "APOB","APOD","APOJ","CLU","LDLR","LRP1","LRP2","LRP8","VLDLR","SCARB1","NPC1","NPC2","STARD3",
     "OSBPL1A","PLTP","CETP","SOAT1","MVK","FDFT1","HMGCR","SREBF1","SREBF2","SCAP","INSIG1"},
 "MYELIN_OL_reference": g("GOBP_ENSHEATHMENT_OF_NEURONS")|g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOCC_MYELIN_SHEATH"),
 "SYNAPSE_reference": g("GOCC_SYNAPSE"),
}
rows=[]
for (ctype,ctx),sub in df.groupby(["CellType","CytoContext"]):
    if len(sub)<2000: continue
    for sn,ss in SETS.items():
        inset=sub.Symbol.isin(ss)
        n_in=int(inset.sum())
        if n_in<50: continue
        sig_in=int(sub.sig[inset].sum()); sig_out=int(sub.sig[~inset].sum())
        n_out=int((~inset).sum())
        if sig_in+sig_out==0: continue
        tab=[[sig_in,n_in-sig_in],[sig_out,n_out-sig_out]]
        orr,p=stats.fisher_exact(tab,alternative="two-sided")
        d_sig=sub.delta[inset&sub.sig].dropna()
        rows.append({"cell":ctype,"ctx":ctx,"set":sn,"n_bins":n_in,"sig":sig_in,
                     "rate_in":sig_in/n_in,"rate_out":sig_out/max(n_out,1),"OR":orr,"p":p,
                     "mean_delta_sig":float(d_sig.mean()) if len(d_sig) else np.nan,
                     "n_hyper":int((d_sig>0).sum()),"n_hypo":int((d_sig<0).sum())})
R=pd.DataFrame(rows)
R=R.sort_values("p")
m=len(R); R["q"]=(R.p*m/np.arange(1,m+1)).cummax().clip(upper=1)
print(f"\ntests: {m} (cell type x cytosine context x gene set)")
print(f"\n--- LIPID sets, OPC and Oligo (the cells of interest) ---")
sel=R[(R.cell.isin(["OPC","Oligo"]))&(R["set"].str.contains("LIPID|CHOL|LIPOPROTEIN|FATTY|ABC"))]
print(sel[["cell","ctx","set","n_bins","sig","rate_in","rate_out","OR","p","q","mean_delta_sig","n_hyper","n_hypo"]].to_string(index=False,float_format=lambda x:f"{x:.4g}"))
print(f"\n--- reference sets, OPC and Oligo ---")
sel2=R[(R.cell.isin(["OPC","Oligo"]))&(R["set"].str.contains("MYELIN|SYNAPSE"))]
print(sel2[["cell","ctx","set","n_bins","sig","rate_in","rate_out","OR","p","q","mean_delta_sig","n_hyper","n_hypo"]].to_string(index=False,float_format=lambda x:f"{x:.4g}"))
print(f"\n--- anything significant at q<0.1 across all cells/sets ---")
s=R[R.q<0.1]
print(s[["cell","ctx","set","n_bins","sig","OR","p","q","mean_delta_sig"]].to_string(index=False,float_format=lambda x:f"{x:.4g}") if len(s) else "  NONE")
R.to_csv(P/"outputs/phase7/lipid_dmr_results.csv",index=False)
print("\nwrote outputs/phase7/lipid_dmr_results.csv")

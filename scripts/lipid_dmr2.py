#!/usr/bin/env python3
"""Gene/promoter-level DMRs are too few (34 total). So: map the 5,302
diagnosis-significant 500bp bins to genes by coordinate, then test lipid
transport vs myelin vs synapse, with direction."""
import gzip
import numpy as np, pandas as pd
from scipy import stats
from pathlib import Path
C=Path("data/asd-corpus"); D=Path("data/dvad-corpus/asd_scope/celltype/GSE298243_Age_DMRs_PrespecifiedRegions.tsv.gz")
P=Path(".")
use=["FeatureType","CellType","CytoContext","Chr","Start","End","Symbol",
     "ReportSignif_Diagnosis","Pred_CTL_Young","Pred_ASD_Young","Pred_CTL_Mid",
     "Pred_ASD_Mid","Pred_CTL_Old","Pred_ASD_Old"]
df=pd.read_csv(D,sep="\t",usecols=use,low_memory=False)
df["sig"]=df.ReportSignif_Diagnosis.astype(str).str.strip().eq("*")
for a in ["Young","Mid","Old"]:
    df[f"d{a}"]=pd.to_numeric(df[f"Pred_ASD_{a}"],errors="coerce")-pd.to_numeric(df[f"Pred_CTL_{a}"],errors="coerce")
df["delta"]=df[["dYoung","dMid","dOld"]].mean(axis=1)
gp=df[df.FeatureType.isin(["gene","promoter"])].copy()
gp["Symbol"]=gp.Symbol.astype(str).str.strip()
print("=== the 34 gene/promoter-level diagnosis-DMRs ===")
h=gp[gp.sig]
print(h[["FeatureType","CellType","CytoContext","Symbol","delta"]].to_string(index=False,float_format=lambda x:f"{x:+.4f}"))
# ---- map significant bins to genes by coordinate
bins=df[df.FeatureType.eq("bin500bp")].copy()
print(f"\nbins: {len(bins):,}; diagnosis-significant: {int(bins.sig.sum()):,}")
genes={}
with gzip.open("data/refs/ensembl115.chr.gff3.gz","rt") as f:
    for line in f:
        if line.startswith("#"): continue
        p=line.split("\t")
        if len(p)<9 or p[2]!="gene": continue
        attr="\t".join(p[8:])
        if "Name=" not in attr: continue
        nm=attr.split("Name=")[1].split(";")[0].strip()
        ch=p[0] if p[0].startswith("chr") else "chr"+p[0]
        genes.setdefault(ch,[]).append((int(p[3]),int(p[4]),nm))
print(f"gene models loaded for {len(genes)} contigs")
def assign(ch,s,e):
    L=genes.get(ch)
    if not L: return None
    for gs,ge,nm in L:
        if s<ge and e>gs: return nm
    return None
idx={}
for ch,L in genes.items():
    L.sort(); idx[ch]=(np.array([x[0] for x in L]),np.array([x[1] for x in L]),[x[2] for x in L])
def assign_fast(ch,s,e):
    if ch not in idx: return None
    st,en,nm=idx[ch]
    j=np.searchsorted(st,e)
    for k in range(max(0,j-40),j):
        if en[k]>s: return nm[k]
    return None
bins["gene"]=[assign_fast(c,s,e) for c,s,e in zip(bins.Chr,bins.Start,bins.End)]
mapped=bins.gene.notna()
print(f"bins mapped to a gene: {mapped.sum():,} ({100*mapped.mean():.1f}%); "
      f"significant+mapped: {int((bins.sig&mapped).sum()):,}")
bins=bins[mapped]
gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
g=lambda t:gmt.get(t,set())
SETS={
 "LIPID_TRANSPORT": g("GOBP_LIPID_TRANSPORT")|g("GOBP_LIPID_LOCALIZATION"),
 "CHOL_STEROL_TRANSPORT": g("GOBP_CHOLESTEROL_TRANSPORT")|g("GOBP_STEROL_TRANSPORT"),
 "LIPID_BIOSYNTHESIS": g("GOBP_LIPID_BIOSYNTHETIC_PROCESS"),
 "LIPID_curated": {"ABCA1","ABCA2","ABCA7","ABCG1","ABCG4","ABCD1","ABCD2","ABCD3","APOE","APOD",
   "CLU","LDLR","LRP1","LRP2","LRP8","VLDLR","SCARB1","NPC1","NPC2","PLTP","SOAT1","HMGCR",
   "SREBF1","SREBF2","SCAP","INSIG1","FDFT1","DHCR7","DHCR24","CYP46A1","SLC27A1","FABP7"},
 "MYELIN_OL": g("GOBP_ENSHEATHMENT_OF_NEURONS")|g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOCC_MYELIN_SHEATH"),
 "SYNAPSE": g("GOCC_SYNAPSE"),
}
rows=[]
for (ct,cx),sub in bins.groupby(["CellType","CytoContext"]):
    if sub.sig.sum()<10: continue
    for sn,ss in SETS.items():
        m=sub.gene.isin(ss); n_in=int(m.sum())
        if n_in<200: continue
        si=int(sub.sig[m].sum()); so=int(sub.sig[~m].sum()); no=int((~m).sum())
        orr,p=stats.fisher_exact([[si,n_in-si],[so,no-so]])
        d=sub.delta[m&sub.sig].dropna()
        rows.append({"cell":ct,"ctx":cx,"set":sn,"bins_in":n_in,"sig_in":si,
            "rate_in":si/n_in,"rate_out":so/no,"OR":orr,"p":p,
            "mean_delta":float(d.mean()) if len(d) else np.nan,
            "hyper":int((d>0).sum()),"hypo":int((d<0).sum())})
R=pd.DataFrame(rows).sort_values("p")
m=len(R); R["q"]=(R.p*m/np.arange(1,m+1)).cummax().clip(upper=1)
print(f"\n=== OPC and Oligo ({m} tests total) ===")
print(R[R.cell.isin(["OPC","Oligo"])][["cell","ctx","set","bins_in","sig_in","rate_in","rate_out","OR","p","q","mean_delta","hyper","hypo"]].to_string(index=False,float_format=lambda x:f"{x:.4g}"))
print(f"\n=== anything q<0.1 ===")
s=R[R.q<0.1]
print(s[["cell","ctx","set","sig_in","OR","p","q","mean_delta","hyper","hypo"]].to_string(index=False,float_format=lambda x:f"{x:.4g}") if len(s) else "  NONE")
R.to_csv(P/"outputs/phase7/lipid_dmr_binlevel.csv",index=False)
print("\nwrote outputs/phase7/lipid_dmr_binlevel.csv")

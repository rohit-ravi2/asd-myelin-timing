#!/usr/bin/env python3
"""Is the global cross-cell-type co-dysregulation (median r=+0.63) biological,
or residual technical variation (RIN, PMI, age, batch)?"""
import numpy as np, pandas as pd
from scipy import stats
from pathlib import Path
V=Path("data/dvad-corpus/asd_scope/velmeshev2019"); P=Path(".")
D=pd.read_csv(P/"outputs/testB/dysregulation_scores.csv")
meta=pd.read_csv(V/"meta.tsv",sep="\t")
cov=(meta.groupby("individual")
        .agg(RIN=("RNA Integrity Number","first"),PMI=("post-mortem interval (hours)","first"),
             age=("age","first"),sex=("sex","first"),cap=("Capbatch","first"),seq=("Seqbatch","first"),
             mito=("RNA mitochondr. percent","mean"),ngene=("genes","mean"))
        .reset_index())
cov["individual"]=cov.individual.astype(str)
D["individual"]=D.individual.astype(str)
D=D.merge(cov,on="individual",how="left")
A=D[D.diagnosis=="ASD"]
print("Per-individual technical covariates vs mean dysregulation (ASD only):")
ind=A.groupby("individual").agg(dys=("dysreg_adj","mean"),RIN=("RIN","first"),PMI=("PMI","first"),
                                age=("age","first"),mito=("mito","first"),ngene=("ngene","first")).dropna()
for c in ["RIN","PMI","age","mito","ngene"]:
    r,p=stats.spearmanr(ind[c],ind.dys)
    print(f"   {c:<7} rho={r:+.3f} p={p:.4f} {'<-- confound' if p<0.05 else ''}")
# residualise dysregulation on technical covariates per cell type, then redo pair correlations
rows=[]
for c,g in D.groupby("celltype"):
    g=g.dropna(subset=["RIN","PMI","mito","ngene"])
    if len(g)<12: continue
    X=np.column_stack([np.ones(len(g)),g.RIN.values,g.PMI.values,g.mito.values,
                       np.log10(g.ngene.values),g.age.values])
    y=g.dysreg_adj.values
    r=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
    gg=g.copy(); gg["dys_tech_adj"]=r; rows.append(gg)
D2=pd.concat(rows)
A2=D2[D2.diagnosis=="ASD"].pivot(index="individual",columns="celltype",values="dys_tech_adj")
cols=list(A2.columns); pr=[]
for i in range(len(cols)):
    for j in range(i+1,len(cols)):
        s=A2[[cols[i],cols[j]]].dropna()
        if len(s)<8: continue
        rr,pp=stats.pearsonr(s[cols[i]],s[cols[j]]); pr.append({"a":cols[i],"b":cols[j],"r":rr,"p":pp,"n":len(s)})
PD=pd.DataFrame(pr).sort_values("r",ascending=False).reset_index(drop=True); PD["rank"]=PD.index+1
print(f"\nAFTER technical adjustment (RIN,PMI,mito,ngene,age):")
print(f"  {len(PD)} pairs; median r = {PD.r.median():+.3f} (was +0.631)")
def show(a,b):
    m=PD[((PD.a==a)&(PD.b==b))|((PD.a==b)&(PD.b==a))]
    if not len(m): return f"  {a} x {b}: n/a"
    r=m.iloc[0]; pct=100*(1-(r['rank']-1)/len(PD))
    return f"  {a:<9} x {b:<18} r={r['r']:+.3f} p={r['p']:.4f} rank {int(r['rank'])}/{len(PD)} ({pct:.0f}th pct)"
for b in ["OPC","Oligodendrocytes"]: print(show("IN-PV",b))
for a in ["IN-SST","IN-VIP"]:
    for b in ["OPC"]: print(show(a,b))
print(show("L2/3","OPC"))
PD.to_csv(P/"outputs/testB/pair_correlations_techadj.csv",index=False)

#!/usr/bin/env python3
"""TEST B — do PV-interneuron and OPC/OL dysregulation co-occur in the SAME
ASD brains?

Pre-registered interpretation (committed before looking):
  H1 supported  = PV-OPC (and/or PV-OL) dysregulation correlation is positive AND
                  ranks high among ALL cell-type pairs (specific coupling).
  NULL          = correlation absent, or present but typical of all pairs
                  (a global 'this brain is disturbed' factor, which is trivial).

Confound control: per-individual cell yield drives pseudobulk noise, which would
make every cell type's 'dysregulation' correlate. Dysregulation scores are
residualised on log(n_cells) before any correlation is computed.
"""
import json
import numpy as np, pandas as pd
from scipy import stats
from pathlib import Path
P=Path("."); OUT=P/"outputs/testB"
z=np.load(OUT/"pseudobulk.npz",allow_pickle=True)
PB,ncell,keys,genes=z["PB"],z["ncell"],z["keys"],z["genes"]
info=pd.read_csv(OUT/"individuals.csv")
info["individual"]=info.individual.astype(str)
info=info.set_index("individual")
ind=np.array([k.split("||")[0] for k in keys]); ct=np.array([k.split("||")[1] for k in keys])
MINCELL=50
ok=ncell>=MINCELL
print(f"groups: {len(keys)} total, {ok.sum()} with >={MINCELL} cells")
cpm=np.zeros_like(PB)
s=PB.sum(axis=0); nz=s>0
cpm[:,nz]=PB[:,nz]/s[nz]*1e6
lg=np.log2(cpm+1)
dx=info.diagnosis.to_dict()
CTS=sorted(set(ct[ok]))
rows=[]
for c in CTS:
    sel=np.where(ok&(ct==c))[0]
    inds=ind[sel]
    diag=np.array([dx.get(i,"NA") for i in inds])
    if (diag=="Control").sum()<5 or (diag=="ASD").sum()<5: continue
    X=lg[:,sel]
    expressed=(X>0).mean(axis=1)>0.5
    X=X[expressed]
    ctrl=X[:,diag=="Control"]
    mu=ctrl.mean(axis=1); sd=ctrl.std(axis=1,ddof=1); sd[sd==0]=np.nan
    Z=(X-mu[:,None])/sd[:,None]
    score=np.nanmean(np.abs(Z),axis=0)
    for i,g in enumerate(sel):
        rows.append({"individual":inds[i],"celltype":c,"diagnosis":diag[i],
                     "dysreg":float(score[i]),"ncell":int(ncell[g]),
                     "ngenes":int(expressed.sum())})
D=pd.DataFrame(rows)
print(f"cell types analysed: {D.celltype.nunique()}; rows {len(D)}")
# residualise on log(ncell) within cell type
D["logn"]=np.log10(D.ncell)
res=[]
for c,g in D.groupby("celltype"):
    X=np.column_stack([np.ones(len(g)),g.logn.values])
    r=g.dysreg.values-X@np.linalg.lstsq(X,g.dysreg.values,rcond=None)[0]
    gg=g.copy(); gg["dysreg_adj"]=r; res.append(gg)
D=pd.concat(res)
print(f"corr(dysreg, log n_cells) before adj = "
      f"{stats.spearmanr(D.dysreg,D.logn)[0]:+.3f}; after = {stats.spearmanr(D.dysreg_adj,D.logn)[0]:+.3f}")
A=D[D.diagnosis=="ASD"].pivot(index="individual",columns="celltype",values="dysreg_adj")
print(f"\nASD individuals x cell types matrix: {A.shape}")
pairs=[]
cols=list(A.columns)
for i in range(len(cols)):
    for j in range(i+1,len(cols)):
        sub=A[[cols[i],cols[j]]].dropna()
        if len(sub)<8: continue
        r,p=stats.pearsonr(sub[cols[i]],sub[cols[j]])
        pairs.append({"a":cols[i],"b":cols[j],"n":len(sub),"r":r,"p":p})
PD=pd.DataFrame(pairs).sort_values("r",ascending=False).reset_index(drop=True)
PD["rank"]=PD.index+1
print(f"\nall cell-type pairs: {len(PD)}; median r = {PD.r.median():+.3f}; "
      f"IQR {PD.r.quantile(.25):+.3f}..{PD.r.quantile(.75):+.3f}")
def show(a,b):
    m=PD[((PD.a==a)&(PD.b==b))|((PD.a==b)&(PD.b==a))]
    if len(m)==0: return f"  {a} x {b}: not testable"
    r=m.iloc[0]
    pct=100*(1-(r['rank']-1)/len(PD))
    return (f"  {a:<12} x {b:<18} r={r['r']:+.3f} p={r['p']:.4f} n={r['n']}  "
            f"rank {int(r['rank'])}/{len(PD)} ({pct:.0f}th pct of all pairs)")
print("\n--- PRE-REGISTERED COMPARISONS ---")
for b in ["OPC","Oligodendrocytes"]:
    print(show("IN-PV",b))
print("\n--- context: other interneuron x OL pairs (should be similar if non-specific) ---")
for a in ["IN-SST","IN-VIP","IN-SV2C"]:
    for b in ["OPC","Oligodendrocytes"]:
        print(show(a,b))
print("\n--- context: excitatory x OL ---")
for a in ["L2/3","L4","L5/6"]:
    print(show(a,"OPC"))
print("\nTOP 8 pairs overall:")
print(PD.head(8)[["a","b","r","p","n"]].to_string(index=False))
D.to_csv(OUT/"dysregulation_scores.csv",index=False); PD.to_csv(OUT/"pair_correlations.csv",index=False)
print("\nwrote outputs/testB/{dysregulation_scores,pair_correlations}.csv")

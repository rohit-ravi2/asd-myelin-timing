#!/usr/bin/env python3
"""TEST B step 1 — pseudobulk Velmeshev 2019 (65,217 genes x 104,559 cells,
238M nonzeros) into individual x cell-type sums by streaming the .mtx."""
import numpy as np, pandas as pd, json
from pathlib import Path
V=Path("data/dvad-corpus/asd_scope/velmeshev2019/raw"); P=Path(".")
OUT=P/"outputs/testB"; OUT.mkdir(parents=True, exist_ok=True)
genes=pd.read_csv(V/"genes.tsv",sep="\t",header=None,names=["ens","sym"])
barcodes=pd.read_csv(V/"barcodes.tsv",header=None,names=["cell"])
meta=pd.read_csv(V/"meta.txt",sep="\t")
print(f"genes={len(genes)} barcodes={len(barcodes)} meta={len(meta)}")
meta=meta.set_index("cell")
m=meta.reindex(barcodes.cell)
print("meta matched to barcode order:", m.diagnosis.notna().sum(), "of", len(m))
grp=(m.individual.astype(str)+"||"+m.cluster.astype(str)).values
keys=sorted(set(grp[pd.notna(m.cluster.values)]))
kidx={k:i for i,k in enumerate(keys)}
cell2g=np.array([kidx.get(g,-1) for g in grp],dtype=np.int32)
print(f"groups (individual x celltype) = {len(keys)}")
PB=np.zeros((len(genes),len(keys)),dtype=np.float32)
ncell=np.zeros(len(keys),dtype=np.int64)
for i,g in enumerate(cell2g):
    if g>=0: ncell[g]+=1
tot=0
reader=pd.read_csv(V/"matrix.mtx",sep=r"\s+",skiprows=3,header=None,
                   names=["g","c","v"],dtype={"g":np.int32,"c":np.int32,"v":np.float32},
                   chunksize=25_000_000)
for ci,ch in enumerate(reader):
    gi=ch.g.values-1; cc=ch.c.values-1; vv=ch.v.values
    gg=cell2g[cc]; ok=gg>=0
    np.add.at(PB,(gi[ok],gg[ok]),vv[ok])
    tot+=len(ch); print(f"  chunk {ci+1}: {tot:,} nonzeros", flush=True)
np.savez_compressed(OUT/"pseudobulk.npz",PB=PB,ncell=ncell,
                    keys=np.array(keys),genes=genes.sym.values)
info=m.reset_index()[["individual","diagnosis","region","age","sex"]].drop_duplicates("individual")
info.to_csv(OUT/"individuals.csv",index=False)
print(f"\nwrote {OUT}/pseudobulk.npz  PB={PB.shape}  total nonzeros {tot:,}")
print(info.diagnosis.value_counts().to_dict())

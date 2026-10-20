#!/usr/bin/env python3
"""When are the surviving gene sets expressed? The oligodendrocyte arm survived
every genetic control. If ASD-associated OL genes peak PRENATALLY that implies
OPC specification; if POSTNATALLY, myelination proper. Different mechanisms,
different intervention windows."""
import csv, json, re
import numpy as np
from pathlib import Path
B=Path("data/asd-corpus/brainspan"); C=Path("data/asd-corpus")
P=Path(__file__).resolve().parent.parent
cols=list(csv.DictReader(open(B/"columns_metadata.csv")))
rows=list(csv.DictReader(open(B/"rows_metadata.csv")))
sym=[r["gene_symbol"] for r in rows]
X=np.loadtxt(B/"expression_matrix.csv", delimiter=",")[:,1:]
print(f"expression matrix {X.shape}  samples={len(cols)}")
def age_pcw(a):
    m=re.match(r"(\d+)\s*pcw",a)
    if m: return float(m.group(1))
    m=re.match(r"(\d+)\s*mos",a)
    if m: return 40+float(m.group(1))*4.3
    m=re.match(r"(\d+)\s*yrs",a)
    if m: return 40+float(m.group(1))*52
    return np.nan
age=np.array([age_pcw(c["age"]) for c in cols])
cortex=np.array([c["structure_acronym"] not in ("CB","CBC","MD","STR","AMY","HIP") for c in cols])
prenatal=(age<40)&cortex; postnatal=(age>=40)&cortex
print(f"cortical samples: prenatal={prenatal.sum()}  postnatal={postnatal.sum()}")
gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
g=lambda t:gmt.get(t,set())
e2s={}
for line in open(C/"tools/magma_geneloc/NCBI37.3.gene.loc"):
    f=line.split()
    if len(f)>=6: e2s[f[0]]=f[5]
se={}
for line in open(P/"annot/asd_genesets.magma.set"):
    f=line.rstrip("\n").split("\t"); se[f[0]]=set(f[1:])
SFARI={e2s[x] for x in se["SFARI_ALL"] if x in e2s}
OL=g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOBP_ENSHEATHMENT_OF_NEURONS")
SETS={"OL_all":OL,"OL_x_SFARI":OL&SFARI,"SFARI_all":SFARI,
      "APOPTOSIS":g("GOBP_APOPTOTIC_PROCESS"),"APOP_x_SFARI":g("GOBP_APOPTOTIC_PROCESS")&SFARI,
      "MYELIN_CC":g("GOCC_MYELIN_SHEATH")}
lgX=np.log2(X+1)
idx={s:i for i,s in enumerate(sym)}
print(f"\n{'set':<16}{'n':>5}{'prenatal':>10}{'postnatal':>11}{'post/pre':>10}  peak window")
out=[]
for sn,ss in SETS.items():
    ii=[idx[s] for s in ss if s in idx]
    if len(ii)<5: continue
    prof=lgX[ii].mean(axis=0)
    pre=prof[prenatal].mean(); post=prof[postnatal].mean()
    # developmental trajectory in cortical samples
    order=np.argsort(age[cortex]); a=age[cortex][order]; v=prof[cortex][order]
    peak=a[np.argmax(np.convolve(v,np.ones(9)/9,mode="same"))]
    pk = f"{peak:.0f} pcw (prenatal)" if peak<40 else f"{(peak-40)/52:.1f} yr (postnatal)"
    print(f"{sn:<16}{len(ii):>5}{pre:>10.2f}{post:>11.2f}{post-pre:>10.2f}  {pk}")
    out.append({"set":sn,"n":len(ii),"prenatal_mean":round(float(pre),3),
                "postnatal_mean":round(float(post),3),"log2_post_minus_pre":round(float(post-pre),3),
                "peak_pcw":float(peak)})
json.dump(out,open(P/"outputs/apoptosis/brainspan_timing.json","w"),indent=1)
print("\nwrote outputs/apoptosis/brainspan_timing.json")

#!/usr/bin/env python3
"""The 13 intersection genes (CTNNB1 GLI3 HRAS MAPK3 MECP2 MTOR NF1 NOTCH1 NTRK2
PTEN SCN2A SOD1) look like pleiotropic developmental signalling hubs, not a
specific apoptosis-myelin module. Pleiotropic genes are annotated to many GO
terms, so they land in ANY two sets' intersection by construction.

Decisive control: match the permutation null on GO-annotation count as well as
LOEUF and CDS length. If the intersection enrichment is pleiotropy, it dies."""
import gzip, json, random
import numpy as np
from pathlib import Path
C=Path("data/asd-corpus"); P=Path(__file__).resolve().parent.parent
e2s={}
for line in open(C/"tools/magma_geneloc/NCBI37.3.gene.loc"):
    f=line.split()
    if len(f)>=6: e2s[f[0]]=f[5]
se={}
for line in open(P/"annot/asd_genesets.magma.set"):
    f=line.rstrip("\n").split("\t"); se[f[0]]=set(f[1:])
tos=lambda s:{e2s[g] for g in s if g in e2s}
LISTS={"SFARI_ALL":tos(se["SFARI_ALL"]),"SFARI_SCORE1":tos(se["SFARI_SCORE1"])}
gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
g=lambda t:gmt.get(t,set())
APOP=g("GOBP_APOPTOTIC_PROCESS"); OL=g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOBP_ENSHEATHMENT_OF_NEURONS")
INTER=APOP&OL
# pleiotropy = number of GO sets each gene belongs to
nterm={}
for t,genes in gmt.items():
    for x in genes: nterm[x]=nterm.get(x,0)+1
loeuf,cds={},{}
with gzip.open(C/"genetics/gnomad_v2.1.1_lof_metrics.by_gene.txt.bgz","rt") as f:
    h=f.readline().rstrip("\n").split("\t"); gi,li,ci=h.index("gene"),h.index("oe_lof_upper"),h.index("cds_length")
    for line in f:
        p=line.rstrip("\n").split("\t")
        try: loeuf[p[gi]]=float(p[li]); cds[p[gi]]=float(p[ci])
        except Exception: pass
UNIV=sorted(set(loeuf)&set(cds)&set(nterm))
print(f"universe with LOEUF+CDS+GO-count: {len(UNIV)}")
lq={x:i for i,x in enumerate(UNIV)}
def dig(vals,k=5):
    a=np.array(vals); return np.digitize(a,np.quantile(a,np.linspace(0,1,k+1)[1:-1]))
lo=dig([loeuf[x] for x in UNIV]); cd=dig([cds[x] for x in UNIV]); nt=dig([nterm[x] for x in UNIV])
strata={}
for i,x in enumerate(UNIV): strata.setdefault((lo[i],cd[i],nt[i]),[]).append(x)
rng=random.Random(777)
def test(genes,sset,n_perm=4000):
    tgt=[x for x in genes if x in lq]; obs=len(set(tgt)&sset); comp={}
    for x in tgt:
        k=(lo[lq[x]],cd[lq[x]],nt[lq[x]]); comp[k]=comp.get(k,0)+1
    null=[]
    for _ in range(n_perm):
        s=set()
        for k,n in comp.items():
            pool=strata.get(k,[])
            if len(pool)>=n: s.update(rng.sample(pool,n))
        null.append(len(s&sset))
    null=np.array(null); e=null.mean()
    return obs,float(e),float(obs/e) if e else float("nan"),float((null>=obs).mean())
print("\nGO-annotation-count MATCHED null (pleiotropy controlled):")
print(f"{'list':<14}{'set':<16}{'hit':>5}{'exp':>8}{'ratio':>7}{'p':>9}")
out=[]
for ln,L in LISTS.items():
    for sn,ss in [("APOPTOSIS",APOP),("OLIGO",OL),("INTERSECTION",INTER)]:
        o,e,r,p=test(L,ss)
        star="***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
        print(f"{ln:<14}{sn:<16}{o:>5}{e:>8.1f}{r:>7.2f}{p:>9.4f} {star}")
        out.append({"list":ln,"set":sn,"hit":o,"exp":round(e,1),"ratio":round(r,2),"p":round(p,4)})
med_int=np.median([nterm.get(x,0) for x in INTER if x in nterm])
med_all=np.median([nterm[x] for x in UNIV])
print(f"\npleiotropy check: median GO-term count — intersection genes {med_int:.0f} vs genome median {med_all:.0f}")
json.dump(out,open(P/"outputs/apoptosis/pleiotropy_control.json","w"),indent=1)

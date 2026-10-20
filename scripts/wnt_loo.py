#!/usr/bin/env python3
"""Leave-one-out on the two surviving results. Everything earlier reduced to
CTNNB1; does this too?"""
import csv, gzip, random
import numpy as np
from pathlib import Path
C=Path("data/asd-corpus"); S=Path("data/dvad-corpus/asd_scope/rare_variant")
P=Path(".")
e2s={}
for line in open(C/"tools/magma_geneloc/NCBI37.3.gene.loc"):
    f=line.split()
    if len(f)>=6: e2s[f[0]]=f[5]
se={}
for line in open(P/"annot/asd_genesets.magma.set"):
    f=line.rstrip("\n").split("\t"); se[f[0]]=set(f[1:])
tos=lambda s:{e2s[g] for g in s if g in e2s}
sat={r["hugoGene"] for r in csv.DictReader(open(S/"Satterstrom2020_102genes.csv")) if r.get("hugoGene")}
fu={r["HGNC"] for r in csv.DictReader(open(S/"Fu2022_72_studywide_genes.csv")) if r.get("HGNC")}
LISTS={"SFARI_ALL":tos(se["SFARI_ALL"]),"SFARI_SCORE1":tos(se["SFARI_SCORE1"]),"GOLD_127":sat|fu}
gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
g=lambda t:gmt.get(t,set())
NEG=g("GOBP_NEGATIVE_REGULATION_OF_OLIGODENDROCYTE_DIFFERENTIATION")
POS=g("GOBP_POSITIVE_REGULATION_OF_OLIGODENDROCYTE_DIFFERENTIATION")
BCAT=g("GOCC_BETA_CATENIN_TCF_COMPLEX")
nterm={}
for t,gs in gmt.items():
    for x in gs: nterm[x]=nterm.get(x,0)+1
loeuf,cds={},{}
with gzip.open(C/"genetics/gnomad_v2.1.1_lof_metrics.by_gene.txt.bgz","rt") as f:
    h=f.readline().rstrip("\n").split("\t"); gi,li,ci=h.index("gene"),h.index("oe_lof_upper"),h.index("cds_length")
    for line in f:
        p=line.rstrip("\n").split("\t")
        try: loeuf[p[gi]]=float(p[li]); cds[p[gi]]=float(p[ci])
        except Exception: pass
UNIV=sorted(set(loeuf)&set(cds)&set(nterm)); gi_={x:i for i,x in enumerate(UNIV)}
def dig(v,k=5):
    a=np.array(v); return np.digitize(a,np.quantile(a,np.linspace(0,1,k+1)[1:-1]))
lo=dig([loeuf[x] for x in UNIV]); cd=dig([cds[x] for x in UNIV]); nt=dig([nterm[x] for x in UNIV])
strata={}
for i,x in enumerate(UNIV): strata.setdefault((lo[i],cd[i],nt[i]),[]).append(x)
def run(genes,seed,n=4000):
    rng=random.Random(seed)
    tgt=set(x for x in genes if x in gi_); comp={}
    for x in tgt:
        i=gi_[x]; k=(lo[i],cd[i],nt[i]); comp[k]=comp.get(k,0)+1
    D=[]
    for _ in range(n):
        s=set()
        for k,c_ in comp.items():
            pool=strata.get(k,[])
            if len(pool)>=c_: s.update(rng.sample(pool,c_))
        D.append(s)
    nN,nP=len(NEG&set(UNIV)),len(POS&set(UNIV))
    obsd=len(tgt&NEG)/nN-len(tgt&POS)/nP
    pd_=float(np.mean([ (len(d&NEG)/nN-len(d&POS)/nP)>=obsd for d in D]))
    ob=len(tgt&BCAT); pb=float(np.mean([len(d&BCAT)>=ob for d in D]))
    eb=float(np.mean([len(d&BCAT) for d in D]))
    return pd_,ob,eb,pb
for ln,L in LISTS.items():
    pd0,ob0,eb0,pb0=run(L,7)
    print(f"\n=== {ln} ===")
    print(f"  full: direction p={pd0:.4f} | beta-catenin/TCF {ob0} vs {eb0:.1f} exp, p={pb0:.4f}")
    bh=sorted(set(x for x in L if x in gi_)&BCAT)
    nh=sorted(set(x for x in L if x in gi_)&NEG)
    print(f"  beta-catenin/TCF genes: {bh}")
    print(f"  LOO on beta-catenin/TCF:")
    for dg in bh:
        _,o,e,p=run(set(L)-{dg},7)
        print(f"     drop {dg:<9} -> {o} vs {e:.1f} exp, p={p:.4f}  {'SURVIVES' if p<0.05 else 'DIES'}")
    print(f"  LOO on direction test (NEG genes {nh}):")
    for dg in nh:
        p,_,_,_=run(set(L)-{dg},7)
        print(f"     drop {dg:<9} -> direction p={p:.4f}  {'SURVIVES' if p<0.05 else 'DIES'}")

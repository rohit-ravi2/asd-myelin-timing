#!/usr/bin/env python3
"""Is APOP x OLIGO a genuine CONVERGENCE, or just two independent enrichments
that happen to overlap? Under independence the intersection enrichment should be
~ the product of the marginals. Convergence requires it to EXCEED that."""
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
APOP=g("GOBP_APOPTOTIC_PROCESS")
OL=g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOBP_ENSHEATHMENT_OF_NEURONS")
loeuf,cds={},{}
with gzip.open(C/"genetics/gnomad_v2.1.1_lof_metrics.by_gene.txt.bgz","rt") as f:
    h=f.readline().rstrip("\n").split("\t"); gi,li,ci=h.index("gene"),h.index("oe_lof_upper"),h.index("cds_length")
    for line in f:
        p=line.rstrip("\n").split("\t")
        try: loeuf[p[gi]]=float(p[li]); cds[p[gi]]=float(p[ci])
        except Exception: pass
UNIV=sorted(set(loeuf)&set(cds)); lq={x:i for i,x in enumerate(UNIV)}
lo=np.array([loeuf[x] for x in UNIV]); cd=np.array([cds[x] for x in UNIV])
ld=np.digitize(lo,np.quantile(lo,np.linspace(0,1,11)[1:-1]))
cdd=np.digitize(cd,np.quantile(cd,np.linspace(0,1,11)[1:-1]))
strata={}
for i,x in enumerate(UNIV): strata.setdefault((ld[i],cdd[i]),[]).append(x)
rng=random.Random(99)
def draws(genes,n_perm=4000):
    tgt=[x for x in genes if x in lq]; comp={}
    for x in tgt:
        k=(ld[lq[x]],cdd[lq[x]]); comp[k]=comp.get(k,0)+1
    D=[]
    for _ in range(n_perm):
        s=set()
        for k,n in comp.items():
            pool=strata.get(k,[])
            if len(pool)>=n: s.update(rng.sample(pool,n))
        D.append(s)
    return set(tgt),D
for lname,L in LISTS.items():
    tgt,D=draws(L)
    nA=len(tgt&APOP); nO=len(tgt&OL); nI=len(tgt&APOP&OL)
    eA=np.mean([len(d&APOP) for d in D]); eO=np.mean([len(d&OL) for d in D])
    eI=np.mean([len(d&APOP&OL) for d in D])
    rA,rO,rI=nA/eA,nO/eO,nI/eI
    pI=np.mean([len(d&APOP&OL)>=nI for d in D])
    print(f"\n=== {lname} (n={len(tgt)}) ===")
    print(f"  APOPTOSIS      {nA:>4} obs / {eA:>6.1f} exp  -> {rA:.2f}x")
    print(f"  OLIGO          {nO:>4} obs / {eO:>6.1f} exp  -> {rO:.2f}x")
    print(f"  INTERSECTION   {nI:>4} obs / {eI:>6.1f} exp  -> {rI:.2f}x   (marginal p={pI:.4f})")
    print(f"  --> independence predicts {rA*rO:.2f}x ; observed {rI:.2f}x "
          f"-> {'SUPER-additive (convergence)' if rI>rA*rO else 'NOT super-additive: two independent enrichments overlapping'}")
    # conditional test: intersection enrichment WITHIN the apoptosis-annotated ASD genes
    sub=tgt&APOP
    subd=[d&APOP for d in D]
    oIn=len(sub&OL); eIn=np.mean([len(s&OL) for s in subd])
    pIn=np.mean([len(s&OL)>=oIn for s in subd])
    print(f"  conditional: among ASD genes already apoptosis-annotated, OL hits "
          f"{oIn} vs {eIn:.1f} exp -> {oIn/eIn if eIn else float('nan'):.2f}x, p={pIn:.4f}")
    inter=sorted(tgt&APOP&OL); print(f"  genes: {inter}")
    if len(inter)<=20:
        surv=0
        for dgene in inter:
            t2=tgt-{dgene}
            o2=len(t2&APOP&OL); p2=np.mean([len(d&APOP&OL)>=o2 for d in D])
            if p2<0.05: surv+=1
        print(f"  leave-one-out: {surv}/{len(inter)} removals still p<0.05")

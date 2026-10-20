#!/usr/bin/env python3
"""REVERSE APPROACH: don't test pre-chosen sets. Scan EVERY GO term against the
ASD gene list with the pleiotropy-controlled null, and see what actually comes
out on top. One set of matched draws is reused across all terms, so this is fast
and the multiple-testing correction is honest."""
import gzip, json, random
import numpy as np
from scipy import sparse
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
TARGETS={"SFARI_ALL":tos(se["SFARI_ALL"]),"SFARI_SCORE1":tos(se["SFARI_SCORE1"])}

gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
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
UNIV=sorted(set(loeuf)&set(cds)&set(nterm))
gi_={x:i for i,x in enumerate(UNIV)}; n=len(UNIV)
TERMS=[t for t,gs in gmt.items() if 10 <= len(gs & set(UNIV)) <= 500]
print(f"universe={n} genes, scanning {len(TERMS)} GO terms (size 10-500)")
rowi,coli=[],[]
for ti,t in enumerate(TERMS):
    for x in gmt[t]:
        if x in gi_: rowi.append(ti); coli.append(gi_[x])
M=sparse.csr_matrix((np.ones(len(rowi),dtype=np.float32),(rowi,coli)),shape=(len(TERMS),n))

def dig(v,k=5):
    a=np.array(v); return np.digitize(a,np.quantile(a,np.linspace(0,1,k+1)[1:-1]))
lo=dig([loeuf[x] for x in UNIV]); cd=dig([cds[x] for x in UNIV]); nt=dig([nterm[x] for x in UNIV])
strata={}
for i,x in enumerate(UNIV): strata.setdefault((lo[i],cd[i],nt[i]),[]).append(i)
rng=random.Random(2024)
NPERM=2000
results={}
for tname,TG in TARGETS.items():
    idx=[gi_[x] for x in TG if x in gi_]
    comp={}
    for i in idx:
        k=(lo[i],cd[i],nt[i]); comp[k]=comp.get(k,0)+1
    v=np.zeros(n,dtype=np.float32); v[idx]=1
    obs=np.asarray(M@v).ravel()
    NULL=np.zeros((len(TERMS),NPERM),dtype=np.float32)
    for j in range(NPERM):
        w=np.zeros(n,dtype=np.float32)
        for k,c_ in comp.items():
            pool=strata.get(k,[])
            if len(pool)>=c_: w[rng.sample(pool,c_)]=1
        NULL[:,j]=np.asarray(M@w).ravel()
    exp=NULL.mean(axis=1)
    pv=(NULL>=obs[:,None]).mean(axis=1)
    order=np.argsort(pv)
    m=len(TERMS); q=np.minimum.accumulate((pv[order]*m/np.arange(1,m+1))[::-1])[::-1]
    qq=np.empty(m); qq[order]=np.clip(q,0,1)
    allrows=[{"term":TERMS[i],"rank":int(r+1),"hit":int(obs[i]),"exp":round(float(exp[i]),1),"ratio":round(float(obs[i]/exp[i]),2) if exp[i]>0 else None,"p":round(float(pv[i]),5),"q":round(float(qq[i]),4)} for r,i in enumerate(order)]
    results[tname+"_FULL"]=allrows
    rows=[{"term":TERMS[i],"hit":int(obs[i]),"exp":round(float(exp[i]),1),
           "ratio":round(float(obs[i]/exp[i]),2) if exp[i]>0 else None,
           "p":round(float(pv[i]),5),"q":round(float(qq[i]),4)} for i in order[:40]]
    results[tname]=rows
    print(f"\n{'='*100}\nTOP 25 GO TERMS — {tname} (n={len(idx)}), pleiotropy-controlled, BH over {m} terms\n{'='*100}")
    print(f"{'rank':<5}{'term':<62}{'hit':>5}{'exp':>7}{'ratio':>7}{'q':>9}")
    for r,row in enumerate(rows[:25],1):
        print(f"{r:<5}{row['term'][:60]:<62}{row['hit']:>5}{row['exp']:>7.1f}{row['ratio']:>7.2f}{row['q']:>9.4f}")
    print(f"  (significant at q<0.05: {int((qq<0.05).sum())} terms)")
json.dump(results,open(P/"outputs/apoptosis/unbiased_go_scan_full.json","w"))
print("\nwrote outputs/apoptosis/unbiased_go_scan.json")

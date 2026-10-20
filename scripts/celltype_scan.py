#!/usr/bin/env python3
"""Which cell type do ASD risk genes belong to? Unbiased scan across all 34
Velmeshev 2019 human-cortex cell-type marker sets, pleiotropy-controlled."""
import gzip, glob, os, random
import numpy as np
from pathlib import Path
C=Path("data/asd-corpus"); S="data/dvad-corpus/asd_scope/velmeshev2019/markers"
P=Path(".")
e2s={}
for line in open(C/"tools/magma_geneloc/NCBI37.3.gene.loc"):
    f=line.split()
    if len(f)>=6: e2s[f[0]]=f[5]
se={}
for line in open(P/"annot/asd_genesets.magma.set"):
    f=line.rstrip("\n").split("\t"); se[f[0]]=set(f[1:])
tos=lambda s:{e2s[g] for g in s if g in e2s}
LISTS={"SFARI_SCORE1":tos(se["SFARI_SCORE1"]),"SFARI_ALL":tos(se["SFARI_ALL"])}
cells={}
for f in sorted(glob.glob(S+"/*.tsv.gz")):
    n=os.path.basename(f).replace(".tsv.gz","")
    s=set()
    with gzip.open(f,"rt") as fh:
        fh.readline()
        for line in fh:
            p=line.split("\t")
            if len(p)>1: s.add(p[1])
    if len(s)>=25: cells[n]=s
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
UNIV=sorted(set(loeuf)&set(cds)&set(nterm)); gi_={x:i for i,x in enumerate(UNIV)}
def dig(v,k=5):
    a=np.array(v); return np.digitize(a,np.quantile(a,np.linspace(0,1,k+1)[1:-1]))
lo=dig([loeuf[x] for x in UNIV]); cd=dig([cds[x] for x in UNIV]); nt=dig([nterm[x] for x in UNIV])
strata={}
for i,x in enumerate(UNIV): strata.setdefault((lo[i],cd[i],nt[i]),[]).append(x)
for ln,L in LISTS.items():
    rng=random.Random(808)
    tgt=set(x for x in L if x in gi_); comp={}
    for x in tgt:
        i=gi_[x]; k=(lo[i],cd[i],nt[i]); comp[k]=comp.get(k,0)+1
    D=[]
    for _ in range(4000):
        s=set()
        for k,c_ in comp.items():
            pool=strata.get(k,[])
            if len(pool)>=c_: s.update(rng.sample(pool,c_))
        D.append(s)
    rows=[]
    for cn,cs in cells.items():
        o=len(tgt&cs); e=float(np.mean([len(d&cs) for d in D]))
        p=float(np.mean([len(d&cs)>=o for d in D]))
        rows.append((p,cn,o,e,(o/e if e else float('nan')),len(cs)))
    rows.sort()
    m=len(rows); q=[]
    for i,(p,*_ ) in enumerate(rows,1): q.append(min(1,p*m/i))
    for i in range(m-2,-1,-1): q[i]=min(q[i],q[i+1])
    print(f"\n{'='*86}\n{ln} (n={len(tgt)}) — cell-type marker enrichment, pleiotropy-controlled, BH over {m}\n{'='*86}")
    print(f"{'cell type':<22}{'nmark':>6}{'hit':>5}{'exp':>7}{'ratio':>7}{'p':>9}{'q':>8}")
    for i,(p,cn,o,e,r,nm) in enumerate(rows[:12]):
        star="***" if q[i]<0.001 else "**" if q[i]<0.01 else "*" if q[i]<0.05 else ""
        print(f"{cn:<22}{nm:>6}{o:>5}{e:>7.1f}{r:>7.2f}{p:>9.4f}{q[i]:>8.3f} {star}")

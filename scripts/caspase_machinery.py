#!/usr/bin/env python3
"""Literature-directed refinement. Erturk/Wang/Sheng 2014 (J Neurosci 34:1672)
showed local caspase-3 prunes spines WITHOUT cell death, spatially confined by
IAP and proteasome 'molecular brakes'. So the relevant gene set is not the
1,972-gene GO apoptosis term but the executioner caspases + their brakes."""
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
SETS={
 "CASPASE_EXECUTIONER": {"CASP3","CASP6","CASP7","CASP9","APAF1","CYCS","BAX","BAK1","BID","DIABLO"},
 "IAP_BRAKES":          {"XIAP","BIRC2","BIRC3","BIRC5","BIRC6","BIRC7","NAIP"},
 "CASPASE_PLUS_BRAKES": {"CASP3","CASP6","CASP7","CASP9","APAF1","CYCS","BAX","BAK1","BID","DIABLO",
                         "XIAP","BIRC2","BIRC3","BIRC5","BIRC6","BIRC7","NAIP"},
 "COMPLEMENT_PRUNING":  {"C1QA","C1QB","C1QC","C3","C4A","C4B","CR3","ITGAM","CX3CR1","TREM2","TYROBP","MERTK","GAS6"},
 "GO_APOPTOSIS_broad":  g("GOBP_APOPTOTIC_PROCESS"),
 "OLIGO_functional":    g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOBP_ENSHEATHMENT_OF_NEURONS"),
}
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
UNIV=sorted(set(loeuf)&set(cds)&set(nterm)); lq={x:i for i,x in enumerate(UNIV)}
def dig(v,k=5):
    a=np.array(v); return np.digitize(a,np.quantile(a,np.linspace(0,1,k+1)[1:-1]))
lo=dig([loeuf[x] for x in UNIV]); cd=dig([cds[x] for x in UNIV]); nt=dig([nterm[x] for x in UNIV])
strata={}
for i,x in enumerate(UNIV): strata.setdefault((lo[i],cd[i],nt[i]),[]).append(x)
rng=random.Random(4242)
def test(genes,sset,n=4000):
    tgt=[x for x in genes if x in lq]; obs=len(set(tgt)&sset); comp={}
    for x in tgt:
        k=(lo[lq[x]],cd[lq[x]],nt[lq[x]]); comp[k]=comp.get(k,0)+1
    null=[]
    for _ in range(n):
        s=set()
        for k,c_ in comp.items():
            pool=strata.get(k,[])
            if len(pool)>=c_: s.update(rng.sample(pool,c_))
        null.append(len(s&sset))
    null=np.array(null); e=null.mean()
    return obs,float(e),(obs/e if e else float("nan")),float((null>=obs).mean())
print("Pleiotropy-controlled (LOEUF x CDS x GO-annotation-count matched):\n")
print(f"{'list':<14}{'set':<24}{'hit':>5}{'exp':>7}{'ratio':>7}{'p':>9}")
out=[]
for ln,L in LISTS.items():
    for sn,ss in SETS.items():
        o,e,r,p=test(L,ss)
        star="***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
        print(f"{ln:<14}{sn:<24}{o:>5}{e:>7.1f}{r:>7.2f}{p:>9.4f} {star}")
        out.append({"list":ln,"set":sn,"hit":o,"exp":round(e,1),"ratio":round(r,2),"p":round(p,4)})
    print()
for ln,L in LISTS.items():
    hits=sorted(set(x for x in L if x in lq)&SETS["CASPASE_PLUS_BRAKES"])
    hc=sorted(set(x for x in L if x in lq)&SETS["COMPLEMENT_PRUNING"])
    print(f"{ln}: caspase/brake hits={hits}  complement hits={hc}")
json.dump(out,open(P/"outputs/apoptosis/caspase_machinery.json","w"),indent=1)

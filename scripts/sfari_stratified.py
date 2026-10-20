#!/usr/bin/env python3
"""Does the apoptosis x oligodendrocyte convergence survive at scale?
Gold-standard n=127 gave 4 driving genes. SFARI_ALL (~1142) is ~9x larger; if
the convergence is a general ASD mechanism it must appear here robustly."""
import gzip, json, random
import numpy as np
from pathlib import Path
C=Path("data/asd-corpus"); P=Path(__file__).resolve().parent.parent

e2s={}
for line in open(C/"tools/magma_geneloc/NCBI37.3.gene.loc"):
    f=line.split()
    if len(f)>=6: e2s[f[0]]=f[5]
sets_entrez={}
for line in open(P/"annot/asd_genesets.magma.set"):
    f=line.rstrip("\n").split("\t"); sets_entrez[f[0]]=set(f[1:])
tosym=lambda s:{e2s[g] for g in s if g in e2s}
SFARI_ALL=tosym(sets_entrez["SFARI_ALL"]); SFARI1=tosym(sets_entrez["SFARI_SCORE1"])

gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt","c2.cp.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
g=lambda t:gmt.get(t,set())
APOP=g("GOBP_APOPTOTIC_PROCESS")
OL=g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOBP_ENSHEATHMENT_OF_NEURONS")
MTOR={"PTEN","TSC1","TSC2","NF1","MTOR","AKT1","AKT2","AKT3","RPTOR","RICTOR","PIK3CA",
      "PIK3R2","DEPDC5","STRADA","RHEB","RPS6KB1","EIF4EBP1","NPRL2","NPRL3","FLCN"}
TESTS={"APOPTOSIS":APOP,"OLIGO_functional":OL,"APOP_x_OLIGO":APOP&OL,"MTORopathy":MTOR,
       "PRUNING_AUTOPHAGY":g("GOBP_MACROAUTOPHAGY")|g("GOBP_SYNAPSE_PRUNING"),
       "CHROMATIN_posctrl":g("GOBP_CHROMATIN_REMODELING"),"SYNAPSE_posctrl":g("GOCC_SYNAPSE")}

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
rng=random.Random(5)
def mp(genes,sset,n_perm=3000):
    tgt=[x for x in genes if x in lq]; obs=len(set(tgt)&sset)
    comp={}
    for x in tgt:
        k=(ld[lq[x]],cdd[lq[x]]); comp[k]=comp.get(k,0)+1
    null=[]
    for _ in range(n_perm):
        s=set()
        for k,n in comp.items():
            pool=strata.get(k,[])
            if len(pool)>=n: s.update(rng.sample(pool,n))
        null.append(len(s&sset))
    null=np.array(null); return obs,float(null.mean()),float((null>=obs).mean())

rows=[]
for lname,L in [("SFARI_ALL",SFARI_ALL),("SFARI_SCORE1",SFARI1)]:
    print(f"\n=== {lname} (n={len([x for x in L if x in lq])}) ===")
    print(f"  {'set':<22}{'hit':>5}{'exp':>8}{'ratio':>7}{'p':>9}")
    for s,ss in TESTS.items():
        o,e,p=mp(L,ss); rows.append({"list":lname,"set":s,"hit":o,"exp":round(e,1),
                                     "ratio":round(o/e,2) if e else None,"p":round(p,4)})
        star="***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
        print(f"  {s:<22}{o:>5}{e:>8.1f}{(o/e if e else 0):>7.2f}{p:>9.4f} {star}")
json.dump(rows,open(P/"outputs/apoptosis/sfari_stratified.json","w"),indent=1)
print("\nwrote outputs/apoptosis/sfari_stratified.json")

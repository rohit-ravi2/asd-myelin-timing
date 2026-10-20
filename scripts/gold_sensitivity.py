#!/usr/bin/env python3
"""Leave-one-gene-out sensitivity + fixed controls + BH. The APOP x OLIGO result
rests on 4 genes; if it dies when any one is removed it is not a finding."""
import csv, gzip, json, random
import numpy as np
from pathlib import Path
C=Path("data/asd-corpus"); S=Path("data/dvad-corpus/asd_scope/rare_variant")
P=Path(__file__).resolve().parent.parent
gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt","c2.cp.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
g=lambda t: gmt.get(t,set())
loeuf,cds={},{}
with gzip.open(C/"genetics/gnomad_v2.1.1_lof_metrics.by_gene.txt.bgz","rt") as f:
    h=f.readline().rstrip("\n").split("\t"); gi,li,ci=h.index("gene"),h.index("oe_lof_upper"),h.index("cds_length")
    for line in f:
        p=line.rstrip("\n").split("\t")
        try: loeuf[p[gi]]=float(p[li]); cds[p[gi]]=float(p[ci])
        except Exception: pass
UNIV=sorted(set(loeuf)&set(cds))
lq={x:i for i,x in enumerate(UNIV)}
lo=np.array([loeuf[x] for x in UNIV]); cd=np.array([cds[x] for x in UNIV])
ld=np.digitize(lo,np.quantile(lo,np.linspace(0,1,11)[1:-1]))
cdd=np.digitize(cd,np.quantile(cd,np.linspace(0,1,11)[1:-1]))
strata={}
for i,x in enumerate(UNIV): strata.setdefault((ld[i],cdd[i]),[]).append(x)
rng=random.Random(23)
def mp(genes,sset,n_perm=5000):
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
    null=np.array(null)
    return obs,float(null.mean()),float((null>=obs).mean())

sat={r["hugoGene"] for r in csv.DictReader(open(S/"Satterstrom2020_102genes.csv")) if r.get("hugoGene")}
fu={r["HGNC"] for r in csv.DictReader(open(S/"Fu2022_72_studywide_genes.csv")) if r.get("HGNC")}
UNION=sat|fu
APOP=g("GOBP_APOPTOTIC_PROCESS")
OL_full=g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOBP_ENSHEATHMENT_OF_NEURONS")|g("GOCC_MYELIN_SHEATH")
OL_func=g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOBP_ENSHEATHMENT_OF_NEURONS")  # drop CC localisation
MTOR={"PTEN","TSC1","TSC2","NF1","MTOR","AKT1","AKT2","AKT3","RPTOR","RICTOR","PIK3CA",
      "PIK3R2","DEPDC5","STRADA","RHEB","RPS6KB1","EIF4EBP1","NPRL2","NPRL3","FLCN"}
SETS={"APOP_x_OLIGO_full":APOP&OL_full,"APOP_x_OLIGO_functional":APOP&OL_func,
      "OLIGO_full":OL_full,"MTORopathy":MTOR,
      "CHROMATIN_REMODELING_posctrl":g("GOBP_CHROMATIN_REMODELING"),
      "SYNAPSE_posctrl":g("GOCC_SYNAPSE")}
print("=== headline + FIXED controls (UNION, n=%d) ==="%len(UNION))
res={}
for s,ss in SETS.items():
    o,e,p=mp(UNION,ss); res[s]=p
    print(f"  {s:<30} hit={o:>3} exp={e:>5.1f} ratio={o/e if e else float('nan'):>5.2f} p={p:.4f}")
print("\n=== LEAVE-ONE-GENE-OUT: APOP x OLIGO (functional, SCN2A-artefact removed) ===")
driving=sorted(UNION & SETS["APOP_x_OLIGO_functional"])
print(f"  driving genes: {driving}")
for drop in driving:
    o,e,p=mp(UNION-{drop}, SETS["APOP_x_OLIGO_functional"])
    print(f"    drop {drop:<8} -> hit={o} exp={e:.1f} p={p:.4f}  {'SURVIVES' if p<0.05 else 'DIES'}")
print("\n=== LEAVE-ONE-OUT: mTORopathy ===")
for drop in sorted(UNION & MTOR):
    o,e,p=mp(UNION-{drop}, MTOR)
    print(f"    drop {drop:<8} -> hit={o} exp={e:.1f} p={p:.4f}  {'SURVIVES' if p<0.05 else 'DIES'}")
json.dump(res, open(P/"outputs/apoptosis/sensitivity.json","w"), indent=1)

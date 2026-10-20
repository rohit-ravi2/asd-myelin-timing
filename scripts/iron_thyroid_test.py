#!/usr/bin/env python3
"""Is there GENETIC predisposition in iron handling or thyroid-hormone signalling?
If maternal iron deficiency were the route to OL dysfunction, genetic variation in
iron transport/storage should also confer risk. Same pleiotropy-controlled null
that killed the apoptosis arm."""
import csv, gzip, random
import numpy as np
from pathlib import Path
C=Path("data/asd-corpus"); S=Path("data/dvad-corpus/asd_scope/rare_variant"); P=Path(".")
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
SETS={
 "IRON_curated": {"TFRC","TF","TFR2","SLC11A2","SLC40A1","HAMP","HFE","FTH1","FTL","SLC25A37",
                  "SLC25A28","STEAP3","ACO1","IREB2","CYBRD1","HEPH","CP","MCOLN1","FXN","ABCB7"},
 "IRON_GO": g("GOBP_IRON_ION_HOMEOSTASIS")|g("GOBP_IRON_ION_TRANSPORT")|g("GOBP_CELLULAR_IRON_ION_HOMEOSTASIS"),
 "THYROID_curated": {"DIO1","DIO2","DIO3","THRA","THRB","SLC16A2","SLC7A5","SLC7A8","TRHR","TSHR",
                     "TG","TPO","TSHB","PAX8","NKX2-1","FOXE1","IYD","SLC5A5","SERPINA7","TTR"},
 "THYROID_GO": g("GOBP_THYROID_HORMONE_GENERATION")|g("GOBP_CELLULAR_RESPONSE_TO_THYROID_HORMONE_STIMULUS")|g("GOBP_THYROID_HORMONE_METABOLIC_PROCESS"),
 "CHOLINE_FOLATE_B12": {"CHKA","CHKB","PCYT1A","PCYT1B","CHDH","BHMT","PEMT","SLC44A1","SLC5A7",
                        "MTHFR","MTR","MTRR","FOLR1","FOLR2","SLC19A1","DHFR","TCN2","MMACHC","MMAA","MMAB"},
 "OLIGO_positive_control": g("GOBP_ENSHEATHMENT_OF_NEURONS")|g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION"),
 "SYNAPSE_positive_control": g("GOCC_SYNAPSE"),
}
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
    rng=random.Random(606)
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
    print(f"\n=== {ln} (n={len(tgt)}) — pleiotropy-controlled ===")
    print(f"  {'set':<26}{'nset':>5}{'hit':>5}{'exp':>7}{'ratio':>7}{'p':>9}")
    for sn,ss in SETS.items():
        ss2=ss&set(UNIV)
        o=len(tgt&ss2); e=float(np.mean([len(d&ss2) for d in D])); p=float(np.mean([len(d&ss2)>=o for d in D]))
        star="***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
        print(f"  {sn:<26}{len(ss2):>5}{o:>5}{e:>7.1f}{(o/e if e else 0):>7.2f}{p:>9.4f} {star}")
        if o and sn in ("IRON_curated","IRON_GO","THYROID_curated","THYROID_GO","CHOLINE_FOLATE_B12"):
            print(f"      genes: {sorted(tgt&ss2)}")

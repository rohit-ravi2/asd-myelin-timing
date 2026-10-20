#!/usr/bin/env python3
"""Are the gold-standard ASD genes (Satterstrom 102, Fu 72) enriched for
apoptosis / oligodendrocyte / mTOR annotation, against a background matched on
LOEUF constraint and CDS length?

This is better powered than another burden run: it uses the highest-confidence
ASD gene list rather than genome-wide association signal.
"""
import csv, gzip, json, random
import numpy as np
from pathlib import Path
from scipy import stats
C = Path("data/asd-corpus"); S = Path("data/dvad-corpus/asd_scope/rare_variant")
P = Path(__file__).resolve().parent.parent

# ---- gold standard
sat = {r["hugoGene"] for r in csv.DictReader(open(S/"Satterstrom2020_102genes.csv")) if r.get("hugoGene")}
fu  = {r["HGNC"] for r in csv.DictReader(open(S/"Fu2022_72_studywide_genes.csv")) if r.get("HGNC")}
GOLD = {"SATTERSTROM102": sat, "FU72": fu, "UNION": sat | fu}

# ---- constraint + length for matching
loeuf, cds = {}, {}
with gzip.open(C/"genetics/gnomad_v2.1.1_lof_metrics.by_gene.txt.bgz","rt") as f:
    h=f.readline().rstrip("\n").split("\t"); gi,li,ci=h.index("gene"),h.index("oe_lof_upper"),h.index("cds_length")
    for line in f:
        p=line.rstrip("\n").split("\t")
        try: loeuf[p[gi]]=float(p[li]); cds[p[gi]]=float(p[ci])
        except Exception: pass
UNIV = sorted(set(loeuf) & set(cds))
print(f"universe with LOEUF+CDS: {len(UNIV)}")

# ---- gene sets
gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt","c2.cp.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
def g(t): return gmt.get(t,set())

APOP = g("GOBP_APOPTOTIC_PROCESS")
PRUNE = g("GOBP_MACROAUTOPHAGY") | g("GOBP_SYNAPSE_PRUNING") | g("GOBP_REGULATION_OF_SYNAPSE_PRUNING")
OL = g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION") | g("GOBP_ENSHEATHMENT_OF_NEURONS") | g("GOCC_MYELIN_SHEATH")
MTOR_CUR = {"PTEN","TSC1","TSC2","NF1","MTOR","AKT1","AKT2","AKT3","RPTOR","RICTOR","PIK3CA",
            "PIK3R2","DEPDC5","STRADA","RHEB","RPS6KB1","EIF4EBP1","NPRL2","NPRL3","FLCN"}
SETS = {
 "APOPTOSIS":        APOP,
 "PRUNING_AUTOPHAGY":PRUNE,
 "APOP_x_PRUNING":   APOP & PRUNE,
 "OLIGO_MYELIN":     OL,
 "APOP_x_OLIGO":     APOP & OL,
 "MTOR_KEGG":        g("KEGG_MTOR_SIGNALING_PATHWAY"),
 "MTORopathy_CURATED": MTOR_CUR,
 "CHROMATIN_POSCTRL":g("GOBP_CHROMATIN_ORGANIZATION"),
 "SYNAPTIC_POSCTRL": g("GOCC_SYNAPSE"),
}
rng = random.Random(23)
def matched_null(genes, n_perm=5000):
    """Sample background genes matched on LOEUF decile x CDS decile."""
    lq = {x: i for i,x in enumerate(UNIV)}
    lo = np.array([loeuf[x] for x in UNIV]); cd = np.array([cds[x] for x in UNIV])
    ld = np.digitize(lo, np.quantile(lo, np.linspace(0,1,11)[1:-1]))
    cdd= np.digitize(cd, np.quantile(cd, np.linspace(0,1,11)[1:-1]))
    strata={}
    for i,x in enumerate(UNIV): strata.setdefault((ld[i],cdd[i]),[]).append(x)
    tgt=[x for x in genes if x in lq]
    comp={}
    for x in tgt: comp[(ld[lq[x]],cdd[lq[x]])] = comp.get((ld[lq[x]],cdd[lq[x]]),0)+1
    draws=[]
    for _ in range(n_perm):
        s=set()
        for k,n in comp.items():
            pool=strata.get(k,[])
            if len(pool)>=n: s.update(rng.sample(pool,n))
        draws.append(s)
    return tgt, draws

print(f"\n{'gold list':<16}{'annotation set':<22}{'n_gold':>7}{'hit':>5}{'exp':>7}{'OR':>6}{'matched p':>11}")
rows=[]
for gname, gl in GOLD.items():
    tgt, draws = matched_null(gl)
    for sname, sset in SETS.items():
        obs = len(set(tgt) & sset)
        null = np.array([len(d & sset) for d in draws])
        exp = null.mean()
        p = (null >= obs).mean()
        orr = (obs/exp) if exp>0 else float("nan")
        rows.append({"gold":gname,"set":sname,"n_gold":len(tgt),"hit":obs,
                     "exp":round(float(exp),2),"ratio":round(float(orr),2),"matched_p":round(float(p),4)})
        star = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
        print(f"{gname:<16}{sname:<22}{len(tgt):>7}{obs:>5}{exp:>7.1f}{orr:>6.2f}{p:>11.4f} {star}")
json.dump(rows, open(P/"outputs/apoptosis/goldstandard_enrichment.json","w"), indent=1)
print("\nwrote outputs/apoptosis/goldstandard_enrichment.json")

#!/usr/bin/env python3
"""Direction-aware test. The claim is not 'OL genes are enriched' (already known)
but that ASD genes preferentially hit NEGATIVE regulators of oligodendrocyte
differentiation -> fewer mature OLs -> hypomyelination.

So the statistic is the RATE DIFFERENCE: hits_NEG/|NEG| - hits_POS/|POS|,
tested against the same pleiotropy-matched null. Enrichment of both directions
equally would mean no directional claim survives."""
import csv, gzip, json, random
import numpy as np
from pathlib import Path
C=Path("data/asd-corpus"); S=Path("data/dvad-corpus/asd_scope/rare_variant")
P=Path(__file__).resolve().parent.parent
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
# curated, literature-based (stated explicitly as OUR curation)
NEG_CUR={"LINGO1","ID2","ID4","HES1","HES5","NOTCH1","JAG1","SOX5","SOX6","GPR17",
         "RTN4","RTN4R","DLX1","DLX2","TCF7L2","CTNNB1","WNT3A","NKX2-2","BMP4","BMP2","SIRT2"}
POS_CUR={"OLIG1","OLIG2","SOX10","MYRF","ZEB2","YY1","ASCL1","QKI","NRG1","ERBB3",
         "MBP","PLP1","MAG","MOG","CNP","MOBP","SMARCA4","HDAC1","HDAC2","TCF7L1","ZNF191"}
PAIRS={
 "OL_differentiation_GO":(g("GOBP_NEGATIVE_REGULATION_OF_OLIGODENDROCYTE_DIFFERENTIATION"),
                          g("GOBP_POSITIVE_REGULATION_OF_OLIGODENDROCYTE_DIFFERENTIATION")),
 "MYELINATION_GO":       (g("GOBP_NEGATIVE_REGULATION_OF_MYELINATION"),
                          g("GOBP_POSITIVE_REGULATION_OF_MYELINATION")),
 "CURATED_literature":   (NEG_CUR, POS_CUR),
}
SOLO={"CANONICAL_WNT":g("GOBP_CANONICAL_WNT_SIGNALING_PATHWAY"),
      "NEG_REG_CANONICAL_WNT":g("GOBP_NEGATIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY"),
      "BETA_CATENIN_TCF":g("GOCC_BETA_CATENIN_TCF_COMPLEX")}
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
rng=random.Random(31337)
def draws(genes,n=5000):
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
    return tgt,D
res=[]
for ln,L in LISTS.items():
    tgt,D=draws(L)
    print(f"\n{'='*94}\n{ln}  (n={len(tgt)})\n{'='*94}")
    for sn,ss in SOLO.items():
        o=len(tgt&ss); e=np.mean([len(d&ss) for d in D]); p=np.mean([len(d&ss)>=o for d in D])
        print(f"  {sn:<28} hit={o:>3} exp={e:>5.1f} ratio={o/e if e else float('nan'):>5.2f} p={p:.4f}")
        res.append({"list":ln,"test":sn,"hit":o,"exp":round(float(e),1),"p":round(float(p),4)})
    for pn,(NEG,POS) in PAIRS.items():
        nN,nP=len(NEG&set(UNIV)),len(POS&set(UNIV))
        if nN<3 or nP<3: continue
        oN,oP=len(tgt&NEG),len(tgt&POS)
        obs=oN/nN - oP/nP
        null=np.array([len(d&NEG)/nN - len(d&POS)/nP for d in D])
        p=float((null>=obs).mean())
        eN=np.mean([len(d&NEG) for d in D]); eP=np.mean([len(d&POS) for d in D])
        print(f"  {pn:<28} NEG {oN}/{nN} (exp {eN:.1f})  POS {oP}/{nP} (exp {eP:.1f})")
        print(f"  {'':<28}   rate-difference test (NEG>POS): obs={obs:+.4f}  p={p:.4f} "
              f"{'<-- DIRECTIONAL' if p<0.05 else ''}")
        res.append({"list":ln,"test":pn+"_direction","neg_hit":oN,"neg_n":nN,"pos_hit":oP,
                    "pos_n":nP,"obs_rate_diff":round(float(obs),4),"p":round(p,4)})
        if oN: print(f"  {'':<28}   NEG genes hit: {sorted(tgt&NEG)}")
json.dump(res,open(P/"outputs/apoptosis/wnt_direction.json","w"),indent=1)
print("\nwrote outputs/apoptosis/wnt_direction.json")

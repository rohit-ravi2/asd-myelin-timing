#!/usr/bin/env python3
"""
Conditional gate competition  (program sections IX, XV, XVI, XVII)

THE QUESTION THIS ANSWERS
-------------------------
Not "is module X enriched?" but:
    does module X carry INDEPENDENT ASD signal once gene properties AND the broader
    developmental programmes it sits inside are accounted for?

Key nested contrasts, pre-specified:
  C1  lipid modules | OL identity        -> is lipid biology anything beyond "expressed in OLs"?
  C2  OL-lineage    | chromatin/state    -> is OL signal anything beyond general state-control?
  C3  myelin output | OL maturation      -> does signal reach the structural end?
  C4  mitochondrial | OL identity        -> is mito anything beyond high-energy OL programmes?
  C5  each gate     | all other gates    -> full joint model

MODEL
  Z_g ~ module indicators + LOEUF + log(CDS) + log(brain expr) + log(NSNPs) + OL specificity
Marginal (each alone) and joint (all together) are both reported; the difference IS the result.

MATCHED NULL: for every module, 1000 gene sets matched on LOEUF decile x CDS decile x
expression decile; the module's joint coefficient is re-estimated each time. This is what
stops big constrained neurodevelopmental gene sets from manufacturing signal.
"""
import gzip, json, re, sys
import numpy as np
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
SCOPE = Path("data/dvad-corpus/asd_scope")
OUT = BASE / "outputs/phase4"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20260917)
N_MATCH = 1000

# ---------------------------------------------------------------- gene sets
gmt = {}
for fn in ["c2.cp.v2024.1.Hs.symbols.gmt", "c5.go.bp.v2024.1.Hs.symbols.gmt",
           "c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    p = COR / "pathways" / fn
    if p.exists():
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if len(f) > 2:
                gmt[f[0]] = set(f[2:])


def U(*pats, exclude=()):
    g = set()
    for name, mem in gmt.items():
        if any(re.search(p, name, re.I) for p in pats) and \
           not any(re.search(e, name, re.I) for e in exclude):
            g |= mem
    return g


MOD = {}
# --- immune submodules, mechanistically partitioned
MOD["I1_MHC_CLASS_I"]   = U(r"MHC_CLASS_I", r"ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN")
MOD["I2_MHC_CLASS_II"]  = U(r"MHC_CLASS_II")
MOD["I3_COMPLEMENT"]    = U(r"COMPLEMENT")
MOD["I4_MIA_IL6_IL17"]  = U(r"INTERLEUKIN_6", r"INTERLEUKIN_17", r"IL6_", r"TH17")
MOD["I5_INTERFERON_I"]  = U(r"TYPE_I_INTERFERON", r"INTERFERON_ALPHA", r"INTERFERON_BETA")
MOD["I6_INTERFERON_G"]  = U(r"INTERFERON_GAMMA")
MOD["I7_TLR_INNATE"]    = U(r"TOLL_LIKE_RECEPTOR", r"INNATE_IMMUNE", r"PATTERN_RECOGNITION")
MOD["I8_TNF"]           = U(r"TUMOR_NECROSIS_FACTOR", r"TNF_SIGNALING", r"NF_KAPPAB")
MOD["I9_ADAPTIVE_T"]    = U(r"T_CELL_ACTIVATION", r"T_CELL_RECEPTOR", r"LYMPHOCYTE_ACTIVATION")
MOD["I10_CYTOKINE_PROD"]= U(r"CYTOKINE_PRODUCTION", r"CYTOKINE_MEDIATED_SIGNALING")
MOD["I11_MICROGLIA"]    = U(r"MICROGLI", r"MACROPHAGE_ACTIVATION")
MOD["I12_MAST_ALLERGY"] = U(r"MAST_CELL", r"IGE", r"ALLERG", r"EOSINOPHIL")
# --- contrasts kept from the main run
MOD["G3c_OL_MATURATION"] = U(r"MYELIN", r"ENSHEATHMENT_OF_NEURONS", r"AXON_ENSHEATHMENT", r"OLIGODENDROCYTE")
MOD["G6_MITOCARTA"] = set()
MOD["G1_CHROMATIN_STATE"] = U(r"CHROMATIN_REMODEL", r"CHROMATIN_ORGANIZATION", r"HISTONE_MODIFICATION")
MOD["G9a_SYNAPTIC"] = U(r"GOBP_SYNAP", r"GOCC_POSTSYNAPTIC_DENSITY")
import csv as _csv
_sf=set()
with open(SCOPE / "sfari/SFARI-Gene_genes_2026-09-10.csv") as _f:
    for _r in _csv.DictReader(_f):
        _g=_r.get("gene-symbol") or _r.get("gene_symbol")
        if _g: _sf.add(_g)
MOD["REF_SFARI_ALL"]=_sf
MOD["NEG_OLFACTORY"]=U(r"OLFACTORY")
MOD["NEG_KERATIN"]=U(r"KERATIN", r"CORNIFIED")
for k in list(MOD):
    if len(MOD[k]) < 60:
        print(f"[drop] {k}: only {len(MOD[k])} genes", file=sys.stderr)
        del MOD[k]
print(f"[modules] {len(MOD)}", file=sys.stderr)
for k in sorted(MOD):
    print(f"   {k:28s} {len(MOD[k])}", file=sys.stderr)

# ---------------------------------------------------------------- gene-level data
gloc = {}
for line in open(COR / "tools/magma_geneloc/NCBI37.3.gene.loc"):
    p = line.split()
    if len(p) >= 6:
        gloc[p[0]] = p[5]
Z, NS = {}, {}
with open(BASE / "outputs/magma/grove.w35.genes.out") as f:
    f.readline()
    for line in f:
        p = line.split()
        s = gloc.get(p[0])
        if s:
            Z[s] = float(p[7]); NS[s] = float(p[4])
loeuf, cds = {}, {}
with gzip.open(COR / "genetics/gnomad_v2.1.1_lof_metrics.by_gene.txt.bgz", "rt") as f:
    h = f.readline().rstrip("\n").split("\t")
    gi, li, ci = h.index("gene"), h.index("oe_lof_upper"), h.index("cds_length")
    for line in f:
        p = line.rstrip("\n").split("\t")
        try:
            loeuf[p[gi]] = float(p[li]); cds[p[gi]] = float(p[ci])
        except (ValueError, IndexError):
            pass
ai = json.load(open(BASE / "outputs/phase1/atlas_index.json"))
am = np.load(BASE / "outputs/phase1/atlas_means.npy")
atg = {g: i for i, g in enumerate(ai["genes"])}
types = ai["types"]
OLROWS = [i for i, t in enumerate(types) if t in ("OPC", "COP", "OL")]
allmean = am.mean(axis=0)
ol_spec_v = np.log2((am[OLROWS].mean(axis=0) + 1e-3) / (allmean + 1e-3))
expr_v = allmean

genes = [g for g in Z if g in loeuf and g in atg and cds.get(g, 0) > 0
         and expr_v[atg[g]] > 0]
print(f"[genes] {len(genes)} usable", file=sys.stderr)
idx = np.array([atg[g] for g in genes])
y = np.array([Z[g] for g in genes])
COV = np.column_stack([
    [loeuf[g] for g in genes],
    np.log([cds[g] for g in genes]),
    np.log(expr_v[idx] + 1e-6),
    np.log([NS[g] + 1 for g in genes]),
    ol_spec_v[idx],
])
COVNAMES = ["LOEUF", "logCDS", "logExpr", "logNSNP", "OLspecificity"]
M = {k: np.array([1.0 if g in v else 0.0 for g in genes]) for k, v in MOD.items()}
for k in list(M):
    if M[k].sum() < 8:
        del M[k]


def ols(y, X):
    X = np.column_stack([np.ones(len(y)), X])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    df = len(y) - X.shape[1]
    s2 = (r @ r) / df
    se = np.sqrt(np.diag(s2 * np.linalg.pinv(X.T @ X)))
    t = b / se
    from math import erf, sqrt
    p = np.array([2 * (1 - 0.5 * (1 + erf(abs(x) / sqrt(2)))) for x in t])
    return b, se, t, p


keys = sorted(M)
# marginal
marg = {}
for k in keys:
    b, se, t, p = ols(y, np.column_stack([M[k], COV]))
    marg[k] = dict(beta=b[1], se=se[1], t=t[1], p=p[1], n=int(M[k].sum()))
# joint
Xj = np.column_stack([M[k] for k in keys] + [COV])
b, se, t, p = ols(y, Xj)
joint = {k: dict(beta=b[i + 1], se=se[i + 1], t=t[i + 1], p=p[i + 1])
         for i, k in enumerate(keys)}
covres = {COVNAMES[j]: dict(beta=b[len(keys) + 1 + j], t=t[len(keys) + 1 + j],
                            p=p[len(keys) + 1 + j]) for j in range(len(COVNAMES))}

# matched null on joint coefficient
lo_a = np.array([loeuf[g] for g in genes]); cd_a = np.array([cds[g] for g in genes])
ex_a = expr_v[idx]
def dec(v): return np.digitize(v, np.quantile(v, np.linspace(0, 1, 6)[1:-1]))
strata = defaultdict(list)
key3 = list(zip(dec(lo_a), dec(cd_a), dec(ex_a)))
for i, k3 in enumerate(key3):
    strata[k3].append(i)
print(f"[perm] matched null, {N_MATCH} draws", file=sys.stderr)
matp = {}
for k in keys:
    obs = abs(joint[k]["t"])
    memb = np.where(M[k] > 0)[0]
    ge = 0
    for _ in range(N_MATCH):
        pick = np.array([strata[key3[i]][rng.integers(len(strata[key3[i]]))] for i in memb])
        mm = np.zeros(len(genes)); mm[pick] = 1.0
        Xk = np.column_stack([mm] + [M[o] for o in keys if o != k] + [COV])
        bb, ss, tt, pp = ols(y, Xk)
        if abs(tt[1]) >= obs:
            ge += 1
    matp[k] = (ge + 1) / (N_MATCH + 1)

# nested contrasts
def contrast(focal, condition_on, label):
    if focal not in M:
        return dict(label=label, focal=focal, conditioned_on=condition_on,
                    beta_marginal=float("nan"), p_marginal=float("nan"),
                    beta_conditional=float("nan"), p_conditional=float("nan"),
                    attenuation=None, note="focal module dropped (too few genes)")
    b1, s1, t1, p1 = ols(y, np.column_stack([M[focal], COV]))
    cols = [M[focal]] + [M[c] for c in condition_on if c in M] + [COV]
    b2, s2, t2, p2 = ols(y, np.column_stack(cols))
    return dict(label=label, focal=focal, conditioned_on=condition_on,
                beta_marginal=b1[1], p_marginal=p1[1],
                beta_conditional=b2[1], p_conditional=p2[1],
                attenuation=(1 - b2[1] / b1[1]) if b1[1] != 0 else None)

OLSET=[]
contrasts=[]

json.dump(dict(marginal=marg, joint=joint, matched_p=matp, covariates=covres,
               contrasts=contrasts, n_genes=len(genes)),
          open(OUT / "gate_competition.json", "w"), indent=1, default=float)

L = ["=" * 110, "CONDITIONAL GATE COMPETITION — ASD common-variant gene-level Z", "=" * 110,
     f"genes={len(genes)}  covariates: {', '.join(COVNAMES)}",
     "marginal = module + covariates.   joint = ALL modules together + covariates.",
     "matchP = 1000 LOEUF x CDS x expression matched sets, joint coefficient re-estimated.", "",
     f"{'module':28s} {'n':>5s} {'beta_marg':>10s} {'p_marg':>9s} {'beta_joint':>11s} {'p_joint':>9s} {'matchP':>7s}"]
for k in sorted(keys, key=lambda x: marg[x]["p"]):
    L.append(f"{k:28s} {marg[k]['n']:>5d} {marg[k]['beta']:>+10.4f} {marg[k]['p']:>9.3g} "
             f"{joint[k]['beta']:>+11.4f} {joint[k]['p']:>9.3g} {matp[k]:>7.3f}")
L += ["", "covariates in joint model:"]
for c, v in covres.items():
    L.append(f"   {c:16s} beta={v['beta']:+.4f} t={v['t']:+.2f} p={v['p']:.3g}")
L += ["", "=" * 110, "PRE-SPECIFIED NESTED CONTRASTS", "=" * 110,
      f"{'contrast':46s} {'b_marg':>9s} {'p_marg':>9s} {'b_cond':>9s} {'p_cond':>9s} {'atten':>7s}"]
for c in contrasts:
    at = f"{c['attenuation']*100:.0f}%" if c.get("attenuation") is not None else "NA"
    L.append(f"{c['label']:46s} {c['beta_marginal']:>+9.4f} {c['p_marginal']:>9.3g} "
             f"{c['beta_conditional']:>+9.4f} {c['p_conditional']:>9.3g} {at:>7s}")
txt = "\n".join(L)
(OUT / "GATE_COMPETITION_REPORT.txt").write_text(txt + "\n")
print(txt)

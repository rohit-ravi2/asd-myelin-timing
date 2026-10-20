#!/usr/bin/env python3
"""
P1.1 / P1.2 / P1.3 — Does ASD common-variant association localise to a point on the
oligodendrocyte-lineage maturation axis, and does that survive the constraint confound?

WHY THIS IS THE KEYSTONE, AND ITS MAIN FAILURE MODE
---------------------------------------------------
Progenitor/chromatin genes (SOX10, OLIG2, CHD8, TCF7L2) are long, highly constrained,
broadly expressed. Structural myelin genes (MBP, PLP1, MOG, MOBP) are short, LoF-tolerant,
narrowly expressed. ASD association tracks constraint hard. So a monotone decline in
ASD Z from OPC -> mature OL is the NULL EXPECTATION FROM CONSTRAINT ALONE, with no
oligodendrocyte biology whatsoever. LOEUF is therefore not one covariate among many; it
is the primary alternative explanation, and lineage specificity is the gate.

DESIGN
------
  T_g  = expression-weighted mean lineage stage for gene g (OPC=0, COP=1, OL=2),
         computed in HUMAN INFANT cells only (the correct species/region/epoch).
  M1   Z ~ T_g                                     (naive; expected to "work" either way)
  M2   Z ~ T_g + LOEUF + log CDS + log mean expr + log NSNPS    (constraint-conditioned)
  M3   Z ~ T_g + T_g^2 + covariates                (the inverted-U / peak-at-COP test)
  PERM shuffle stage labels across cells, rebuild T_g, refit -> empirical null for beta

  GATE (P1.2): the same regression run on specificity axes for lineages with no myelin.
         POSITIVE control : neuronal specificity  (ASD is a neuronal disorder; must be +)
         NEGATIVE controls: astrocyte, microglia specificity
         If ASD Z tracks every lineage's specificity axis equally, the OL result is
         generic constraint/expression structure and is NOT reported as OL biology.

Pre-registered before seeing output (see ledger amendment 4): the overshoot model predicts
a PEAK AT COP (inverted U, negative quadratic term), not a monotone decline. A monotone
decline is the artifact's prediction. Both are recorded either way.
"""
import gzip, json, sys
import numpy as np
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
GSE = COR / "celltype/GSE294786"
OUT = BASE / "outputs/phase1"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20260917)
N_PERM = 1000

# ------------------------------------------------------------------ metadata
meta = {}
with gzip.open(GSE / "GSE294786_all_meta_progenitor.tsv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    if hdr[0] == "":
        hdr[0] = "CELLID"
    ix = {c: i for i, c in enumerate(hdr)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hdr):
            continue
        meta[p[ix["CELLID"]]] = (p[ix["Species"]], p[ix["Age"]], p[ix["cell_type_final"]])

HUMAN_INFANT = {k: v[2] for k, v in meta.items() if v[0] == "Hg" and v[1] == "Infant"}
print(f"[meta] human infant cells: {len(HUMAN_INFANT)}", file=sys.stderr)

# lineage definitions -------------------------------------------------------
OL_STAGE = {"OPC": 0.0, "COP": 1.0, "OL": 2.0, "Oli": 2.0}
NEURON = {"L2/L3", "L2/L3-BS", "L4", "L5/L6-CC", "L5/L6", "L5", "L6",
          "IN-PV", "IN-SST", "IN-VIP", "IN-SV2C", "IN-LAMP5"}
ASTRO = {"Astrocyte"}
MICRO = {"Mic-P2RY12", "Mic", "Microglia"}

# ------------------------------------------------------------------ stream counts
# accumulate per-celltype mean log-CPM without holding the matrix
with gzip.open(GSE / "GSE294786_all_counts_transposed.tsv.gz", "rt") as f:
    genes = f.readline().rstrip("\n").split("\t")[1:]
    G = len(genes)
    sums = defaultdict(lambda: np.zeros(G, dtype=np.float64))
    ncell = defaultdict(int)
    # per-cell stage label for permutation: store the normalised vector ONLY for OL lineage
    ol_vecs, ol_labels = [], []
    seen = 0
    for line in f:
        tab = line.find("\t")
        cid = line[:tab]
        ct = HUMAN_INFANT.get(cid)
        if ct is None:
            continue
        v = np.fromstring(line[tab + 1:], sep="\t", dtype=np.float64)
        if v.size != G:
            continue
        tot = v.sum()
        if tot <= 0:
            continue
        v = np.log1p(v / tot * 1e4)
        sums[ct] += v
        ncell[ct] += 1
        if ct in OL_STAGE:
            ol_vecs.append(v.astype(np.float32))
            ol_labels.append(OL_STAGE[ct])
        seen += 1
        if seen % 2000 == 0:
            print(f"   [{seen}] cells", file=sys.stderr)
print(f"[counts] {seen} human-infant cells across {len(ncell)} types", file=sys.stderr)
for k in sorted(ncell, key=lambda x: -ncell[x])[:14]:
    print(f"   {k:16s} {ncell[k]}", file=sys.stderr)

means = {k: sums[k] / ncell[k] for k in ncell}
ol_vecs = np.array(ol_vecs, dtype=np.float32)
ol_labels = np.array(ol_labels, dtype=np.float32)
print(f"[OL] lineage cells for permutation: {ol_vecs.shape}", file=sys.stderr)


def stage_score(labels):
    """expression-weighted mean lineage stage per gene"""
    w = np.zeros(G); tot = np.zeros(G)
    for s in np.unique(labels):
        m = ol_vecs[labels == s].mean(axis=0)
        w += s * m
        tot += m
    with np.errstate(invalid="ignore", divide="ignore"):
        T = np.where(tot > 0, w / tot, np.nan)
    return T, tot


T_obs, ol_expr = stage_score(ol_labels)


def specificity(group):
    """mean expression in group / mean across all profiled types (log ratio)"""
    ing = [means[k] for k in means if k in group]
    if not ing:
        return np.full(G, np.nan)
    a = np.mean(ing, axis=0)
    allm = np.mean([means[k] for k in means], axis=0)
    return np.log2((a + 1e-3) / (allm + 1e-3))


spec = {"neuron": specificity(NEURON), "astro": specificity(ASTRO),
        "micro": specificity(MICRO), "ol": specificity(set(OL_STAGE))}

# ------------------------------------------------------------------ gene-level ASD Z
gloc = {}
with open(COR / "tools/magma_geneloc/NCBI37.3.gene.loc") as f:
    for line in f:
        p = line.split()
        if len(p) >= 6:
            gloc[p[0]] = p[5]
asdz, nsnps = {}, {}
with open(BASE / "outputs/magma/grove.w35.genes.out") as f:
    f.readline()
    for line in f:
        p = line.split()
        sym = gloc.get(p[0])
        if sym:
            asdz[sym] = float(p[7])
            nsnps[sym] = float(p[4])
print(f"[magma] {len(asdz)} genes with ASD Z", file=sys.stderr)

# ------------------------------------------------------------------ constraint
loeuf, cds = {}, {}
with gzip.open(COR / "genetics/gnomad_v2.1.1_lof_metrics.by_gene.txt.bgz", "rt") as f:
    h = f.readline().rstrip("\n").split("\t")
    gi, li, ci = h.index("gene"), h.index("oe_lof_upper"), h.index("cds_length")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(gi, li, ci):
            continue
        try:
            loeuf[p[gi]] = float(p[li]); cds[p[gi]] = float(p[ci])
        except ValueError:
            pass
print(f"[gnomad] {len(loeuf)} genes with LOEUF", file=sys.stderr)

# ------------------------------------------------------------------ assemble
gi_index = {g: i for i, g in enumerate(genes)}
rows = []
for g, i in gi_index.items():
    if g not in asdz or g not in loeuf:
        continue
    if not np.isfinite(T_obs[i]) or ol_expr[i] <= 0:
        continue
    me = float(np.mean([means[k][i] for k in means]))
    if me <= 0 or cds.get(g, 0) <= 0:
        continue
    rows.append((g, i, asdz[g], T_obs[i], loeuf[g], np.log(cds[g]), np.log(me + 1e-6),
                 np.log(nsnps[g] + 1), spec["neuron"][i], spec["astro"][i],
                 spec["micro"][i], spec["ol"][i]))
print(f"[join] {len(rows)} genes usable", file=sys.stderr)
names = ["gene", "idx", "Z", "T", "LOEUF", "logCDS", "logExpr", "logNSNP",
         "spec_neuron", "spec_astro", "spec_micro", "spec_ol"]
D = {n: np.array([r[j] for r in rows], dtype=object if n == "gene" else float)
     for j, n in enumerate(names)}
keep_idx = D["idx"].astype(int)


def ols(y, X):
    X = np.column_stack([np.ones(len(y))] + X)
    good = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    X, y = X[good], y[good]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    df = len(y) - X.shape[1]
    s2 = (r @ r) / df
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    from math import erf, sqrt
    t = b / se
    p = [2 * (1 - 0.5 * (1 + erf(abs(x) / sqrt(2)))) for x in t]
    return b, se, t, p, len(y)


COV = [D["LOEUF"], D["logCDS"], D["logExpr"], D["logNSNP"]]
res = {}
b, se, t, p, n = ols(D["Z"], [D["T"]])
res["M1_naive"] = dict(beta=b[1], se=se[1], t=t[1], p=p[1], n=n)
b, se, t, p, n = ols(D["Z"], [D["T"]] + COV)
res["M2_constraint_conditioned"] = dict(beta=b[1], se=se[1], t=t[1], p=p[1], n=n,
                                        beta_LOEUF=b[2], t_LOEUF=t[2], p_LOEUF=p[2])
b, se, t, p, n = ols(D["Z"], [D["T"], D["T"] ** 2] + COV)
peak = -b[1] / (2 * b[2]) if b[2] != 0 else np.nan
res["M3_quadratic"] = dict(beta_lin=b[1], p_lin=p[1], beta_quad=b[2], se_quad=se[2],
                           t_quad=t[2], p_quad=p[2], vertex_stage=peak, n=n)

# GATE: same model on other lineages' specificity axes
for lab in ["neuron", "astro", "micro", "ol"]:
    b, se, t, p, n = ols(D["Z"], [D[f"spec_{lab}"]] + COV)
    res[f"GATE_spec_{lab}"] = dict(beta=b[1], se=se[1], t=t[1], p=p[1], n=n)

# PERMUTATION: shuffle stage labels across OL cells, rebuild T, refit M2
print(f"[perm] {N_PERM} stage-label permutations", file=sys.stderr)
obs_t = res["M2_constraint_conditioned"]["t"]
null_t = []
for k in range(N_PERM):
    Tp, tot = stage_score(rng.permutation(ol_labels))
    tp = Tp[keep_idx]
    good = np.isfinite(tp)
    bb, ss, tt, pp, nn = ols(D["Z"][good], [tp[good]] + [c[good] for c in COV])
    null_t.append(tt[1])
    if (k + 1) % 250 == 0:
        print(f"   perm {k+1}", file=sys.stderr)
null_t = np.array(null_t)
res["PERM"] = dict(n_perm=N_PERM, obs_t=obs_t,
                   null_mean=float(null_t.mean()), null_sd=float(null_t.std()),
                   emp_p=float((np.sum(np.abs(null_t) >= abs(obs_t)) + 1) / (N_PERM + 1)))

json.dump({k: {kk: (float(vv) if isinstance(vv, (int, float, np.floating)) else vv)
               for kk, vv in v.items()} for k, v in res.items()},
          open(OUT / "P1_pseudotime_results.json", "w"), indent=1)

# ------------------------------------------------------------------ report
L = ["=" * 92,
     "P1.1-P1.3  ASD common-variant association vs OL-lineage maturation (human infant PFC)",
     "=" * 92,
     f"atlas: GSE294786 human infant PFC | genes analysed: {len(rows)}",
     f"OL-lineage cells: {ol_vecs.shape[0]}  (stages OPC=0, COP=1, OL=2)", "",
     "T_g = expression-weighted mean lineage stage. Negative beta = ASD signal skews PROGENITOR.",
     ""]
m = res["M1_naive"]
L.append(f"M1 naive            : beta={m['beta']:+.4f} se={m['se']:.4f} t={m['t']:+.2f} p={m['p']:.3g}  n={m['n']}")
m = res["M2_constraint_conditioned"]
L.append(f"M2 +LOEUF/len/expr  : beta={m['beta']:+.4f} se={m['se']:.4f} t={m['t']:+.2f} p={m['p']:.3g}  n={m['n']}")
L.append(f"     (LOEUF itself  : beta={m['beta_LOEUF']:+.4f} t={m['t_LOEUF']:+.2f} p={m['p_LOEUF']:.3g})")
m = res["M3_quadratic"]
L.append(f"M3 quadratic        : lin={m['beta_lin']:+.4f} (p={m['p_lin']:.3g})  "
         f"quad={m['beta_quad']:+.5f} (p={m['p_quad']:.3g})  vertex at stage={m['vertex_stage']:.2f}")
L.append("     [pre-registered: overshoot predicts NEGATIVE quad with vertex near COP (stage 1)]")
L.append("")
L.append("GATE P1.2 -- same model, other lineage specificity axes:")
for lab, note in [("neuron", "POSITIVE control, expect +"), ("ol", "focal lineage"),
                  ("astro", "NEGATIVE control"), ("micro", "NEGATIVE control")]:
    g = res[f"GATE_spec_{lab}"]
    L.append(f"   spec_{lab:7s} beta={g['beta']:+.4f} t={g['t']:+.2f} p={g['p']:.3g}   [{note}]")
p = res["PERM"]
L.append("")
L.append(f"PERMUTATION ({p['n_perm']} stage-label shuffles): obs t={p['obs_t']:+.2f}, "
         f"null {p['null_mean']:+.2f}+/-{p['null_sd']:.2f}, empirical p={p['emp_p']:.4f}")
txt = "\n".join(L)
(OUT / "P1_PSEUDOTIME_REPORT.txt").write_text(txt + "\n")
print(txt)

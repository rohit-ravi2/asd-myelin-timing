#!/usr/bin/env python3
"""
§8-§16 — INDEPENDENT OL-ATLAS REPLICATION of the ASD developmental-state test.

Atlas 2: GSE155488 (Perlman/Fuzik-type human OL-lineage scRNA-seq, PMID 31958186).
  second-trimester fetal (n=4 pooled, 360 cells) | 2 yr x2 (3,438) | 13 yr (3,087)
  | 58 yr + 62 yr (3,512).  Total 10,397 OL-LINEAGE cells.
  This atlas has the MATURE END that GSE294786 lacks (66 mature OL there).
  CAVEAT (§8): surgical tissue, so donors are not neurologically normal.

§10 requires several maturation scores, and CONCORDANCE matters more than one p-value:
  A  pseudotime-correlation score : corr(gene, cell maturation score)
  B  stage-centroid score         : expression-weighted mean stage
  C  mature-vs-OPC contrast       : log2(MOL/OPC)
  D  nonlinear peak stage         : which stage bin maximises ASD Z

§6/§11/§12 gates are mandatory: covariate conditioning, matched-null permutation, and
matched negative-control lineages are NOT available inside an OL-only atlas, so the
lineage gate is carried over from atlas 1 and this run tests replication of the OL result
itself.
"""
import gzip, json, sys
import numpy as np
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
GSE = COR / "celltype/GSE155488"
OUT = BASE / "outputs/phase6"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20260917)
N_PERM = 1000

SAMPLES = [
    ("GSE155488_All_fetal_OL_lineage_subset.csv.gz", "fetal", 0),
    ("GSM4875773_HA_773_OL_lineage_subset.csv.gz", "2yr", 2),
    ("GSM4875774_HA_776_OL_lineage_subset.csv.gz", "2yr", 2),
    ("GSM4875775_HA_772_OL_lineage_subset.csv.gz", "13yr", 13),
    ("GSM4875776_HA_738_OL_lineage_subset.csv.gz", "adult", 62),
    ("GSM4875777_OPK_363_OL_lineage_subset.csv.gz", "adult", 58),
]

# canonical human OL-lineage stage markers
STAGE = {
    "OPC":  "PDGFRA CSPG4 PTPRZ1 LHFPL3 SOX6 OLIG1 DSCAM GPR17".split(),
    "COP":  "GPR17 BCAS1 ENPP6 TMEM235 NKX2-2 SIRT2 FYN".split(),
    "NFOL": "BCAS1 ENPP6 ITPR2 CNP TCF7L2 NKX6-2".split(),
    "MFOL": "MAL MOG OPALIN CTPS1 CNP PLP1".split(),
    "MOL":  "PLP1 MBP MOBP CLDN11 APOD KLK6 ASPA CARNS1 SLC5A11".split(),
}
STAGE_IDX = {"OPC": 0.0, "COP": 1.0, "NFOL": 2.0, "MFOL": 3.0, "MOL": 4.0}

# ------------------------------------------------------------------ load
cells, gene_union = [], None
mats = []
for fn, agelab, age in SAMPLES:
    p = GSE / fn
    if not p.exists():
        print(f"[skip] {fn}", file=sys.stderr)
        continue
    with gzip.open(p, "rt") as f:
        rows, gnames = [], []
        for line in f:
            i = line.find(",")
            gnames.append(line[:i].strip().strip('"'))
            rows.append(np.fromstring(line[i + 1:], sep=",", dtype=np.float32))
    M = np.vstack(rows)                    # genes x cells
    mats.append((M, gnames, agelab, age))
    print(f"[load] {fn}: {M.shape[0]} genes x {M.shape[1]} cells ({agelab})", file=sys.stderr)

common = set(mats[0][1])
for M, g, _, _ in mats[1:]:
    common &= set(g)
common = sorted(common)
ci = {g: i for i, g in enumerate(common)}
print(f"[genes] {len(common)} shared across samples", file=sys.stderr)

blocks, ages = [], []
for M, g, agelab, age in mats:
    idx = np.array([g.index(x) for x in common])
    sub = M[idx]
    tot = sub.sum(axis=0); tot[tot <= 0] = 1
    blocks.append(np.log1p(sub / tot * 1e4).astype(np.float32))
    ages += [age] * sub.shape[1]
X = np.hstack(blocks)                       # genes x cells, log-CPM
ages = np.array(ages)
print(f"[matrix] {X.shape[0]} genes x {X.shape[1]} cells", file=sys.stderr)

# ------------------------------------------------------------------ cell maturation score
def score_cells(markers):
    ix = [ci[m] for m in markers if m in ci]
    return X[ix].mean(axis=0) if ix else np.zeros(X.shape[1])

sc = {s: score_cells(m) for s, m in STAGE.items()}
for s in sc:                                 # z-score each stage score across cells
    v = sc[s]
    sc[s] = (v - v.mean()) / (v.std() + 1e-9)
S = np.vstack([sc[s] for s in ["OPC", "COP", "NFOL", "MFOL", "MOL"]])
assign = np.argmax(S, axis=0)
w = np.exp(S - S.max(axis=0))
w /= w.sum(axis=0)
cell_mat = (w * np.array([[0.], [1.], [2.], [3.], [4.]])).sum(axis=0)  # soft stage 0-4

print("[stages] cell counts by argmax stage:", file=sys.stderr)
for k, s in enumerate(["OPC", "COP", "NFOL", "MFOL", "MOL"]):
    print(f"   {s:5s} {int((assign==k).sum()):6d}", file=sys.stderr)
print("[stages] by donor age:", file=sys.stderr)
for a in sorted(set(ages)):
    m = ages == a
    print(f"   age {a:>3}: n={m.sum():5d}  mean stage={cell_mat[m].mean():.2f}", file=sys.stderr)

# ------------------------------------------------------------------ gene-level scores
expr = X.mean(axis=1)
keep = expr > 0.01
# A: correlation of gene expression with cell maturation score
cm = (cell_mat - cell_mat.mean()) / (cell_mat.std() + 1e-9)
Xc = X - X.mean(axis=1, keepdims=True)
sdg = X.std(axis=1) + 1e-9
T_corr = (Xc @ cm) / (X.shape[1] * sdg)
# B: expression-weighted mean stage (centroid)
stage_mean = np.vstack([X[:, assign == k].mean(axis=1) if (assign == k).any()
                        else np.zeros(X.shape[0]) for k in range(5)])
tot = stage_mean.sum(axis=0) + 1e-9
T_cent = (stage_mean * np.arange(5)[:, None]).sum(axis=0) / tot
# C: mature-vs-OPC contrast
mol = stage_mean[4] + 1e-3; opc = stage_mean[0] + 1e-3
T_contrast = np.log2(mol / opc)

# ------------------------------------------------------------------ ASD gene Z + covariates
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
    a, b, c = h.index("gene"), h.index("oe_lof_upper"), h.index("cds_length")
    for line in f:
        p = line.rstrip("\n").split("\t")
        try:
            loeuf[p[a]] = float(p[b]); cds[p[a]] = float(p[c])
        except (ValueError, IndexError):
            pass

use = [g for g in common if g in Z and g in loeuf and cds.get(g, 0) > 0
       and keep[ci[g]] and np.isfinite(T_cent[ci[g]])]
gi = np.array([ci[g] for g in use])
y = np.array([Z[g] for g in use])
COV = np.column_stack([[loeuf[g] for g in use], np.log([cds[g] for g in use]),
                       np.log(expr[gi] + 1e-6), np.log([NS[g] + 1 for g in use])])
print(f"[join] {len(use)} genes usable", file=sys.stderr)


def ols(y, X):
    X = np.column_stack([np.ones(len(y))] + list(X))
    ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    X, y = X[ok], y[ok]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b; df = len(y) - X.shape[1]
    se = np.sqrt(np.diag((r @ r) / df * np.linalg.pinv(X.T @ X)))
    t = b / se
    from math import erf, sqrt
    p = [2 * (1 - 0.5 * (1 + erf(abs(x) / sqrt(2)))) for x in t]
    return b, se, t, p, len(y)


res = {}
for nm, vec in [("A_pseudotime_corr", T_corr), ("B_stage_centroid", T_cent),
                ("C_mature_vs_OPC", T_contrast)]:
    v = vec[gi]
    b, se, t, p, n = ols(y, [v])
    res[nm + "_naive"] = dict(beta=b[1], t=t[1], p=p[1], n=n)
    b, se, t, p, n = ols(y, [v] + [COV[:, j] for j in range(COV.shape[1])])
    res[nm + "_conditioned"] = dict(beta=b[1], se=se[1], t=t[1], p=p[1], n=n)

# D: nonlinear — ASD Z by centroid-stage bin
bins = np.digitize(T_cent[gi], [0.5, 1.5, 2.5, 3.5])
resid = y - np.column_stack([np.ones(len(y))] + [COV[:, j] for j in range(COV.shape[1])]) @ \
    np.linalg.lstsq(np.column_stack([np.ones(len(y))] + [COV[:, j] for j in range(COV.shape[1])]),
                    y, rcond=None)[0]
res["D_stage_bins"] = {["OPC", "COP", "NFOL", "MFOL", "MOL"][k]:
                       dict(n=int((bins == k).sum()),
                            mean_resid_Z=float(resid[bins == k].mean()) if (bins == k).any() else None)
                       for k in range(5)}
# quadratic
b, se, t, p, n = ols(y, [T_cent[gi], T_cent[gi] ** 2] + [COV[:, j] for j in range(COV.shape[1])])
res["D_quadratic"] = dict(beta_lin=b[1], p_lin=p[1], beta_quad=b[2], p_quad=p[2],
                          vertex=float(-b[1] / (2 * b[2])) if b[2] else None)

# permutation on the centroid score (shuffle cell stage assignment)
obs_t = res["B_stage_centroid_conditioned"]["t"]
null = []
for _ in range(N_PERM):
    pa = rng.permutation(assign)
    sm = np.vstack([X[:, pa == k].mean(axis=1) if (pa == k).any() else np.zeros(X.shape[0])
                    for k in range(5)])
    tt_ = (sm * np.arange(5)[:, None]).sum(axis=0) / (sm.sum(axis=0) + 1e-9)
    bb, ss, tt, pp, nn = ols(y, [tt_[gi]] + [COV[:, j] for j in range(COV.shape[1])])
    null.append(tt[1])
null = np.array(null)
res["PERM"] = dict(n_perm=N_PERM, obs_t=obs_t, null_mean=float(null.mean()),
                   null_sd=float(null.std()),
                   emp_p=float((np.sum(np.abs(null) >= abs(obs_t)) + 1) / (N_PERM + 1)))

json.dump(res, open(OUT / "atlas2_replication.json", "w"), indent=1, default=float)

L = ["=" * 100, "§8-16 INDEPENDENT OL-ATLAS REPLICATION — GSE155488", "=" * 100,
     f"{X.shape[1]} OL-lineage cells, {len(use)} genes with ASD Z + covariates",
     "stage assignment (argmax of marker score):"]
for k, s in enumerate(["OPC", "COP", "NFOL", "MFOL", "MOL"]):
    L.append(f"   {s:5s} {int((assign==k).sum()):6d} cells")
L.append("")
L.append(f"{'maturation score':26s} {'beta_naive':>11s} {'p_naive':>10s} {'beta_cond':>11s} {'p_cond':>10s}")
for nm in ["A_pseudotime_corr", "B_stage_centroid", "C_mature_vs_OPC"]:
    a_, c_ = res[nm + "_naive"], res[nm + "_conditioned"]
    L.append(f"{nm:26s} {a_['beta']:>+11.4f} {a_['p']:>10.3g} {c_['beta']:>+11.4f} {c_['p']:>10.3g}")
L.append("")
L.append("D. ASD Z (covariate-residualised) by stage bin:")
for s, d in res["D_stage_bins"].items():
    if d["mean_resid_Z"] is not None:
        L.append(f"   {s:5s} n={d['n']:>6d}  mean resid Z = {d['mean_resid_Z']:+.4f}")
q = res["D_quadratic"]
L.append(f"   quadratic: lin={q['beta_lin']:+.4f} (p={q['p_lin']:.3g}) "
         f"quad={q['beta_quad']:+.5f} (p={q['p_quad']:.3g}) vertex={q['vertex']}")
pm = res["PERM"]
L.append("")
L.append(f"PERMUTATION ({pm['n_perm']}): obs t={pm['obs_t']:+.2f}, null {pm['null_mean']:+.2f}"
         f"+/-{pm['null_sd']:.2f}, empirical p={pm['emp_p']:.4f}")
L.append("")
L.append("CAVEAT (§8): surgical tissue; donors not neurologically normal. Fetal n=360 pooled.")
txt = "\n".join(L)
(OUT / "ATLAS2_REPLICATION_REPORT.txt").write_text(txt + "\n")
print(txt)

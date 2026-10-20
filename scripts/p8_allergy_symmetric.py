#!/usr/bin/env python3
"""
SYMMETRIC TEST — does IMMUNE liability localise to the OL lineage where ASD does not?

The file formerly mislabelled ASD_Grove2019.gz is a large, well-powered allergic-disease
GWAS (N=242,569, 81 GWS loci: TSLP, IL33, IL1RL1/IL18R1, ORMDL3/GSDMB, SMAD3, TLR1,
EMSY/LRRC32, MHC-II; h2 Z=20.6). That makes it an ideal IMMUNE-TRAIT comparison trait.

Three outcomes, each informative:
  (a) allergy localises to OL but ASD does not  -> the immune->OL genetic bridge is real
      and ASD simply enters the system elsewhere.
  (b) neither localises to OL                   -> no genetic form of the immune->OL bridge
      exists in either trait; any immune->OL effect must be non-genetic / exposure-driven.
  (c) both localise                             -> shared route.

Identical model to P1.2 (atlas 1 cell-type specificity, gene-property conditioned), run on
both traits so the ONLY difference is the GWAS.
"""
import gzip, json, sys
import numpy as np
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
P1 = BASE / "outputs/phase1"
OUT = BASE / "outputs/phase8"
OUT.mkdir(parents=True, exist_ok=True)

ai = json.load(open(P1 / "atlas_index.json"))
am = np.load(P1 / "atlas_means.npy")
types, genes = ai["types"], ai["genes"]
ti = {t: i for i, t in enumerate(types)}
gi = {g: i for i, g in enumerate(genes)}
allmean = am.mean(axis=0)

NEUR = [t for t in types if t.startswith(("L2", "L4", "L5", "L6", "IN-"))]
GROUPS = {"ol": [t for t in types if t in ("OPC", "COP", "OL")],
          "opc": [t for t in types if t == "OPC"],
          "cop": [t for t in types if t == "COP"],
          "neuron": NEUR,
          "astro": [t for t in types if t == "Astrocyte"],
          "micro": [t for t in types if t.startswith("Mic")]}
spec = {}
for k, ts in GROUPS.items():
    if not ts:
        continue
    a = np.mean([am[ti[t]] for t in ts], axis=0)
    spec[k] = np.log2((a + 1e-3) / (allmean + 1e-3))

gloc = {}
for line in open(COR / "tools/magma_geneloc/NCBI37.3.gene.loc"):
    p = line.split()
    if len(p) >= 6:
        gloc[p[0]] = p[5]
loeuf, cds = {}, {}
with gzip.open(COR / "genetics/gnomad_v2.1.1_lof_metrics.by_gene.txt.bgz", "rt") as f:
    h = f.readline().rstrip("\n").split("\t")
    a_, b_, c_ = h.index("gene"), h.index("oe_lof_upper"), h.index("cds_length")
    for line in f:
        p = line.rstrip("\n").split("\t")
        try:
            loeuf[p[a_]] = float(p[b_]); cds[p[a_]] = float(p[c_])
        except (ValueError, IndexError):
            pass


def load_z(path):
    Z, NS = {}, {}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split()
            s = gloc.get(p[0])
            if s:
                Z[s] = float(p[7]); NS[s] = float(p[4])
    return Z, NS


def ols(y, cols):
    X = np.column_stack([np.ones(len(y))] + list(cols))
    ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    X, y = X[ok], y[ok]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b; df = len(y) - X.shape[1]
    se = np.sqrt(np.diag((r @ r) / df * np.linalg.pinv(X.T @ X)))
    t = b / se
    from math import erf, sqrt
    p = [2 * (1 - 0.5 * (1 + erf(abs(x) / sqrt(2)))) for x in t]
    return b, se, t, p, len(y)


TRAITS = {"ASD": BASE / "outputs/magma/grove.w35.genes.out",
          "ALLERGY": BASE / "outputs/magma/allergy.w35.genes.out"}
res = {}
for tname, path in TRAITS.items():
    Z, NS = load_z(path)
    use = [g for g in genes if g in Z and g in loeuf and cds.get(g, 0) > 0
           and allmean[gi[g]] > 0]
    idx = np.array([gi[g] for g in use])
    y = np.array([Z[g] for g in use])
    COV = [np.array([loeuf[g] for g in use]), np.log([cds[g] for g in use]),
           np.log(allmean[idx] + 1e-6), np.log([NS[g] + 1 for g in use])]
    res[tname] = dict(n_genes=len(use))
    for k in spec:
        b, se, t, p, n = ols(y, [spec[k][idx]] + COV)
        res[tname][k] = dict(beta=b[1], se=se[1], t=t[1], p=p[1], n=n)
    print(f"[{tname}] {len(use)} genes", file=sys.stderr)

json.dump(res, open(OUT / "allergy_symmetric.json", "w"), indent=1, default=float)

L = ["=" * 92, "SYMMETRIC TEST — cell-type localisation of ASD vs ALLERGY liability",
     "=" * 92,
     "identical model, gene-property conditioned (LOEUF, logCDS, logExpr, logNSNPs);",
     "atlas 1 (GSE294786 human infant PFC) cell-type specificity as predictor.",
     f"ASD n={res['ASD']['n_genes']} genes | ALLERGY n={res['ALLERGY']['n_genes']} genes", "",
     f"{'lineage':10s} {'ASD beta':>10s} {'ASD t':>8s} {'ASD p':>10s}   "
     f"{'ALL beta':>10s} {'ALL t':>8s} {'ALL p':>10s}"]
for k in ["neuron", "ol", "opc", "cop", "astro", "micro"]:
    if k not in res["ASD"]:
        continue
    a, b = res["ASD"][k], res["ALLERGY"][k]
    L.append(f"{k:10s} {a['beta']:>+10.4f} {a['t']:>+8.2f} {a['p']:>10.3g}   "
             f"{b['beta']:>+10.4f} {b['t']:>+8.2f} {b['p']:>10.3g}")
L += ["", "READING:",
      "  (a) allergy OL-positive & ASD OL-null -> immune->OL genetic bridge exists, ASD enters elsewhere",
      "  (b) both OL-null                      -> NO genetic form of the immune->OL bridge in either trait;",
      "                                           any immune->OL effect must be exposure-driven, not encoded",
      "  (c) both OL-positive                  -> shared route"]
txt = "\n".join(L)
(OUT / "ALLERGY_SYMMETRIC_REPORT.txt").write_text(txt + "\n")
print(txt)

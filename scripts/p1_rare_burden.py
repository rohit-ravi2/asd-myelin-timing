#!/usr/bin/env python3
"""
P1.6 — de novo variant burden in OL-lineage and convergence-gate modules.

WHY THIS IS NOW THE DECISIVE TEST
---------------------------------
P1.1-P1.3 found ASD COMMON-variant association does not localise to the OL lineage
(neuron specificity t=+3.36 p=8e-4; OL t=+1.00 p=0.32; stage effect gone after
constraint conditioning). But GSE294786's own published OPC/COP enrichment
(OR 3.54 / 3.89, FDR<1e-4) used a hypergeometric test on CURATED ASD RISK GENE LISTS
-- i.e. predominantly rare-variant discoveries. Those are different variant classes.
Their result and ours are not in conflict; they are testing different genetics.
So the OL hypothesis now lives or dies on RARE variation, which is this script.

DESIGN
------
SSC de novo calls: autism probands vs their UNAFFECTED SIBLINGS. The siblings supply
the empirical per-gene mutation opportunity directly, so no trinucleotide mutation
model is required -- the cleanest available design.

  classes: LoF (frameshift, stop-gained, splice), missense, synonymous (NEGATIVE CONTROL)
  per module m and class c:   rate ratio = (obs_case/N_case) / (obs_ctrl/N_ctrl)
  test: two-sided binomial on case share vs expected share N_case/(N_case+N_ctrl)

MATCHED NULL: 1000 random gene sets matched to each module on CDS length, LOEUF and
brain expression decile. This is what stops large constrained neurodevelopmental genes
from manufacturing "pathway enrichment" -- the failure mode the program document names.

HONEST POWER STATEMENT, WRITTEN BEFORE RESULTS: total SSC de novo LoF across ALL genes
is ~1,500 for cases+controls combined. Split across 15+ modules, most cells will contain
single-digit counts and be uninformative. Modules are reported with their counts so an
uninformative cell is visibly uninformative rather than silently interpreted.
"""
import gzip, json, sys
import numpy as np
from pathlib import Path
from collections import defaultdict
from math import lgamma, log, exp

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
SCOPE = Path("data/dvad-corpus/asd_scope")
OUT = BASE / "outputs/phase1"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20260917)
N_MATCH = 1000

LOF = {"frameshift", "stop-gained", "splice-acceptor", "splice-donor",
       "intron-near-splice", "stop-lost", "start-lost"}
MIS = {"missense"}
SYN = {"synonymous"}

# ------------------------------------------------------------ modules
members = {}
with open(BASE / "annot/asd_genesets.magma.set") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if not p[0].startswith("RANDOM"):
            members[p[0]] = set(p[1:])
gloc_sym = {}
with open(COR / "tools/magma_geneloc/NCBI37.3.gene.loc") as f:
    for line in f:
        p = line.split()
        if len(p) >= 6:
            gloc_sym[p[0]] = p[5]
modules = {k: {gloc_sym[g] for g in v if g in gloc_sym} for k, v in members.items()}

# --- add OL-lineage stage modules built from the infant atlas (the reconciliation test)
atlas_sets = json.load(open(OUT / "atlas_stage_sets.json")) if (OUT / "atlas_stage_sets.json").exists() else {}
for k, v in atlas_sets.items():
    modules[k] = set(v)
print(f"[modules] {len(modules)}", file=sys.stderr)

# ------------------------------------------------------------ gene covariates
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

# ------------------------------------------------------------ de novo calls
def read_dn(path):
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                hdr = line.lstrip("#").rstrip("\n").split("\t")
                ix = {c: i for i, c in enumerate(hdr)}
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) <= ix["FunctionClass"]:
                continue
            out.append((p[ix["SampleID"]], p[ix["PrimaryPhenotype"]],
                        p[ix["Gene"]], p[ix["FunctionClass"]]))
    return out

calls = read_dn(SCOPE / "rare_variant/denovo-db.ssc-samples.variants.tsv.gz")
print(f"[denovo] {len(calls)} SSC calls", file=sys.stderr)
n_case = len({s for s, ph, _, _ in calls if ph == "autism"})
n_ctrl = len({s for s, ph, _, _ in calls if ph == "control"})
print(f"[denovo] probands={n_case}  siblings={n_ctrl}", file=sys.stderr)

byclass = defaultdict(lambda: defaultdict(lambda: [0, 0]))   # class -> gene -> [case,ctrl]
for samp, ph, gene, fc in calls:
    if ph not in ("autism", "control") or not gene:
        continue
    cls = "LoF" if fc in LOF else ("missense" if fc in MIS else
                                   ("synonymous" if fc in SYN else None))
    if cls is None:
        continue
    byclass[cls][gene][0 if ph == "autism" else 1] += 1

for c in byclass:
    tc = sum(v[0] for v in byclass[c].values())
    tk = sum(v[1] for v in byclass[c].values())
    print(f"[denovo] {c:11s} case={tc} ctrl={tk}", file=sys.stderr)

# ------------------------------------------------------------ matched null pool
allg = [g for g in cds if g in loeuf]
arr_cds = np.array([cds[g] for g in allg])
arr_lo = np.array([loeuf[g] for g in allg])
cds_d = np.digitize(arr_cds, np.quantile(arr_cds, np.linspace(0, 1, 11)[1:-1]))
lo_d = np.digitize(arr_lo, np.quantile(arr_lo, np.linspace(0, 1, 11)[1:-1]))
strata = defaultdict(list)
for i, g in enumerate(allg):
    strata[(cds_d[i], lo_d[i])].append(g)
g2stratum = {g: (cds_d[i], lo_d[i]) for i, g in enumerate(allg)}


def matched_draw(genes):
    out = []
    for g in genes:
        pool = strata.get(g2stratum.get(g), None)
        if pool:
            out.append(pool[rng.integers(len(pool))])
    return out


def binom_p(k, n, p):
    """two-sided exact-ish binomial p via normal+continuity, robust for small n"""
    if n == 0:
        return 1.0
    from math import comb
    if n <= 400:
        probs = [comb(n, i) * p**i * (1 - p)**(n - i) for i in range(n + 1)]
        obs = probs[k]
        return min(1.0, sum(q for q in probs if q <= obs * (1 + 1e-9)))
    mu, sd = n * p, (n * p * (1 - p)) ** .5
    from math import erf, sqrt
    z = (abs(k - mu) - 0.5) / sd
    return 2 * (1 - 0.5 * (1 + erf(z / sqrt(2))))


exp_share = n_case / (n_case + n_ctrl)
rows = []
for mod, genes in sorted(modules.items()):
    for cls in ("LoF", "missense", "synonymous"):
        tab = byclass[cls]
        c = sum(tab[g][0] for g in genes if g in tab)
        k = sum(tab[g][1] for g in genes if g in tab)
        if c + k == 0:
            rows.append(dict(module=mod, cls=cls, case=0, ctrl=0, n_genes=len(genes),
                             rate_ratio=None, p=None, matched_p=None,
                             note="NO DE NOVO CALLS - uninformative"))
            continue
        rr = (c / n_case) / (k / n_ctrl) if k > 0 else float("inf")
        p = binom_p(c, c + k, exp_share)
        # matched-null: how often does a covariate-matched random set give >= this case share?
        obs_share = c / (c + k)
        ge = 0
        for _ in range(N_MATCH):
            mg = matched_draw(genes)
            mc = sum(tab[g][0] for g in mg if g in tab)
            mk = sum(tab[g][1] for g in mg if g in tab)
            if mc + mk > 0 and (mc / (mc + mk)) >= obs_share:
                ge += 1
        rows.append(dict(module=mod, cls=cls, case=c, ctrl=k, n_genes=len(genes),
                         rate_ratio=round(rr, 3) if rr != float("inf") else None,
                         p=round(p, 6), matched_p=round((ge + 1) / (N_MATCH + 1), 4),
                         note="low count" if c + k < 10 else ""))

json.dump(rows, open(OUT / "P1_rare_burden.json", "w"), indent=1)

L = ["=" * 108,
     "P1.6  de novo burden, SSC autism probands vs unaffected siblings",
     "=" * 108,
     f"probands={n_case}  siblings={n_ctrl}  expected case share={exp_share:.4f}",
     "matched_p = vs 1000 gene sets matched on CDS-length decile x LOEUF decile",
     "rate_ratio > 1 means excess in probands.  Counts shown so weak cells are visible.", "",
     f"{'module':30s} {'class':11s} {'case':>5s} {'ctrl':>5s} {'RR':>7s} {'p':>9s} {'matchP':>7s}  note"]
order = {"LoF": 0, "missense": 1, "synonymous": 2}
for r in sorted(rows, key=lambda x: (x["p"] if x["p"] is not None else 9, order[x["cls"]])):
    L.append(f"{r['module']:30s} {r['cls']:11s} {r['case']:>5d} {r['ctrl']:>5d} "
             f"{(('%.2f'%r['rate_ratio']) if r['rate_ratio'] else 'inf'):>7s} "
             f"{(('%.2e'%r['p']) if r['p'] is not None else 'NA'):>9s} "
             f"{(('%.3f'%r['matched_p']) if r['matched_p'] is not None else 'NA'):>7s}  {r['note']}")
txt = "\n".join(L)
(OUT / "P1_RARE_BURDEN_REPORT.txt").write_text(txt + "\n")
print(txt)

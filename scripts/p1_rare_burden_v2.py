#!/usr/bin/env python3
"""
P1.6 v2 — de novo burden, CARRIER-LEVEL.

WHY v2 EXISTS
-------------
v1 counted VARIANTS and tested them as independent binomial draws. They are not:
de novo calls cluster within samples, so one hypermutated proband can manufacture a
module-level "enrichment". v1's global negative control was clean (synonymous RR=0.99,
missense 1.03, LoF 1.40 = the known genome-wide ASD excess), but several MODULE-level
synonymous cells came out strongly case-enriched (OPC_VELM 58v9, GATE6 12v1), which is
biologically impossible for neutral variants and is the signature of exactly this problem.

v2 fixes it by counting CARRIERS: each sample contributes at most 1 per module per class.
That makes the unit of analysis the individual, which is the correct unit for burden.

Also adds:
  - BH FDR across all module x class tests
  - the synonymous row for every module, retained as a per-module negative control
  - carrier counts printed so weak cells stay visibly weak
"""
import gzip, json, sys
import numpy as np
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
SCOPE = Path("data/dvad-corpus/asd_scope")
OUT = BASE / "outputs/phase1"
rng = np.random.default_rng(20260917)
N_MATCH = 1000

LOF = {"frameshift", "stop-gained", "splice-acceptor", "splice-donor",
       "intron-near-splice", "stop-lost", "start-lost"}

# ---------------- modules
members = {}
with open(BASE / "annot/asd_genesets.magma.set") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if not p[0].startswith("RANDOM"):
            members[p[0]] = set(p[1:])
gsym = {}
with open(COR / "tools/magma_geneloc/NCBI37.3.gene.loc") as f:
    for line in f:
        p = line.split()
        if len(p) >= 6:
            gsym[p[0]] = p[5]
modules = {k: {gsym[g] for g in v if g in gsym} for k, v in members.items()}
modules.update({k: set(v) for k, v in json.load(open(OUT / "atlas_stage_sets.json")).items()})

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

# ---------------- carriers: sample -> class -> set(genes)
carr = defaultdict(lambda: defaultdict(set))
pheno = {}
with gzip.open(SCOPE / "rare_variant/denovo-db.ssc-samples.variants.tsv.gz", "rt") as f:
    for line in f:
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            h = line.lstrip("#").rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(h)}
            continue
        p = line.rstrip("\n").split("\t")
        ph = p[ix["PrimaryPhenotype"]]
        if ph not in ("autism", "control"):
            continue
        s, g, fc = p[ix["SampleID"]], p[ix["Gene"]], p[ix["FunctionClass"]]
        pheno[s] = ph
        cls = "LoF" if fc in LOF else ("missense" if fc == "missense" else
                                       ("synonymous" if fc == "synonymous" else None))
        if cls and g:
            carr[s][cls].add(g)

cases = [s for s, p in pheno.items() if p == "autism"]
ctrls = [s for s, p in pheno.items() if p == "control"]
NC, NK = len(cases), len(ctrls)
print(f"[carriers] probands={NC} siblings={NK}", file=sys.stderr)

# hypermutator diagnostic
for cls in ("LoF", "missense", "synonymous"):
    cc = sorted((len(carr[s][cls]) for s in cases), reverse=True)[:5]
    kk = sorted((len(carr[s][cls]) for s in ctrls), reverse=True)[:5]
    print(f"   {cls:11s} max genes/sample  case={cc}  ctrl={kk}", file=sys.stderr)

# ---------------- matched-null strata
allg = [g for g in cds if g in loeuf]
a_cds = np.array([cds[g] for g in allg]); a_lo = np.array([loeuf[g] for g in allg])
cd = np.digitize(a_cds, np.quantile(a_cds, np.linspace(0, 1, 11)[1:-1]))
ld = np.digitize(a_lo, np.quantile(a_lo, np.linspace(0, 1, 11)[1:-1]))
strata = defaultdict(list)
for i, g in enumerate(allg):
    strata[(cd[i], ld[i])].append(g)
g2s = {g: (cd[i], ld[i]) for i, g in enumerate(allg)}


def fisher(a, b, c, d):
    """two-sided Fisher exact"""
    from math import lgamma
    def lc(n, k): return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    tot = lc(a + b, a) + lc(c + d, c) - lc(n, a + c)
    obs = tot
    s = 0.0
    for i in range(0, min(a + b, a + c) + 1):
        j, k, l = a + b - i, a + c - i, d - (a - i)
        if j < 0 or k < 0 or l < 0:
            continue
        pr = lc(a + b, i) + lc(c + d, k) - lc(n, a + c)
        if pr <= obs + 1e-9:
            s += np.exp(pr)
    return min(1.0, s)


rows = []
for mod, genes in sorted(modules.items()):
    for cls in ("LoF", "missense", "synonymous"):
        a = sum(1 for s in cases if carr[s][cls] & genes)      # case carriers
        c = sum(1 for s in ctrls if carr[s][cls] & genes)      # ctrl carriers
        b, d = NC - a, NK - c
        if a + c == 0:
            rows.append(dict(module=mod, cls=cls, case_carr=0, ctrl_carr=0,
                             n_genes=len(genes), OR=None, p=None, matched_p=None,
                             note="NO CARRIERS - uninformative"))
            continue
        orr = ((a / NC) / (c / NK)) if c > 0 else float("inf")
        p = fisher(a, b, c, d)
        ge = 0
        for _ in range(N_MATCH):
            mg = {strata[g2s[g]][rng.integers(len(strata[g2s[g]]))]
                  for g in genes if g in g2s}
            ma = sum(1 for s in cases if carr[s][cls] & mg)
            mc = sum(1 for s in ctrls if carr[s][cls] & mg)
            if ma + mc > 0 and (ma / NC) / max(mc / NK, 1e-12) >= orr:
                ge += 1
        rows.append(dict(module=mod, cls=cls, case_carr=a, ctrl_carr=c,
                         n_genes=len(genes),
                         OR=round(orr, 3) if orr != float("inf") else None,
                         p=p, matched_p=(ge + 1) / (N_MATCH + 1),
                         note="low count" if a + c < 10 else ""))

# BH across all informative tests
inf = [r for r in rows if r["p"] is not None]
inf.sort(key=lambda x: x["p"])
m = len(inf)
for i, r in enumerate(inf, 1):
    r["q"] = min(1.0, r["p"] * m / i)
for i in range(m - 2, -1, -1):
    inf[i]["q"] = min(inf[i]["q"], inf[i + 1]["q"])

json.dump(rows, open(OUT / "P1_rare_burden_v2.json", "w"), indent=1, default=str)

L = ["=" * 116,
     "P1.6 v2  de novo CARRIER burden, SSC autism probands vs unaffected siblings",
     "=" * 116,
     f"probands={NC}  siblings={NK}   unit = individual (each sample counts once per module x class)",
     "OR>1 = excess carriers among probands.  matched_p vs 1000 CDS-length x LOEUF matched sets.",
     "synonymous row = per-module NEGATIVE CONTROL: it must be null or the module's result is suspect.",
     "",
     f"{'module':34s} {'class':11s} {'caseC':>6s} {'ctrlC':>6s} {'OR':>6s} {'p':>10s} {'q(BH)':>8s} {'matchP':>7s} note"]
for r in sorted(inf, key=lambda x: x["p"]):
    L.append(f"{r['module']:34s} {r['cls']:11s} {r['case_carr']:>6d} {r['ctrl_carr']:>6d} "
             f"{(('%.2f'%r['OR']) if r['OR'] else 'inf'):>6s} {r['p']:>10.2e} {r['q']:>8.3g} "
             f"{r['matched_p']:>7.3f} {r['note']}")
non = [r for r in rows if r["p"] is None]
if non:
    L.append("")
    L.append("uninformative (no carriers): " + ", ".join(f"{r['module']}/{r['cls']}" for r in non))
txt = "\n".join(L)
(OUT / "P1_RARE_BURDEN_V2_REPORT.txt").write_text(txt + "\n")
print(txt)

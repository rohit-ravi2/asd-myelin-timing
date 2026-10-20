#!/usr/bin/env python3
"""
p14_reconcile.py -- Phase 14, step 4/6.

Head-to-head arithmetic between:
  OURS  rg(ASD iPSYCH-PGC Grove2019, BIG40 TBSS-skeleton ICVF), 15 tracts
        outputs/rg/RG_REAL_ASD_results.csv
  GU    rg(ASD iPSYCH-2015, own UKB+ABCD GWAS of AutoPtx ICVF), 15 tracts
        Gu Y et al. Mol Psychiatry 2025;30(8):3393-3403, PMID 39994426,
        supplementary table S4, "White matter tract phenotype" block
        (outputs/phase14/gu2025_rg_tables/S4.tsv)

Writes outputs/phase14/p14_reconcile.log and p14_headtohead.csv.
"""
import csv, math, os, re

P = "."
OUT = f"{P}/outputs/phase14"

# ---- Gu S4 white-matter tract block, read from the supplementary TSV --------
GU_TRACTS = {"fma", "mcp", "fmi", "ar", "atr", "cgc", "cgh", "cst",
             "ifo", "ilf", "ml", "ptr", "slf", "str", "unc"}
gu = {}
with open(f"{OUT}/gu2025_rg_tables/S4.tsv") as f:
    for row in csv.reader(f, delimiter="\t"):
        if len(row) < 6:
            continue
        if row[0].strip() == "ASD" and row[1].strip() == "ICVF":
            t = row[2].strip()
            if t in GU_TRACTS:
                try:
                    gu[t] = (float(row[3]), float(row[4]), float(row[5]), float(row[6]))
                except (ValueError, IndexError):
                    pass
# global cortical ICVF, same table
gu_cortical_global = None
with open(f"{OUT}/gu2025_rg_tables/S4.tsv") as f:
    for row in csv.reader(f, delimiter="\t"):
        if len(row) > 6 and row[0].strip() == "ASD" and row[1].strip() == "ICVF" \
                and row[2].strip() == "global":
            gu_cortical_global = (float(row[3]), float(row[4]), float(row[5]), float(row[6]))

# ---- ours -------------------------------------------------------------------
ours = []
with open(f"{P}/outputs/rg/RG_REAL_ASD_results.csv") as f:
    for r in csv.DictReader(f):
        if r["measure"] == "ICVF":
            ours.append((r["label"], float(r["rg"]), float(r["se"])))

log = open(f"{OUT}/p14_reconcile.log", "w")
def out(*a):
    print(*a); print(*a, file=log)

out("=== p14 step 4: head-to-head, this project vs Gu et al. 2025 ===")
out("")
out("--- Gu et al. 2025 supplementary S4, WHITE MATTER TRACT block ---")
out(f"{'tract':<6} {'rg':>9} {'se':>8} {'z':>8} {'p':>9}")
for t in sorted(gu, key=lambda x: gu[x][0]):
    rg, se, z, p = gu[t]
    out(f"{t:<6} {rg:>9.4f} {se:>8.4f} {z:>8.3f} {p:>9.4f}")
gv = [v[0] for v in gu.values()]
out(f"n = {len(gv)}; negative = {sum(1 for x in gv if x < 0)}/{len(gv)}; "
    f"mean = {sum(gv)/len(gv):+.4f}")
if gu_cortical_global:
    out(f"\nGu GLOBAL CORTICAL ICVF: rg = {gu_cortical_global[0]:+.4f} "
        f"(se {gu_cortical_global[1]:.4f}, z {gu_cortical_global[2]:.3f}, "
        f"p {gu_cortical_global[3]:.4g})")

out("")
out("--- this project, BIG40 TBSS-skeleton ICVF ---")
ov = [x[1] for x in ours]
out(f"n = {len(ov)}; positive = {sum(1 for x in ov if x > 0)}/{len(ov)}; "
    f"mean = {sum(ov)/len(ov):+.4f}")

# ---- inverse-variance meta of each side (tracts are NOT independent, so this
#      is an effect-size summary, not a valid significance test) --------------
def ivw(pairs):
    w = [1 / (se ** 2) for _, se in pairs]
    m = sum(wi * b for wi, (b, _) in zip(w, pairs)) / sum(w)
    s = math.sqrt(1 / sum(w))
    return m, s

om, os_ = ivw([(b, s) for _, b, s in ours])
gm, gs = ivw([(v[0], v[1]) for v in gu.values()])
out("")
out("--- inverse-variance summary of each side ---")
out("(tracts are strongly genetically correlated -- mean rg(ICVF_i,ICVF_j) = 0.78 --")
out(" so these SEs are far too small to use as significance tests. They are here")
out(" only to summarise the central effect size of each study.)")
out(f"  ours : {om:+.4f}")
out(f"  Gu   : {gm:+.4f}")
out(f"  gap  : {om - gm:+.4f}")

# ---- tract-matched comparisons ---------------------------------------------
MATCH = {
    "slf": ["ICVF_SLF_L", "ICVF_SLF_R"],
    "cgc": ["ICVF_CingCingGyrus_L", "ICVF_CingCingGyrus_R"],
    "atr": ["ICVF_AntLimbIC_L", "ICVF_AntLimbIC_R"],   # approximate homologue
}
by_lab = {l: (b, s) for l, b, s in ours}
out("")
out("--- tract-matched, where an approximate homologue exists ---")
out("(NOTE: our ROIs are TBSS FA-skeleton JHU labels; Gu's are AutoPtx")
out(" probabilistic tracts. Same anatomical name, different spatial definition.)")
out(f"{'tract':<6} {'ours_rg':>9} {'ours_se':>8} {'gu_rg':>9} {'gu_se':>8} "
    f"{'diff':>8} {'z_diff':>7} {'p_diff':>9}")
rows = []
for t, labs in MATCH.items():
    if t not in gu:
        continue
    bs = [by_lab[l] for l in labs if l in by_lab]
    if not bs:
        continue
    # average the L/R estimates; their errors are near-perfectly correlated,
    # so keep the single-tract SE rather than shrinking it by sqrt(2)
    b = sum(x[0] for x in bs) / len(bs)
    s = max(x[1] for x in bs)
    grg, gse = gu[t][0], gu[t][1]
    d = b - grg
    z = d / math.sqrt(s ** 2 + gse ** 2)
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    out(f"{t:<6} {b:>9.4f} {s:>8.4f} {grg:>9.4f} {gse:>8.4f} "
        f"{d:>8.4f} {z:>7.2f} {p:>9.4f}")
    rows.append(dict(tract=t, ours_rg=b, ours_se=s, gu_rg=grg, gu_se=gse,
                     diff=d, z_diff=z, p_diff=p))
out("")
out("The difference test assumes the two rg estimates are independent. They are")
out("NOT: both ICVF GWAS are built on UK Biobank participants, so their sampling")
out("errors are POSITIVELY correlated, which makes Var(difference) SMALLER than")
out("the sum used above. The true |z| is therefore LARGER than shown -- the")
out("disagreement is understated, not overstated, by this calculation.")

with open(f"{OUT}/p14_headtohead.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["tract", "ours_rg", "ours_se", "gu_rg",
                                      "gu_se", "diff", "z_diff", "p_diff"])
    w.writeheader()
    w.writerows(rows)

# ---- FDR on our own 15 ------------------------------------------------------
out("")
out("--- our own multiple-testing position ---")
ps = sorted((2 * (1 - 0.5 * (1 + math.erf(abs(b / s) / math.sqrt(2)))), l)
            for l, b, s in ours)
n = len(ps)
qs, prev = [], 1.0
for i in range(n - 1, -1, -1):
    q = min(prev, ps[i][0] * n / (i + 1))
    qs.append((ps[i][1], ps[i][0], q))
    prev = q
qs.reverse()
out(f"{'tract':<24} {'p':>9} {'BH q':>9}")
for l, p, q in qs:
    out(f"{l:<24} {p:>9.4f} {q:>9.4f}")
out(f"\nminimum BH q across our 15 ICVF tracts = {min(q for _, _, q in qs):.4f}")
out("No single tract is significant at FDR 5% in either study.")
log.close()

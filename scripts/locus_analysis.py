#!/usr/bin/env python3
"""
Test C (APOE) + Test D (SREBF2) locus-level analysis in ASD GWAS.

Angles:
  C1. Direct lookup of the APOE-defining SNPs rs429358 (e4) / rs7412 (e2).
  C2. Regional association scan across the APOE locus (tight + broad windows).
  C3. Calibration: is the APOE region's min-p unusual vs MATCHED random regions
      drawn genome-wide with the same SNP count? (This is the honest test --
      "min-p = 0.002 somewhere in a 2 Mb window" is not evidence by itself.)
  D1. Same treatment for SREBF2, LRP1, MOG and every curated cholesterol gene.
  D2. Rank all cholesterol/myelin genes by regional min-p, report the full table.

Everything is reported: nulls, nominals, and misses alike.
"""
import gzip, json, random, sys, bisect
from pathlib import Path
from collections import defaultdict

random.seed(20260916)
OUT = Path("./outputs/locus")
OUT.mkdir(parents=True, exist_ok=True)
# dataset selectable: real ASD, or the mislabelled allergy file as a
# positive control that the pipeline detects genuine signal when present.
DATASETS = {
 "ASD_real":  dict(path="data/asd-corpus/genetics/ASD_Grove2019_REAL.gz",
                   cols=dict(CHR="CHR", BP="BP", P="P", RS="SNP", EFF="OR", SE="SE",
                             EA="A1", OA="A2"), sep=None),
 "ALLERGY_mislabelled": dict(path="data/dvad-corpus/asd_myelin/ASD_Grove2019.gz",
                   cols=dict(CHR="CHR", BP="BP", P="PVALUE", RS="RS_ID", EFF="BETA",
                             SE="SE", EA="EFFECT_ALLELE", OA="OTHER_ALLELE"), sep=None),
}
WHICH = sys.argv[1] if len(sys.argv) > 1 else "ASD_real"
DS = DATASETS[WHICH]
GROVE = DS["path"]
GLOC = "data/asd-corpus/tools/magma_geneloc/NCBI37.3.gene.loc"

# ------------------------------------------------------------------ gene coords
genes = {}
with open(GLOC) as f:
    for line in f:
        p = line.split()
        if len(p) >= 6:
            genes[p[5]] = (p[1], int(p[2]), int(p[3]))

TARGETS = {
    # name: (flank_kb, role)
    "APOE": (100, "CONTRAST-primary"),
    "TOMM40": (100, "CONTRAST"),
    "APOC1": (100, "CONTRAST"),
    "SREBF2": (100, "PRIMARY-finemap"),
    "SREBF1": (100, "PRIMARY"),
    "LRP1": (100, "PRIMARY-finemap"),
    "MOG": (100, "PRIMARY-finemap"),
    "HMGCR": (100, "PRIMARY"), "HMGCS1": (100, "PRIMARY"),
    "DHCR7": (100, "PRIMARY-syndromic"), "DHCR24": (100, "PRIMARY"),
    "NPC1": (100, "PRIMARY-syndromic"), "NPC2": (100, "PRIMARY"),
    "ABCA1": (100, "PRIMARY"), "ABCG1": (100, "PRIMARY"), "ABCA2": (100, "PRIMARY"),
    "LDLR": (100, "PRIMARY"), "CYP46A1": (100, "PRIMARY"),
    "SQLE": (100, "PRIMARY"), "FDFT1": (100, "PRIMARY"), "LSS": (100, "PRIMARY"),
    "MVK": (100, "PRIMARY"), "MVD": (100, "PRIMARY"), "IDI1": (100, "PRIMARY"),
    "EBP": (100, "PRIMARY"), "SC5D": (100, "PRIMARY"), "MSMO1": (100, "PRIMARY"),
    "CYP51A1": (100, "PRIMARY"), "NSDHL": (100, "PRIMARY"),
    "SCAP": (100, "PRIMARY"), "INSIG1": (100, "PRIMARY"), "INSIG2": (100, "PRIMARY"),
    "SOAT1": (100, "PRIMARY"),
    # myelin structural
    "MBP": (100, "MYELIN"), "PLP1": (100, "MYELIN"), "MAG": (100, "MYELIN"),
    "CNP": (100, "MYELIN"), "MOBP": (100, "MYELIN"), "MAL": (100, "MYELIN"),
    "SOX10": (100, "MYELIN"), "OLIG1": (100, "MYELIN"), "OLIG2": (100, "MYELIN"),
    "MYRF": (100, "MYELIN"), "QKI": (100, "MYELIN"), "UGT8": (100, "MYELIN"),
    "PMP22": (100, "MYELIN"), "CLDN11": (100, "MYELIN"), "GJC2": (100, "MYELIN"),
    # positive controls: known ASD GWAS loci from Grove 2019
    "MACROD2": (100, "POSCTRL-asd"), "XRN2": (100, "POSCTRL-asd"),
    "KIZ": (100, "POSCTRL-asd"), "KMT2E": (100, "POSCTRL-asd"),
    "SORCS3": (100, "POSCTRL-asd"), "NEGR1": (100, "POSCTRL-asd"),
    "PTBP2": (100, "POSCTRL-asd"), "CADPS": (100, "POSCTRL-asd"),
}
DEFINING = {"rs429358": "APOE e4-defining", "rs7412": "APOE e2-defining"}

# ------------------------------------------------------------------ read GWAS
print(f"[read] {WHICH}: {GROVE}", file=sys.stderr)
bychr = defaultdict(list)      # chr -> list of (bp, p, rsid, beta, se)
defining_hits = {}
n = 0
with gzip.open(GROVE, "rt") as f:
    hdr = f.readline().split()
    ix = {c: i for i, c in enumerate(hdr)}
    C = DS["cols"]
    for line in f:
        p = line.split()
        try:
            pv = float(p[ix[C["P"]]])
        except (ValueError, IndexError):
            continue
        c, bp = p[ix[C["CHR"]]], int(p[ix[C["BP"]]])
        rs = p[ix[C["RS"]]]
        bychr[c].append((bp, pv, rs, p[ix[C["EFF"]]], p[ix[C["SE"]]]))
        if rs in DEFINING:
            defining_hits[rs] = dict(chr=c, bp=bp, p=pv, beta=p[ix[C["EFF"]]],
                                     se=p[ix[C["SE"]]], ea=p[ix[C["EA"]]],
                                     oa=p[ix[C["OA"]]])
        n += 1
print(f"[read] {n} SNPs", file=sys.stderr)
for c in bychr:
    bychr[c].sort()
pos = {c: [x[0] for x in v] for c, v in bychr.items()}


def region(c, lo, hi):
    if c not in bychr:
        return []
    a = bisect.bisect_left(pos[c], lo)
    b = bisect.bisect_right(pos[c], hi)
    return bychr[c][a:b]


# ------------------------------------------------------------------ C1 defining SNPs
res = {"C1_apoe_defining_snps": {}}
for rs, lbl in DEFINING.items():
    if rs in defining_hits:
        h = defining_hits[rs]
        res["C1_apoe_defining_snps"][rs] = dict(label=lbl, **h)
    else:
        res["C1_apoe_defining_snps"][rs] = dict(label=lbl, status="NOT_IN_GWAS")

# ------------------------------------------------------------------ C2/D1 regions
def scan(name, flank_kb):
    if name not in genes:
        return None
    c, s, e = genes[name]
    lo, hi = s - flank_kb * 1000, e + flank_kb * 1000
    snps = region(c, lo, hi)
    if not snps:
        return dict(gene=name, chr=c, start=s, end=e, n_snps=0, status="NO_SNPS")
    best = min(snps, key=lambda x: x[1])
    return dict(gene=name, chr=c, start=s, end=e, window_kb=flank_kb,
                n_snps=len(snps), min_p=best[1], min_p_rsid=best[2],
                min_p_bp=best[0], beta=best[3], se=best[4],
                n_p_lt_1e3=sum(1 for x in snps if x[1] < 1e-3),
                n_p_lt_1e4=sum(1 for x in snps if x[1] < 1e-4),
                n_p_lt_5e8=sum(1 for x in snps if x[1] < 5e-8))


regions = {}
for g, (fk, role) in TARGETS.items():
    r = scan(g, fk)
    if r:
        r["role"] = role
        regions[g] = r
    else:
        regions[g] = dict(gene=g, status="NOT_IN_GENE_LOC", role=role)

# ------------------------------------------------------------------ C3 calibration
# For each observed region, draw N_MATCH random genome regions with a similar SNP
# count and ask how often their min-p beats the observed one. This converts
# "min-p = X" into an empirical regional p-value.
N_MATCH = 2000
allchr = [c for c in bychr if c.isdigit()]


def empirical_region_p(n_snps, obs_min_p, n_draw=N_MATCH):
    if n_snps < 5:
        return None
    better = 0
    for _ in range(n_draw):
        c = random.choice(allchr)
        v = bychr[c]
        if len(v) <= n_snps:
            continue
        i = random.randrange(0, len(v) - n_snps)
        m = min(x[1] for x in v[i:i + n_snps])
        if m <= obs_min_p:
            better += 1
    return (better + 1) / (n_draw + 1)


print("[calib] empirical regional p-values ...", file=sys.stderr)
for g, r in regions.items():
    if r.get("n_snps", 0) >= 5:
        r["empirical_regional_p"] = empirical_region_p(r["n_snps"], r["min_p"])
        r["n_matched_draws"] = N_MATCH

res["C2_D1_regions"] = regions

# ------------------------------------------------------------------ write
with open(OUT / f"locus_results_{WHICH}.json", "w") as f:
    json.dump(res, f, indent=1)

# ------------------------------------------------------------------ report
lines = []
lines.append("=" * 100)
lines.append(f"TEST C + D: APOE / cholesterol / myelin locus analysis  [dataset: {WHICH}]")
lines.append("=" * 100)
lines.append("")
lines.append("--- C1: APOE-defining SNPs -------------------------------------------")
for rs, d in res["C1_apoe_defining_snps"].items():
    if d.get("status") == "NOT_IN_GWAS":
        lines.append(f"  {rs:12s} {d['label']:22s} NOT PRESENT in Grove2019 sumstats")
    else:
        lines.append(f"  {rs:12s} {d['label']:22s} chr{d['chr']}:{d['bp']}  "
                     f"P={d['p']:.4g}  beta={d['beta']}  se={d['se']}  "
                     f"({d['ea']}/{d['oa']})")
lines.append("")
lines.append("--- C2/D1: regional scan (gene +/-100kb), sorted by min-p ------------")
lines.append(f"{'gene':12s} {'role':20s} {'chr':>3s} {'nSNP':>6s} {'min_p':>10s} "
             f"{'emp_regP':>9s} {'p<1e-3':>7s} {'p<1e-4':>7s} {'GWS':>4s}  lead")
ok = [r for r in regions.values() if r.get("n_snps", 0) > 0]
for r in sorted(ok, key=lambda x: x["min_p"]):
    ep = r.get("empirical_regional_p")
    lines.append(f"{r['gene']:12s} {r['role']:20s} {r['chr']:>3s} {r['n_snps']:>6d} "
                 f"{r['min_p']:>10.3g} {('%.3f'%ep) if ep else 'NA':>9s} "
                 f"{r['n_p_lt_1e3']:>7d} {r['n_p_lt_1e4']:>7d} {r['n_p_lt_5e8']:>4d}  "
                 f"{r['min_p_rsid']}")
missing = [r["gene"] for r in regions.values() if r.get("n_snps", 0) == 0]
if missing:
    lines.append("")
    lines.append(f"  no SNPs / not in gene.loc: {', '.join(missing)}")

txt = "\n".join(lines)
(OUT / f"TESTCD_locus_report_{WHICH}.txt").write_text(txt + "\n")
print(txt)

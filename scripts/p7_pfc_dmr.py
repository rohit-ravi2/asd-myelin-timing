#!/usr/bin/env python3
"""
§17-§22 — ASD PFC single-nucleus methylation (GSE298243), 2.5M tested 500bp bins.

WHY THIS DESIGN
---------------
§18 forbids "map DMR to nearest gene and call it enrichment". The table already contains
every TESTED bin, so the correct background is the tested set itself, stratified by the
things that drive DMR detection: CELL TYPE, CYTOSINE CONTEXT (CpG vs CpA -- non-CpG
methylation is a major neuronal/glial feature), and CHROMOSOME.

For each (cell type x module) we ask: among bins annotated to module genes, is the rate of
ASD-diagnosis-significant bins higher than among matched non-module bins from the SAME cell
type and cytosine context? Fisher exact + BH.

§21 asks where the burden sits across cell types -- that is the primary output, because if
immune-linked epigenomic change sits in microglia rather than OL, the mechanism differs.

NOTE: this dataset is postmortem and its ASD effects are small relative to age effects.
Absence of enrichment here does NOT refute a prenatal regulatory mechanism (§stated in
ledger). It constrains the postmortem PFC endpoint only.
"""
import gzip, json, sys
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
SCOPE = Path("data/dvad-corpus/asd_scope")
DMR = SCOPE / "celltype/GSE298243_Age_DMRs_PrespecifiedRegions.tsv.gz"
OUT = BASE / "outputs/phase7"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------- modules (reuse the frozen panel + immune submodules)
import re
gmt = {}
for fn in ["c2.cp.v2024.1.Hs.symbols.gmt", "c5.go.bp.v2024.1.Hs.symbols.gmt",
           "c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    p = COR / "pathways" / fn
    if p.exists():
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if len(f) > 2:
                gmt[f[0]] = set(f[2:])


def U(*pats):
    g = set()
    for n, m in gmt.items():
        if any(re.search(p, n, re.I) for p in pats):
            g |= m
    return g


MOD = {
    "OL_MATURATION": U(r"MYELIN", r"ENSHEATHMENT_OF_NEURONS", r"OLIGODENDROCYTE"),
    "CHOL_SYNTHESIS": U(r"CHOLESTEROL_BIOSYNTH", r"STEROL_BIOSYNTH", r"MEVALONATE"),
    "LIPID_TRANSPORT": U(r"LIPID_TRANSPORT", r"CHOLESTEROL_TRANSPORT", r"LIPOPROTEIN_PARTICLE"),
    "SPHINGO_VLCFA": U(r"SPHINGOLIPID", r"CERAMIDE", r"VERY_LONG_CHAIN"),
    "MITOCHONDRIAL": U(r"MITOCHONDRI"),
    "CHROMATIN_STATE": U(r"CHROMATIN_REMODEL", r"CHROMATIN_ORGANIZATION", r"HISTONE_MODIFICATION"),
    "SYNAPTIC": U(r"GOBP_SYNAP", r"GOCC_POSTSYNAPTIC_DENSITY"),
    "INTERFERON": U(r"INTERFERON"),
    "NFKB_TNF": U(r"NF_KAPPAB", r"TUMOR_NECROSIS_FACTOR"),
    "IL6_STAT3": U(r"INTERLEUKIN_6", r"STAT3"),
    "IL17_TH17": U(r"INTERLEUKIN_17", r"TH17"),
    "COMPLEMENT": U(r"COMPLEMENT"),
    "MICROGLIA_FN": U(r"MICROGLI", r"MACROPHAGE_ACTIVATION"),
    "NEG_OLFACTORY": U(r"OLFACTORY"),
    "NEG_KERATIN": U(r"KERATIN", r"CORNIFIED"),
}
import csv as _csv
sf = set()
with open(SCOPE / "sfari/SFARI-Gene_genes_2026-09-10.csv") as f:
    for r in _csv.DictReader(f):
        g = r.get("gene-symbol") or r.get("gene_symbol")
        if g:
            sf.add(g)
MOD["REF_SFARI"] = sf
print(f"[modules] {len(MOD)}", file=sys.stderr)

# ---------------- single streaming pass
# counts[(celltype, cyto, module)] = [sig_in, tot_in, sig_out, tot_out]
counts = defaultdict(lambda: [0, 0, 0, 0])
ct_tot = defaultdict(lambda: [0, 0])          # celltype -> [sig, tot]
ctx_tot = defaultdict(lambda: [0, 0])
n = 0
with gzip.open(DMR, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    iCT, iCX, iSym = ix["CellType"], ix["CytoContext"], ix["Symbol"]
    iDx = ix["ReportSignif_Diagnosis"]
    iInt = ix["ReportSignif_Interaction"]
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= iInt:
            continue
        n += 1
        ct, cx = p[iCT], p[iCX]
        sym = p[iSym].strip()
        sig = 1 if p[iDx].strip() == "*" else 0
        ct_tot[ct][0] += sig; ct_tot[ct][1] += 1
        ctx_tot[cx][0] += sig; ctx_tot[cx][1] += 1
        if not sym:
            continue
        for m, genes in MOD.items():
            k = (ct, cx, m)
            if sym in genes:
                counts[k][0] += sig; counts[k][1] += 1
            else:
                counts[k][2] += sig; counts[k][3] += 1
        if n % 500000 == 0:
            print(f"   [{n}] rows", file=sys.stderr)
print(f"[scan] {n} rows", file=sys.stderr)


def fisher(a, b, c, d):
    from math import lgamma, exp
    def lc(nn, kk): return lgamma(nn + 1) - lgamma(kk + 1) - lgamma(nn - kk + 1)
    N = a + b + c + d
    obs = lc(a + b, a) + lc(c + d, c) - lc(N, a + c)
    s = 0.0
    for i in range(0, min(a + b, a + c) + 1):
        k = a + c - i
        if k < 0 or k > c + d:
            continue
        pr = lc(a + b, i) + lc(c + d, k) - lc(N, a + c)
        if pr <= obs + 1e-9:
            s += exp(pr)
    return min(1.0, s)


rows = []
for (ct, cx, m), (si, ti, so, to) in counts.items():
    if ti < 50 or si + so == 0:
        continue
    ri, ro = si / ti, so / to
    orr = (ri / ro) if ro > 0 else float("inf")
    p = fisher(si, ti - si, so, to - so)
    rows.append(dict(cell_type=ct, cyto=cx, module=m, sig_in=si, tot_in=ti,
                     rate_in=ri, rate_out=ro, OR=orr, p=p))
rows.sort(key=lambda x: x["p"])
mm = len(rows)
for i, r in enumerate(rows, 1):
    r["q"] = min(1.0, r["p"] * mm / i)
for i in range(mm - 2, -1, -1):
    rows[i]["q"] = min(rows[i]["q"], rows[i + 1]["q"])

json.dump(dict(rows=rows, celltype_totals={k: v for k, v in ct_tot.items()},
               ctx_totals={k: v for k, v in ctx_tot.items()}, n_rows=n),
          open(OUT / "pfc_dmr_enrichment.json", "w"), indent=1, default=float)

L = ["=" * 108, "§17-22 ASD PFC single-nucleus methylation — diagnosis-DMR enrichment",
     "=" * 108,
     f"{n} tested 500bp bins. Background = tested bins of the SAME cell type and cytosine",
     "context. This is not nearest-gene mapping.", "",
     "§21 PRIMARY — ASD-diagnosis DMR rate by cell type:",
     f"{'cell type':16s} {'sig':>8s} {'tested':>10s} {'rate':>8s}"]
for ct, (s, t) in sorted(ct_tot.items(), key=lambda x: -(x[1][0] / max(x[1][1], 1))):
    L.append(f"{ct:16s} {s:>8d} {t:>10d} {s/max(t,1):>8.5f}")
L += ["", "by cytosine context:"]
for cx, (s, t) in sorted(ctx_tot.items()):
    L.append(f"   {cx:6s} sig={s} tested={t} rate={s/max(t,1):.5f}")
L += ["", "§18-20 module enrichment (BH q<0.1 shown; full table in JSON):",
      f"{'cell':12s} {'cyto':5s} {'module':18s} {'sigIn':>6s} {'totIn':>7s} "
      f"{'rateIn':>8s} {'rateOut':>8s} {'OR':>6s} {'p':>9s} {'q':>9s}"]
shown = 0
for r in rows:
    if r["q"] < 0.1:
        L.append(f"{r['cell_type']:12s} {r['cyto']:5s} {r['module']:18s} {r['sig_in']:>6d} "
                 f"{r['tot_in']:>7d} {r['rate_in']:>8.5f} {r['rate_out']:>8.5f} "
                 f"{r['OR']:>6.2f} {r['p']:>9.2e} {r['q']:>9.2e}")
        shown += 1
if not shown:
    L.append("   (none at q<0.1)")
L += ["", "CAVEAT: postmortem PFC; ASD effects small relative to age effects. A null here",
      "does NOT refute a prenatal regulatory mechanism -- it constrains the adult endpoint."]
txt = "\n".join(L)
(OUT / "PFC_DMR_REPORT.txt").write_text(txt + "\n")
print(txt)

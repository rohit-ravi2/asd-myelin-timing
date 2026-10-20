#!/usr/bin/env python3
"""
Does the strongest version survive? — characterising the OPC epigenetic signal.

The reformulated model rests on OPC being a convergence node perturbed from upstream, with
an immune/metabolic insult at the OPC->COP window as the candidate cause. That model makes
CHECKABLE predictions about the CONTENT of the OPC DMR signal:

 P1  It should be OPC-SPECIFIC, not a generic OL-lineage phenomenon
     (already suggestive: OPC 0.0217 vs Oligo 0.0130 -- tested formally here).
 P2  If an immune insult drives it, OPC DMR genes should be enriched for immune-response
     programmes RELATIVE TO THE OPC-TESTED BACKGROUND (the previous run tested modules
     pooled across cell types and found nothing at q<0.1; this tests within OPC).
 P3  If instead it is developmental/state-control, chromatin and neurodevelopmental
     programmes should dominate -- which would favour Pipeline B over the immune model.
 P4  Direction matters. A coherent insult should give a directional bias
     (hyper- or hypo-methylation), not a symmetric scatter.
 P5  SFARI overlap. The table carries its own SFARI flag; if OPC DMRs are enriched for
     SFARI genes that argues the signal is ASD-genetic-adjacent rather than exposure-driven.

Background is always the OPC-tested CpG bins themselves, so gene-set results are not
confounded by which bins were testable in OPC.
"""
import gzip, json, re, sys
from pathlib import Path
from collections import defaultdict, Counter
from math import lgamma, exp

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
SCOPE = Path("data/dvad-corpus/asd_scope")
DMR = SCOPE / "celltype/GSE298243_Age_DMRs_PrespecifiedRegions.tsv.gz"
OUT = BASE / "outputs/phase13"
OUT.mkdir(parents=True, exist_ok=True)

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
    "IMMUNE_cytokine": U(r"CYTOKINE_MEDIATED", r"CYTOKINE_PRODUCTION", r"INFLAMMATORY_RESPONSE"),
    "IMMUNE_interferon": U(r"INTERFERON"),
    "IMMUNE_NFKB_TNF": U(r"NF_KAPPAB", r"TUMOR_NECROSIS_FACTOR"),
    "IMMUNE_complement": U(r"COMPLEMENT"),
    "IMMUNE_innate_TLR": U(r"TOLL_LIKE_RECEPTOR", r"INNATE_IMMUNE"),
    "PURINE_adenosine": U(r"PURINE_", r"ADENOSINE", r"NUCLEOTIDE_METABOL"),
    "OL_maturation": U(r"MYELIN", r"OLIGODENDROCYTE", r"ENSHEATHMENT"),
    "LIPID_synthesis": U(r"CHOLESTEROL_BIOSYNTH", r"STEROL_BIOSYNTH", r"SPHINGOLIPID"),
    "MITOCHONDRIAL": U(r"MITOCHONDRI"),
    "CHROMATIN_state": U(r"CHROMATIN_REMODEL", r"CHROMATIN_ORGANIZATION", r"HISTONE_MODIFICATION"),
    "NEURODEV": U(r"NEURON_DIFFERENTIATION", r"NERVOUS_SYSTEM_DEVELOPMENT", r"AXON_GUIDANCE"),
    "SYNAPTIC": U(r"GOBP_SYNAP", r"GOCC_POSTSYNAPTIC_DENSITY"),
    "CELL_CYCLE": U(r"CELL_CYCLE", r"MITOTIC"),
    "NEG_OLFACTORY": U(r"OLFACTORY"),
    "NEG_KERATIN": U(r"KERATIN", r"CORNIFIED"),
}
print(f"[modules] {len(MOD)}", file=sys.stderr)

# ------------------------------------------------------------------ scan
cell_lineage = {"OPC", "Oligo"}
recs = defaultdict(list)          # celltype -> list of (sym, sig, coefEst, sfari)
gene_sig = defaultdict(lambda: [0, 0])   # OPC only: gene -> [sig, tested]
dir_counts = defaultdict(lambda: [0, 0])  # celltype -> [hyper, hypo] among significant
n = 0
with gzip.open(DMR, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    iCT, iCX, iSym, iDx = ix["CellType"], ix["CytoContext"], ix["Symbol"], ix["ReportSignif_Diagnosis"]
    iCoef = ix["coefEst_DiagnosisASD"]
    iSF = ix["SFARI"]
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(iCoef, iSF):
            continue
        if p[iCX] != "CpG":
            continue
        n += 1
        ct = p[iCT]; sym = p[iSym].strip(); sig = 1 if p[iDx].strip() == "*" else 0
        sf = p[iSF].strip()
        try:
            coef = float(p[iCoef])
        except ValueError:
            coef = None
        if sig and coef is not None:
            dir_counts[ct][0 if coef > 0 else 1] += 1
        if ct == "OPC" and sym:
            gene_sig[sym][0] += sig; gene_sig[sym][1] += 1
        if sig and sym:
            recs[ct].append((sym, coef, sf))
print(f"[scan] {n} CpG bins", file=sys.stderr)


def fisher(a, b, c, d):
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


# ---------------- P2/P3: module enrichment WITHIN OPC, background = OPC-tested bins
rows = []
for m, genes in MOD.items():
    si = sum(v[0] for g, v in gene_sig.items() if g in genes)
    ti = sum(v[1] for g, v in gene_sig.items() if g in genes)
    so = sum(v[0] for g, v in gene_sig.items() if g not in genes)
    to = sum(v[1] for g, v in gene_sig.items() if g not in genes)
    if ti < 30 or si + so == 0:
        continue
    rr = (si / ti) / (so / to) if so > 0 else float("inf")
    rows.append(dict(module=m, sig=si, tested=ti, rate=si / ti, rate_bg=so / to,
                     OR=rr, p=fisher(si, ti - si, so, to - so)))
rows.sort(key=lambda x: x["p"])
mm = len(rows)
for i, r in enumerate(rows, 1):
    r["q"] = min(1.0, r["p"] * mm / i)
for i in range(mm - 2, -1, -1):
    rows[i]["q"] = min(rows[i]["q"], rows[i + 1]["q"])

# ---------------- P5 SFARI flag enrichment within OPC
sf_sig = sum(1 for s, c, f in recs["OPC"] if f)
sf_tot = sum(1 for s, c, f in recs["OPC"])

# ---------------- top OPC genes by significant-bin count
top = sorted(((g, v[0], v[1]) for g, v in gene_sig.items() if v[0] > 0),
             key=lambda x: -x[1])[:30]

json.dump(dict(modules=rows, direction={k: v for k, v in dir_counts.items()},
               top_genes=top, sfari_sig=sf_sig, sfari_tot=sf_tot),
          open(OUT / "opc_dmr_character.json", "w"), indent=1, default=float)

L = ["=" * 100, "DOES THE STRONGEST VERSION HOLD? — content of the OPC epigenetic signal",
     "=" * 100,
     "background for every test = the OPC-tested CpG bins themselves", "",
     "P4 — DIRECTION of significant diagnosis effects (hyper : hypo):"]
for ct, (hy, ho) in sorted(dir_counts.items(), key=lambda x: -(x[1][0] + x[1][1])):
    tot = hy + ho
    L.append(f"   {ct:10s} hyper={hy:>5d} hypo={ho:>5d}  ({hy/max(tot,1)*100:.1f}% hyper)")
L += ["", "P2/P3 — module enrichment WITHIN OPC (immune vs developmental):",
      f"{'module':22s} {'sig':>5s} {'tested':>7s} {'rate':>8s} {'bg rate':>8s} {'OR':>6s} {'p':>9s} {'q':>8s}"]
for r in rows:
    L.append(f"{r['module']:22s} {r['sig']:>5d} {r['tested']:>7d} {r['rate']:>8.4f} "
             f"{r['rate_bg']:>8.4f} {r['OR']:>6.2f} {r['p']:>9.3g} {r['q']:>8.3g}")
L += ["", f"P5 — SFARI flag among OPC significant bins: {sf_sig}/{sf_tot} "
      f"({sf_sig/max(sf_tot,1)*100:.1f}%)", "",
      "top OPC genes by number of significant CpG bins:"]
L.append("   " + ", ".join(f"{g}({s}/{t})" for g, s, t in top[:20]))
txt = "\n".join(L)
(OUT / "OPC_DMR_CHARACTER_REPORT.txt").write_text(txt + "\n")
print(txt)

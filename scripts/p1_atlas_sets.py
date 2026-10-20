#!/usr/bin/env python3
"""
Derive OL-lineage stage-specific gene sets from the infant PFC atlas (GSE294786),
plus the neurogenic->gliogenic fate-switch set for gate 6.

These become modules in the rare-variant burden test (P1.6), which is the direct
reconciliation with GSE294786's own published OPC/COP ASD-risk-gene enrichment
(hypergeometric on curated risk-gene lists; OPC OR 3.54, COP OR 3.89).
If that enrichment is real and rare-variant driven, it should reappear here as a
de novo LoF excess in OPC/COP-specific genes. If it does not, their result is
carried by gene-list curation rather than by de novo burden.
"""
import gzip, json, sys
import numpy as np
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parents[1]
GSE = Path("data/asd-corpus/celltype/GSE294786")
OUT = BASE / "outputs/phase1"
OUT.mkdir(parents=True, exist_ok=True)
TOPN = 300

meta = {}
with gzip.open(GSE / "GSE294786_all_meta_progenitor.tsv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    if hdr[0] == "":
        hdr[0] = "CELLID"
    ix = {c: i for i, c in enumerate(hdr)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= len(hdr):
            meta[p[ix["CELLID"]]] = (p[ix["Species"]], p[ix["Age"]], p[ix["cell_type_final"]])
HI = {k: v[2] for k, v in meta.items() if v[0] == "Hg" and v[1] == "Infant"}

with gzip.open(GSE / "GSE294786_all_counts_transposed.tsv.gz", "rt") as f:
    genes = f.readline().rstrip("\n").split("\t")[1:]
    G = len(genes)
    sums = defaultdict(lambda: np.zeros(G)); n = defaultdict(int)
    for line in f:
        tab = line.find("\t")
        ct = HI.get(line[:tab])
        if ct is None:
            continue
        v = np.fromstring(line[tab + 1:], sep="\t")
        if v.size != G:
            continue
        t = v.sum()
        if t <= 0:
            continue
        sums[ct] += np.log1p(v / t * 1e4); n[ct] += 1
means = {k: sums[k] / n[k] for k in n}
print(f"[atlas] {len(means)} cell types", file=sys.stderr)
np.save(OUT / "atlas_means.npy", np.array([means[k] for k in sorted(means)]))
json.dump({"types": sorted(means), "genes": genes}, open(OUT / "atlas_index.json", "w"))

allm = np.mean([means[k] for k in means], axis=0)
sets = {}


def top_specific(types, name, exclude=()):
    a = np.mean([means[t] for t in types if t in means], axis=0)
    other = [means[t] for t in means if t not in types and t not in exclude]
    b = np.mean(other, axis=0)
    score = np.log2((a + 1e-3) / (b + 1e-3))
    score[a < 0.05] = -np.inf          # require real expression in the focal type
    idx = np.argsort(score)[::-1][:TOPN]
    sets[name] = [genes[i] for i in idx]
    print(f"   {name:26s} n={len(sets[name])}  top: {', '.join(sets[name][:8])}", file=sys.stderr)


top_specific({"OPC"}, "ATLAS_OPC_specific")
top_specific({"COP"}, "ATLAS_COP_specific")
top_specific({"OPC", "COP"}, "ATLAS_OLlineage_specific")
top_specific({"Astrocyte"}, "ATLAS_Astrocyte_specific")
top_specific({"Mic-P2RY12"}, "ATLAS_Microglia_specific")
NEUR = {"L2/L3", "L2/L3-BS", "L4", "L5/L6-CC", "L5/L6", "IN-PV", "IN-SST", "IN-VIP", "IN-SV2C"}
top_specific(NEUR, "ATLAS_Neuron_specific")
if "Progenitor" in means:
    top_specific({"Progenitor"}, "ATLAS_Progenitor_specific")

# ---- gate 6: neurogenic -> gliogenic fate switch (curated, hg-symbol)
sets["GATE6_NEUROGENIC_GLIOGENIC_SWITCH"] = """NFIA NFIB NFIX SOX9 SOX8 SOX10 HES1 HES5
HEY1 HEY2 NOTCH1 NOTCH2 RBPJ DLL1 JAG1 STAT3 JAK1 JAK2 EGFR FGFR1 FGFR2 FGFR3 OLIG1 OLIG2
ASCL1 NEUROG1 NEUROG2 NEUROD1 NEUROD2 EOMES PAX6 GLI3 ZBTB20 TCF7L2 RFX4 LHX2 POU3F2 POU3F3
SALL3 ID2 ID4 CHD7 CHD8 EZH2 SUZ12 EED KDM6B SETD2 KMT2A DNMT1 DNMT3A TET1 TET2 TET3
CNTFR LIF IL6ST BMP2 BMP4 BMPR1A SMAD1 SMAD4 CTNNB1 GSK3B""".split()
print(f"   GATE6 switch set n={len(sets['GATE6_NEUROGENIC_GLIOGENIC_SWITCH'])}", file=sys.stderr)

json.dump(sets, open(OUT / "atlas_stage_sets.json", "w"), indent=1)
print(f"[done] {len(sets)} sets -> atlas_stage_sets.json", file=sys.stderr)

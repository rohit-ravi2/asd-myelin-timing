#!/usr/bin/env python3
"""
Does the cytokine->oligodendrocyte arrow have a RECEPTOR substrate in human infant PFC?

The claim under test: "substantial immune responses release inflammatory cytokines that are
known to affect oligodendrocytes." Links 1-2 are solid. This script tests a necessary
condition for the DIRECT form of link 3, in the correct species / region / developmental
window: are the receptors for those cytokines actually expressed on human infant OPC and COP?

If a receptor is absent from OPC/COP at this stage, a direct cytokine->OL-progenitor effect
requires either (a) a different developmental window, (b) an indirect route via microglia /
astrocytes / endothelium, or (c) a non-canonical receptor. Any of those is fine -- but it
changes which experiment is decisive, and it is not the same claim.

Comparison is WITHIN-CELL-TYPE against microglia, which must express these receptors highly.
Microglia therefore act as the positive control for the measurement itself.
"""
import json, sys
import numpy as np
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
P1 = BASE / "outputs/phase1"
OUT = BASE / "outputs/phase4"
OUT.mkdir(parents=True, exist_ok=True)

ai = json.load(open(P1 / "atlas_index.json"))
am = np.load(P1 / "atlas_means.npy")          # types x genes, mean log-CPM, human infant
types = ai["types"]; genes = ai["genes"]
gi = {g: i for i, g in enumerate(genes)}
ti = {t: i for i, t in enumerate(types)}

RECEPTORS = {
    "IL-6":      ["IL6R", "IL6ST"],
    "IL-1beta":  ["IL1R1", "IL1RAP", "IL1R2"],
    "TNF-alpha": ["TNFRSF1A", "TNFRSF1B"],
    "IL-17A":    ["IL17RA", "IL17RC"],
    "IFN-gamma": ["IFNGR1", "IFNGR2"],
    "IFN-I":     ["IFNAR1", "IFNAR2"],
    "IL-18":     ["IL18R1", "IL18RAP"],
    "IL-4/13":   ["IL4R", "IL13RA1"],
    "TGF-beta":  ["TGFBR1", "TGFBR2"],
    "CSF1":      ["CSF1R"],
    "fractalkine": ["CX3CR1"],
    "complement": ["C3AR1", "C5AR1", "CR1", "ITGAM"],
    "adenosine (Xiong route)": ["ADORA1", "ADORA2A", "ADORA2B"],
}
FOCAL = ["OPC", "COP", "OL"]
CTRL = ["Mic-P2RY12", "Astrocyte"]
NEUR = [t for t in types if t.startswith(("L2", "L4", "L5", "L6", "IN-"))]

# expression threshold: a gene is "present" in a type if mean log-CPM exceeds the
# 40th percentile of that type's nonzero distribution (conservative, type-internal)
thr = {}
for t in types:
    v = am[ti[t]]
    nz = v[v > 0]
    thr[t] = float(np.percentile(nz, 40)) if nz.size else 0.0

rows = []
for cyto, recs in RECEPTORS.items():
    for r in recs:
        if r not in gi:
            rows.append(dict(cytokine=cyto, receptor=r, status="NOT_IN_ATLAS"))
            continue
        j = gi[r]
        e = {t: float(am[ti[t]][j]) for t in types}
        neur = float(np.mean([e[t] for t in NEUR])) if NEUR else float("nan")
        rows.append(dict(
            cytokine=cyto, receptor=r, status="ok",
            OPC=e.get("OPC", float("nan")), COP=e.get("COP", float("nan")),
            OL=e.get("OL", float("nan")),
            Microglia=e.get("Mic-P2RY12", float("nan")),
            Astrocyte=e.get("Astrocyte", float("nan")), Neuron_mean=neur,
            present_OPC=bool(e.get("OPC", 0) > thr["OPC"]),
            present_COP=bool(e.get("COP", 0) > thr.get("COP", 0)),
            present_Mic=bool(e.get("Mic-P2RY12", 0) > thr.get("Mic-P2RY12", 0)),
        ))
json.dump(rows, open(OUT / "cytokine_receptor_expression.json", "w"), indent=1, default=str)

L = ["=" * 104,
     "Cytokine-receptor substrate on human infant PFC oligodendrocyte lineage (GSE294786)",
     "=" * 104,
     "mean log-CPM per cell type. 'present' = above that cell type's own 40th-pct of",
     "nonzero expression. Microglia = positive control for the measurement.",
     "",
     f"{'cytokine':22s} {'receptor':10s} {'OPC':>7s} {'COP':>7s} {'OL':>7s} {'Micro':>7s} "
     f"{'Astro':>7s} {'Neuron':>7s}  OPC COP Mic"]
for r in rows:
    if r["status"] != "ok":
        L.append(f"{r['cytokine']:22s} {r['receptor']:10s}  -- not in atlas --")
        continue
    L.append(f"{r['cytokine']:22s} {r['receptor']:10s} {r['OPC']:>7.3f} {r['COP']:>7.3f} "
             f"{r['OL']:>7.3f} {r['Microglia']:>7.3f} {r['Astrocyte']:>7.3f} "
             f"{r['Neuron_mean']:>7.3f}   {'Y' if r['present_OPC'] else '.'}   "
             f"{'Y' if r['present_COP'] else '.'}   {'Y' if r['present_Mic'] else '.'}")

ok = [r for r in rows if r["status"] == "ok"]
nOPC = sum(1 for r in ok if r["present_OPC"])
nMic = sum(1 for r in ok if r["present_Mic"])
L += ["", f"SUMMARY: {nOPC}/{len(ok)} receptors present on OPC; {nMic}/{len(ok)} on microglia.",
      "If microglial presence >> OPC presence, the canonical route to the OL lineage at this",
      "developmental stage is INDIRECT (via microglia/astrocyte), not direct cytokine-on-OPC.",
      "That does not weaken an immune model; it relocates the decisive experiment."]
txt = "\n".join(L)
(OUT / "CYTOKINE_RECEPTOR_REPORT.txt").write_text(txt + "\n")
print(txt)

#!/usr/bin/env python3
"""
§66 + §86 + §87 — programme trajectories along the TRUE OPC->MOL axis (atlas 2, GSE155488).

Three questions that could only be asked once an atlas with full stage representation existed:

§66  RECEPTOR COMPETENCE ALONG MATURATION.
     Does immune sensitivity itself change with OL maturation? A cytokine may act only on a
     narrow developmental state. The program calls this "a major mechanistic prediction."

§86  LIPID AS OUTPUT BOTTLENECK, NOT ENTRY POINT.
     If lipid biology is the myelinating OL's OUTPUT programme rather than an ASD entry
     point, lipid modules should rise steeply and late along the axis (peaking at MFOL/MOL),
     while ASD-implicated programmes (chromatin/state-control) should NOT track that curve.
     This distinguishes "lipid gene mutation -> ASD" from "state disruption -> failure to
     engage the lipid programme."

§87  OL STATE x IMMUNE RESPONSE covariation: do immune-response genes track immature states,
     stress-like states, or mature-OL failure?

Per-stage values are PSEUDOBULK MEANS with per-sample (donor) breakdown so a result driven
by one donor is visible. No causal claim is made from expression trajectories.
"""
import gzip, json, sys
import numpy as np
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
GSE = COR / "celltype/GSE155488"
OUT = BASE / "outputs/phase10"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = [
    ("GSE155488_All_fetal_OL_lineage_subset.csv.gz", "fetal"),
    ("GSM4875773_HA_773_OL_lineage_subset.csv.gz", "2yr_a"),
    ("GSM4875774_HA_776_OL_lineage_subset.csv.gz", "2yr_b"),
    ("GSM4875775_HA_772_OL_lineage_subset.csv.gz", "13yr"),
    ("GSM4875776_HA_738_OL_lineage_subset.csv.gz", "adult_62"),
    ("GSM4875777_OPK_363_OL_lineage_subset.csv.gz", "adult_58"),
]
STAGE = {
    "OPC":  "PDGFRA CSPG4 PTPRZ1 LHFPL3 SOX6 OLIG1 DSCAM".split(),
    "COP":  "GPR17 BCAS1 ENPP6 TMEM235 NKX2-2 SIRT2 FYN".split(),
    "NFOL": "BCAS1 ENPP6 ITPR2 CNP TCF7L2 NKX6-2".split(),
    "MFOL": "MAL MOG OPALIN CTPS1 CNP PLP1".split(),
    "MOL":  "PLP1 MBP MOBP CLDN11 APOD KLK6 ASPA CARNS1".split(),
}
ORDER = ["OPC", "COP", "NFOL", "MFOL", "MOL"]

PROGRAMS = {
    # --- §66 immune receptors
    "RCPT_IL6_classic": ["IL6R"], "RCPT_IL6_gp130": ["IL6ST"],
    "RCPT_IL1": ["IL1R1", "IL1RAP"], "RCPT_TNFR1": ["TNFRSF1A"],
    "RCPT_TNFR2": ["TNFRSF1B"], "RCPT_IL17": ["IL17RA", "IL17RC"],
    "RCPT_IFNg": ["IFNGR1", "IFNGR2"], "RCPT_IFN_I": ["IFNAR1", "IFNAR2"],
    "RCPT_ADORA1": ["ADORA1"], "RCPT_TGFb": ["TGFBR1", "TGFBR2"],
    "RCPT_complement": ["C3AR1", "C5AR1", "CR1"],
    # --- §86 lipid output programmes
    "LIPID_chol_synthesis": "HMGCR HMGCS1 MVK PMVK MVD IDI1 FDPS FDFT1 SQLE LSS CYP51A1 TM7SF2 MSMO1 NSDHL HSD17B7 EBP SC5D DHCR7 DHCR24 SREBF2 SCAP INSIG1".split(),
    "LIPID_sphingo_VLCFA": "SPTLC1 SPTLC2 CERS2 DEGS1 ELOVL1 ELOVL5 ELOVL7 UGT8 GAL3ST1 FA2H SGMS1 ASAH1".split(),
    "LIPID_plasmalogen_perox": "GNPAT AGPS FAR1 TMEM189 ABCD1 ABCD3 ACOX1 PEX5 PEX7 PEX13 HSD17B4".split(),
    "LIPID_uptake": "LDLR LRP1 NPC1 NPC2 ABCA1 ABCG1 SCARB1 APOE".split(),
    # --- myelin structural output
    "MYELIN_structural": "MBP PLP1 MOG MAG CNP MOBP CLDN11 MAL OPALIN ERMN".split(),
    # --- ASD-implicated control programmes (should NOT track the lipid curve if lipid is output)
    "CHROMATIN_state": "CHD8 CHD7 SMARCA4 SMARCA5 KDM6B EZH2 HDAC1 HDAC2 EP300 CREBBP ARID1B SETD5 KMT2A".split(),
    "SFARI_top": "CHD8 SCN2A SYNGAP1 ADNP ARID1B SHANK3 POGZ DYRK1A GRIN2B TBR1 FOXP1 PTEN TCF4 NRXN1 ANK2".split(),
    # --- §87 stress / immune response states
    "STRESS_ISR": "ATF4 DDIT3 EIF2AK3 ERN1 XBP1 HSPA5 ATF6 TRIB3".split(),
    "IMMUNE_response": "STAT1 STAT3 NFKB1 RELA IRF1 IRF9 B2M HLA-A HLA-B HLA-C CD74 SERPINA3".split(),
    "OXPHOS_mito": "NDUFA1 NDUFB2 COX4I1 COX5A COX7C ATP5F1A ATP5F1B UQCRB SDHB CYCS".split(),
}

# ---------------- load, per-sample
per_sample = {}
gene_sets = []
for fn, lab in SAMPLES:
    p = GSE / fn
    if not p.exists():
        continue
    with gzip.open(p, "rt") as f:
        rows, gn = [], []
        for line in f:
            i = line.find(",")
            gn.append(line[:i].strip().strip('"'))
            rows.append(np.fromstring(line[i + 1:], sep=",", dtype=np.float32))
    M = np.vstack(rows)
    tot = M.sum(axis=0); tot[tot <= 0] = 1
    per_sample[lab] = (np.log1p(M / tot * 1e4).astype(np.float32), gn)
    gene_sets.append(set(gn))
    print(f"[load] {lab}: {M.shape}", file=sys.stderr)

common = sorted(set.intersection(*gene_sets))
ci = {g: i for i, g in enumerate(common)}
print(f"[genes] {len(common)}", file=sys.stderr)

# per-sample stage assignment + per-stage pseudobulk
stage_vals = {p: {s: [] for s in ORDER} for p in PROGRAMS}
stage_n = {s: 0 for s in ORDER}
donor_tab = []
for lab, (X, gn) in per_sample.items():
    idx = np.array([gn.index(g) for g in common])
    Xs = X[idx]
    sc = []
    for s in ORDER:
        ix = [ci[m] for m in STAGE[s] if m in ci]
        v = Xs[ix].mean(axis=0) if ix else np.zeros(Xs.shape[1])
        sc.append((v - v.mean()) / (v.std() + 1e-9))
    S = np.vstack(sc)
    assign = np.argmax(S, axis=0)
    for k, s in enumerate(ORDER):
        m = assign == k
        stage_n[s] += int(m.sum())
        if m.sum() < 10:
            continue
        mu = Xs[:, m].mean(axis=1)
        for pname, genes in PROGRAMS.items():
            ix = [ci[g] for g in genes if g in ci]
            if ix:
                stage_vals[pname][s].append(float(mu[ix].mean()))
        donor_tab.append((lab, s, int(m.sum())))

res = {p: {s: (float(np.mean(v)) if v else None, len(v)) for s, v in d.items()}
       for p, d in stage_vals.items()}
json.dump(dict(programs=res, stage_n=stage_n, donors=donor_tab),
          open(OUT / "stage_trajectories.json", "w"), indent=1)

def norm(vals):
    v = [x for x in vals if x is not None]
    if not v or max(v) == min(v):
        return ["  -  "] * len(vals)
    lo, hi = min(v), max(v)
    return [f"{(x-lo)/(hi-lo):5.2f}" if x is not None else "  -  " for x in vals]

L = ["=" * 100, "§66/§86/§87 PROGRAMME TRAJECTORIES ALONG THE OPC->MOL AXIS (atlas 2)",
     "=" * 100,
     "values are min-max normalised across the 5 stages (1.00 = that programme's peak stage);",
     "each stage averaged over donors contributing >=10 cells to it.",
     f"cells per stage: " + "  ".join(f"{s}={stage_n[s]}" for s in ORDER), "",
     f"{'programme':26s} " + "  ".join(f"{s:>5s}" for s in ORDER) + "   peak"]
for pname in PROGRAMS:
    vals = [res[pname][s][0] for s in ORDER]
    nn = norm(vals)
    pk = ORDER[int(np.nanargmax([v if v is not None else -9 for v in vals]))]
    L.append(f"{pname:26s} " + "  ".join(nn) + f"   {pk}")
L += ["", "DONOR CONTRIBUTIONS per stage (guards against single-donor artefacts):"]
for lab, s, n in donor_tab:
    L.append(f"   {lab:10s} {s:5s} n={n}")
txt = "\n".join(L)
(OUT / "STAGE_TRAJECTORY_REPORT.txt").write_text(txt + "\n")
print(txt)

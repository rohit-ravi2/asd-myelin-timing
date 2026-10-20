#!/usr/bin/env python3
"""
Conditional gate competition  (program sections IX, XV, XVI, XVII)

THE QUESTION THIS ANSWERS
-------------------------
Not "is module X enriched?" but:
    does module X carry INDEPENDENT ASD signal once gene properties AND the broader
    developmental programmes it sits inside are accounted for?

Key nested contrasts, pre-specified:
  C1  lipid modules | OL identity        -> is lipid biology anything beyond "expressed in OLs"?
  C2  OL-lineage    | chromatin/state    -> is OL signal anything beyond general state-control?
  C3  myelin output | OL maturation      -> does signal reach the structural end?
  C4  mitochondrial | OL identity        -> is mito anything beyond high-energy OL programmes?
  C5  each gate     | all other gates    -> full joint model

MODEL
  Z_g ~ module indicators + LOEUF + log(CDS) + log(brain expr) + log(NSNPs) + OL specificity
Marginal (each alone) and joint (all together) are both reported; the difference IS the result.

MATCHED NULL: for every module, 1000 gene sets matched on LOEUF decile x CDS decile x
expression decile; the module's joint coefficient is re-estimated each time. This is what
stops big constrained neurodevelopmental gene sets from manufacturing signal.
"""
import gzip, json, os, re, sys
import numpy as np
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
SCOPE = Path("data/dvad-corpus/asd_scope")
OUT = BASE / "outputs/phase2"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20260917)
N_MATCH = 1000

# ---------------------------------------------------------------- gene sets
gmt = {}
for fn in ["c2.cp.v2024.1.Hs.symbols.gmt", "c5.go.bp.v2024.1.Hs.symbols.gmt",
           "c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    p = COR / "pathways" / fn
    if p.exists():
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if len(f) > 2:
                gmt[f[0]] = set(f[2:])


def U(*pats, exclude=()):
    g = set()
    for name, mem in gmt.items():
        if any(re.search(p, name, re.I) for p in pats) and \
           not any(re.search(e, name, re.I) for e in exclude):
            g |= mem
    return g


MOD = {}
# --- gate 1 chromatin / developmental state control
MOD["G1_CHROMATIN_STATE"] = U(r"CHROMATIN_REMODEL", r"CHROMATIN_ORGANIZATION",
                              r"HISTONE_MODIFICATION", r"COVALENT_CHROMATIN")
# --- gate 2 neurogenic -> gliogenic fate switch (curated)
MOD["G2_FATE_SWITCH"] = set("""NFIA NFIB NFIX SOX9 SOX8 SOX10 HES1 HES5 HEY1 HEY2 NOTCH1
NOTCH2 RBPJ DLL1 JAG1 STAT3 JAK1 JAK2 EGFR FGFR1 FGFR2 FGFR3 OLIG1 OLIG2 ASCL1 NEUROG1
NEUROG2 NEUROD1 NEUROD2 EOMES PAX6 GLI3 ZBTB20 TCF7L2 RFX4 LHX2 POU3F2 POU3F3 SALL3 ID2 ID4
CNTFR LIF IL6ST BMP2 BMP4 BMPR1A SMAD1 SMAD4""".split())
# --- gate 3 OL lineage, split by stage
MOD["G3a_OPC_PROLIF"] = set("""PDGFRA CSPG4 SOX2 ASCL1 OLIG1 OLIG2 CCND1 MKI67 NKX2-2 HES5
ID2 ID4 PTPRZ1 LHFPL3 GPR17 DSCAM BCAN NEU4 SHH GLI1 GLI2 CHRDL1""".split())
MOD["G3b_OPC_DIFF_TIMING"] = set("""MYRF SOX10 ZNF488 NKX6-2 TCF7L2 YY1 ZFP191 ZNF191 HDAC1
HDAC2 CHD7 CHD8 SMARCA4 SMARCA5 KDM6B PRDM8 SIRT2 EP300 CREBBP DNMT1 EZH2 ARX NFATC1
GPR17 ENPP6 BMP4 WNT3A CTNNB1 GSK3B""".split())
MOD["G3c_OL_MATURATION"] = set("""CNP MAG MBP PLP1 MOG MOBP CLDN11 OPALIN ERMN ASPA MAL
TF CA2 UGT8 GAL3ST1 FA2H SIRT2 PLLP QKI CNTN2 NFASC NRCAM""".split())
# --- gate 4 mature myelin structural output only
MOD["G4_MYELIN_STRUCTURAL"] = set("PLP1 MBP MOG MAG CNP MOBP CLDN11 MAL PMPa22 PMP22 MPZ OPALIN ERMN".split())
# --- gate 5 lipid, partitioned (section XIV)
MOD["G5a_CHOL_SYNTHESIS"] = set("""SREBF2 SCAP INSIG1 INSIG2 HMGCR HMGCS1 MVK PMVK MVD IDI1
FDPS FDFT1 SQLE LSS CYP51A1 TM7SF2 MSMO1 NSDHL HSD17B7 EBP SC5D DHCR7 DHCR24""".split())
MOD["G5b_CHOL_TRAFFICKING"] = set("""NPC1 NPC2 ABCA1 ABCG1 ABCA2 LDLR LRP1 LRP8 VLDLR APOE
APOA1 SCARB1 STARD3 OSBP OSBPL1A SOAT1 CYP46A1 SORT1 PLTP CETP""".split())
MOD["G5c_SPHINGO_VLCFA"] = set("""SPTLC1 SPTLC2 SPTLC3 CERS1 CERS2 CERS3 CERS4 CERS5 CERS6
DEGS1 DEGS2 ELOVL1 ELOVL2 ELOVL3 ELOVL4 ELOVL5 ELOVL6 ELOVL7 UGT8 GAL3ST1 FA2H SGMS1 SGMS2
SMPD1 ASAH1 GBA ARSA CERT1 COL4A3BP STARD11""".split())
MOD["G5d_PEROXISOME_PLASMALOGEN"] = set("""GNPAT AGPS FAR1 FAR2 TMEM189 ABCD1 ABCD2 ABCD3
ACOX1 PEX1 PEX2 PEX3 PEX5 PEX6 PEX7 PEX10 PEX12 PEX13 PEX14 PEX16 PEX19 PEX26 HSD17B4""".split())
# --- gate 6 mitochondrial (MitoCarta)
try:
    import subprocess
    mc = subprocess.run([os.environ.get("ML_PYTHON", sys.executable), "-c",
        "import pandas as pd;d=pd.ExcelFile('%s').parse('A Human MitoCarta3.0');"
        "print('\\n'.join(map(str,d['Symbol'].dropna().tolist())))"
        % (SCOPE/'mito_metabolic/Human.MitoCarta3.0.xls')],
        capture_output=True, text=True, timeout=300)
    MOD["G6_MITOCARTA"] = {x.strip() for x in mc.stdout.split("\n") if x.strip()}
except Exception as e:
    print(f"[warn] MitoCarta parse failed: {e}", file=sys.stderr)
    MOD["G6_MITOCARTA"] = set()
# --- gate 7 signalling (section XVII)
MOD["G7a_WNT"] = U(r"^WP_WNT", r"REACTOME_SIGNALING_BY_WNT", r"GOBP_WNT_SIGNALING_PATHWAY$")
MOD["G7b_ERK_MAPK"] = U(r"REACTOME_MAPK", r"GOBP_ERK1_AND_ERK2_CASCADE$", r"^WP_MAPK")
MOD["G7c_MTOR_PI3K"] = U(r"MTOR", r"PI3K_AKT", r"REACTOME_PI3K")
MOD["G7d_ER_STRESS_ISR"] = U(r"ENDOPLASMIC_RETICULUM_STRESS", r"UNFOLDED_PROTEIN",
                             r"INTEGRATED_STRESS")
MOD["G7e_NOTCH"] = U(r"REACTOME_SIGNALING_BY_NOTCH$", r"GOBP_NOTCH_SIGNALING_PATHWAY$")
MOD["G7f_THYROID"] = U(r"THYROID_HORMONE")
# --- gate 8 iron / metal (section XXXII)
MOD["G8_IRON"] = set("""TF TFRC TFR2 SLC11A2 SLC40A1 FTH1 FTL CP HEPH STEAP3 ACO1 IREB2
HAMP SLC39A14 SLC39A8 FXN ABCB7 MFRN1 SLC25A37 SLC25A28""".split())
# --- gate 9 neuronal / synaptic / interneuron
MOD["G9a_SYNAPTIC"] = U(r"GOBP_SYNAP", r"GOCC_POSTSYNAPTIC_DENSITY", r"REACTOME_NEURONAL_SYSTEM")
MOD["G9b_INTERNEURON_EI"] = set("""GAD1 GAD2 SLC32A1 PVALB SST VIP NPY CCK CALB1 CALB2 LHX6
DLX1 DLX2 DLX5 DLX6 NKX2-1 ERBB4 NRG1 NRG3 GABRA1 GABRB2 GABRG2 KCNC1 KCNC2 SCN1A SCN1B""".split())
# --- gate 10 activity-dependent myelination (section XVIII)
MOD["G10_ACTIVITY_MYELIN"] = set("""NRG1 NRG3 ERBB3 ERBB4 BDNF NTRK2 GRIA1 GRIA2 GRIA3 GRIA4
GRIN1 GRIN2A GRIN2B GRIN3A P2RY1 P2RX7 ADORA1 ADORA2A GPR17 CX3CL1 CX3CR1 ATP1A2 KCNJ10""".split())
# --- gate 11 immune / microglial niche
MOD["G11_IMMUNE_NICHE"] = set("""TNF TNFRSF1A TNFRSF1B IL1B IL1R1 IL6 IL6R IL6ST IL17A IL17RA
IFNG IFNGR1 IFNGR2 STAT1 STAT3 NFKB1 RELA C1QA C1QB C1QC C3 C3AR1 CR3 ITGAM CX3CR1 TREM2
TYROBP P2RY12 CSF1R IRF8 SPI1 CCL2 CCL3 CCR5 TLR2 TLR4 MYD88""".split())

# --- ASD reference sets
sf_all, sf1 = set(), set()
import csv
with open(SCOPE / "sfari/SFARI-Gene_genes_2026-09-10.csv") as f:
    for r in csv.DictReader(f):
        g = r.get("gene-symbol") or r.get("gene_symbol")
        s = (r.get("gene-score") or r.get("gene_score") or "").strip()
        if g:
            sf_all.add(g)
            if s == "1":
                sf1.add(g)
MOD["REF_SFARI_ALL"] = sf_all
MOD["REF_SFARI_S1"] = sf1
MOD["NEG_OLFACTORY"] = U(r"OLFACTORY")
MOD["NEG_KERATIN"] = U(r"KERATIN", r"CORNIFIED")

for k in list(MOD):
    if len(MOD[k]) < 8:
        print(f"[drop] {k}: only {len(MOD[k])} genes", file=sys.stderr)
        del MOD[k]
print(f"[modules] {len(MOD)}", file=sys.stderr)
for k in sorted(MOD):
    print(f"   {k:28s} {len(MOD[k])}", file=sys.stderr)

# ---------------------------------------------------------------- gene-level data
gloc = {}
for line in open(COR / "tools/magma_geneloc/NCBI37.3.gene.loc"):
    p = line.split()
    if len(p) >= 6:
        gloc[p[0]] = p[5]
Z, NS = {}, {}
with open(BASE / "outputs/magma/grove.w35.genes.out") as f:
    f.readline()
    for line in f:
        p = line.split()
        s = gloc.get(p[0])
        if s:
            Z[s] = float(p[7]); NS[s] = float(p[4])
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
ai = json.load(open(BASE / "outputs/phase1/atlas_index.json"))
am = np.load(BASE / "outputs/phase1/atlas_means.npy")
atg = {g: i for i, g in enumerate(ai["genes"])}
types = ai["types"]
OLROWS = [i for i, t in enumerate(types) if t in ("OPC", "COP", "OL")]
allmean = am.mean(axis=0)
ol_spec_v = np.log2((am[OLROWS].mean(axis=0) + 1e-3) / (allmean + 1e-3))
expr_v = allmean

genes = [g for g in Z if g in loeuf and g in atg and cds.get(g, 0) > 0
         and expr_v[atg[g]] > 0]
print(f"[genes] {len(genes)} usable", file=sys.stderr)
idx = np.array([atg[g] for g in genes])
y = np.array([Z[g] for g in genes])
COV = np.column_stack([
    [loeuf[g] for g in genes],
    np.log([cds[g] for g in genes]),
    np.log(expr_v[idx] + 1e-6),
    np.log([NS[g] + 1 for g in genes]),
    ol_spec_v[idx],
])
COVNAMES = ["LOEUF", "logCDS", "logExpr", "logNSNP", "OLspecificity"]
M = {k: np.array([1.0 if g in v else 0.0 for g in genes]) for k, v in MOD.items()}
for k in list(M):
    if M[k].sum() < 8:
        del M[k]


def ols(y, X):
    X = np.column_stack([np.ones(len(y)), X])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    df = len(y) - X.shape[1]
    s2 = (r @ r) / df
    se = np.sqrt(np.diag(s2 * np.linalg.pinv(X.T @ X)))
    t = b / se
    from math import erf, sqrt
    p = np.array([2 * (1 - 0.5 * (1 + erf(abs(x) / sqrt(2)))) for x in t])
    return b, se, t, p


keys = sorted(M)
# marginal
marg = {}
for k in keys:
    b, se, t, p = ols(y, np.column_stack([M[k], COV]))
    marg[k] = dict(beta=b[1], se=se[1], t=t[1], p=p[1], n=int(M[k].sum()))
# joint
Xj = np.column_stack([M[k] for k in keys] + [COV])
b, se, t, p = ols(y, Xj)
joint = {k: dict(beta=b[i + 1], se=se[i + 1], t=t[i + 1], p=p[i + 1])
         for i, k in enumerate(keys)}
covres = {COVNAMES[j]: dict(beta=b[len(keys) + 1 + j], t=t[len(keys) + 1 + j],
                            p=p[len(keys) + 1 + j]) for j in range(len(COVNAMES))}

# matched null on joint coefficient
lo_a = np.array([loeuf[g] for g in genes]); cd_a = np.array([cds[g] for g in genes])
ex_a = expr_v[idx]
def dec(v): return np.digitize(v, np.quantile(v, np.linspace(0, 1, 6)[1:-1]))
strata = defaultdict(list)
key3 = list(zip(dec(lo_a), dec(cd_a), dec(ex_a)))
for i, k3 in enumerate(key3):
    strata[k3].append(i)
print(f"[perm] matched null, {N_MATCH} draws", file=sys.stderr)
matp = {}
for k in keys:
    obs = abs(joint[k]["t"])
    memb = np.where(M[k] > 0)[0]
    ge = 0
    for _ in range(N_MATCH):
        pick = np.array([strata[key3[i]][rng.integers(len(strata[key3[i]]))] for i in memb])
        mm = np.zeros(len(genes)); mm[pick] = 1.0
        Xk = np.column_stack([mm] + [M[o] for o in keys if o != k] + [COV])
        bb, ss, tt, pp = ols(y, Xk)
        if abs(tt[1]) >= obs:
            ge += 1
    matp[k] = (ge + 1) / (N_MATCH + 1)

# nested contrasts
def contrast(focal, condition_on, label):
    b1, s1, t1, p1 = ols(y, np.column_stack([M[focal], COV]))
    cols = [M[focal]] + [M[c] for c in condition_on if c in M] + [COV]
    b2, s2, t2, p2 = ols(y, np.column_stack(cols))
    return dict(label=label, focal=focal, conditioned_on=condition_on,
                beta_marginal=b1[1], p_marginal=p1[1],
                beta_conditional=b2[1], p_conditional=p2[1],
                attenuation=(1 - b2[1] / b1[1]) if b1[1] != 0 else None)

OLSET = ["G3a_OPC_PROLIF", "G3b_OPC_DIFF_TIMING", "G3c_OL_MATURATION"]
contrasts = [
    contrast("G5a_CHOL_SYNTHESIS", OLSET, "C1 cholesterol synthesis | OL lineage"),
    contrast("G5b_CHOL_TRAFFICKING", OLSET, "C1 cholesterol trafficking | OL lineage"),
    contrast("G5c_SPHINGO_VLCFA", OLSET, "C1 sphingolipid/VLCFA | OL lineage"),
    contrast("G5d_PEROXISOME_PLASMALOGEN", OLSET, "C1 peroxisome/plasmalogen | OL lineage"),
    contrast("G3a_OPC_PROLIF", ["G1_CHROMATIN_STATE"], "C2 OPC prolif | chromatin"),
    contrast("G3b_OPC_DIFF_TIMING", ["G1_CHROMATIN_STATE"], "C2 OPC diff-timing | chromatin"),
    contrast("G4_MYELIN_STRUCTURAL", ["G3c_OL_MATURATION"], "C3 myelin output | OL maturation"),
    contrast("G6_MITOCARTA", OLSET, "C4 mitochondrial | OL lineage"),
    contrast("G2_FATE_SWITCH", ["G1_CHROMATIN_STATE"], "C5 fate-switch | chromatin"),
]

json.dump(dict(marginal=marg, joint=joint, matched_p=matp, covariates=covres,
               contrasts=contrasts, n_genes=len(genes)),
          open(OUT / "gate_competition.json", "w"), indent=1, default=float)

L = ["=" * 110, "CONDITIONAL GATE COMPETITION — ASD common-variant gene-level Z", "=" * 110,
     f"genes={len(genes)}  covariates: {', '.join(COVNAMES)}",
     "marginal = module + covariates.   joint = ALL modules together + covariates.",
     "matchP = 1000 LOEUF x CDS x expression matched sets, joint coefficient re-estimated.", "",
     f"{'module':28s} {'n':>5s} {'beta_marg':>10s} {'p_marg':>9s} {'beta_joint':>11s} {'p_joint':>9s} {'matchP':>7s}"]
for k in sorted(keys, key=lambda x: marg[x]["p"]):
    L.append(f"{k:28s} {marg[k]['n']:>5d} {marg[k]['beta']:>+10.4f} {marg[k]['p']:>9.3g} "
             f"{joint[k]['beta']:>+11.4f} {joint[k]['p']:>9.3g} {matp[k]:>7.3f}")
L += ["", "covariates in joint model:"]
for c, v in covres.items():
    L.append(f"   {c:16s} beta={v['beta']:+.4f} t={v['t']:+.2f} p={v['p']:.3g}")
L += ["", "=" * 110, "PRE-SPECIFIED NESTED CONTRASTS", "=" * 110,
      f"{'contrast':46s} {'b_marg':>9s} {'p_marg':>9s} {'b_cond':>9s} {'p_cond':>9s} {'atten':>7s}"]
for c in contrasts:
    at = f"{c['attenuation']*100:.0f}%" if c["attenuation"] is not None else "NA"
    L.append(f"{c['label']:46s} {c['beta_marginal']:>+9.4f} {c['p_marginal']:>9.3g} "
             f"{c['beta_conditional']:>+9.4f} {c['p_conditional']:>9.3g} {at:>7s}")
txt = "\n".join(L)
(OUT / "GATE_COMPETITION_REPORT.txt").write_text(txt + "\n")
print(txt)

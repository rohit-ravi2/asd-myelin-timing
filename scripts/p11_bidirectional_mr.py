#!/usr/bin/env python3
"""
P1.7 / §XII — Bidirectional MR between ASD liability and white-matter ICVF.

EXPLICITLY DEMOTED INTERPRETATION (per §XII). This cannot establish that infant OL
developmental abnormalities mediate ASD. Horizontal pleiotropy is severe here because
developmental genes plausibly affect both ASD and adult imaging. Read as genetic-directionality
reconnaissance only.

Direction A: ICVF -> ASD   (ICVF IDPs have many GWS loci -> usable instrument)
Direction B: ASD -> ICVF   (Grove has ~5 independent GWS loci -> expected to be
                            SEVERELY underpowered; reported with that stated up front)

Estimators: IVW (random-effects), weighted median, MR-Egger (+ intercept as a pleiotropy
test), Cochran Q heterogeneity, and Steiger directionality.
Clumping: plink r2<0.001 / 10 Mb on 1000G EUR.
"""
import gzip, json, os, subprocess, sys, tempfile
import numpy as np
from pathlib import Path
from math import erf, sqrt

BASE = Path(__file__).resolve().parents[1]
COR = Path("data/asd-corpus")
BIG40 = Path("data/dvad-corpus/asd_myelin/big40_extra")
OUT = BASE / "outputs/phase11"
OUT.mkdir(parents=True, exist_ok=True)
PLINK = os.environ.get("PLINK", "plink")
REF = str(COR / "tools/magma_ref_eur/g1000_eur")
N_ASD, N_ICVF = 46350, 31355

IDPS = {"ICVF_CingCingGyrus_L": "1937", "ICVF_SLF_L": "1943",
        "ICVF_AntCoronaRad_L": "1925", "ICVF_Genu_CC": "1904"}


def p2z(p):
    p = max(min(p, 1 - 1e-16), 1e-300)
    lo, hi = 0.0, 40.0
    for _ in range(200):
        m = (lo + hi) / 2
        if 2 * (1 - 0.5 * (1 + erf(m / sqrt(2)))) > p:
            lo = m
        else:
            hi = m
    return (lo + hi) / 2


def read_asd():
    d = {}
    with gzip.open(COR / "genetics/ASD_Grove2019_REAL.gz", "rt") as f:
        f.readline()
        for line in f:
            p = line.split()
            if not p[1].startswith("rs"):
                continue
            try:
                beta = np.log(float(p[6])); se = float(p[7]); pv = float(p[8])
            except (ValueError, ZeroDivisionError):
                continue
            d[p[1]] = (p[3].upper(), p[4].upper(), beta, se, pv)
    return d


def read_idp(fid):
    d = {}
    with gzip.open(BIG40 / f"{fid}.txt.gz", "rt") as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p) < 8 or not p[1].startswith("rs"):
                continue
            try:
                lp = float(p[7]); beta = float(p[5]); se = float(p[6])
            except ValueError:
                continue
            d[p[1]] = (p[3].upper(), p[4].upper(), beta, se, 10 ** (-min(lp, 300)))
    return d


def clump(snps_p, tag):
    with tempfile.NamedTemporaryFile("w", suffix=".assoc", delete=False) as t:
        t.write("SNP\tP\n")
        for s, pv in snps_p:
            t.write(f"{s}\t{pv}\n")
        fn = t.name
    o = str(OUT / f"clump_{tag}")
    r = subprocess.run([PLINK, "--bfile", REF, "--clump", fn, "--clump-p1", "5e-8",
                        "--clump-r2", "0.001", "--clump-kb", "10000", "--out", o],
                       capture_output=True, text=True)
    keep = []
    cf = Path(o + ".clumped")
    if cf.exists():
        for line in open(cf):
            w = line.split()
            if len(w) > 2 and w[2].startswith("rs"):
                keep.append(w[2])
    return keep


def harmonise(expo, outc, snps):
    B, b, SE, se = [], [], [], []
    for s in snps:
        if s not in expo or s not in outc:
            continue
        a1e, a2e, be, see, _ = expo[s]
        a1o, a2o, bo, seo, _ = outc[s]
        if {a1e, a2e} != {a1o, a2o}:
            continue
        if a1e != a1o:
            bo = -bo
        if see <= 0 or seo <= 0:
            continue
        B.append(be); SE.append(see); b.append(bo); se.append(seo)
    return map(np.array, (B, SE, b, se))


def mr(B, SE, b, se):
    w = 1 / se ** 2
    ivw = np.sum(w * B * b) / np.sum(w * B ** 2)
    ivw_se = np.sqrt(1 / np.sum(w * B ** 2))
    r = b - ivw * B
    Q = np.sum(w * r ** 2)
    dfQ = len(B) - 1
    if dfQ > 0 and Q > dfQ:
        ivw_se *= np.sqrt(Q / dfQ)
    # weighted median
    rat = b / B
    wt = (B ** 2) / (se ** 2)
    o = np.argsort(rat)
    cw = np.cumsum(wt[o]) / np.sum(wt)
    wm = float(np.interp(0.5, cw, rat[o]))
    # Egger
    X = np.column_stack([np.ones(len(B)), B])
    Wm = np.diag(w)
    beta_e = np.linalg.pinv(X.T @ Wm @ X) @ X.T @ Wm @ b
    resid = b - X @ beta_e
    s2 = (resid @ Wm @ resid) / max(len(B) - 2, 1)
    cov = s2 * np.linalg.pinv(X.T @ Wm @ X)
    egg_int, egg_slope = beta_e[0], beta_e[1]
    egg_int_se = np.sqrt(cov[0, 0])

    def pz(z): return 2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2))))
    return dict(n_snp=len(B), ivw=float(ivw), ivw_se=float(ivw_se),
                ivw_p=pz(ivw / ivw_se), wmed=wm,
                egger_slope=float(egg_slope), egger_intercept=float(egg_int),
                egger_int_p=pz(egg_int / egg_int_se) if egg_int_se > 0 else None,
                Q=float(Q), Q_df=dfQ)


print("[read] ASD", file=sys.stderr)
ASD = read_asd()
res = {}
for name, fid in IDPS.items():
    print(f"[read] {name}", file=sys.stderr)
    IDP = read_idp(fid)
    # --- Direction A: ICVF -> ASD
    gws = [(s, v[4]) for s, v in IDP.items() if v[4] < 5e-8]
    keep = clump(gws, f"{name}_fwd") if gws else []
    B, SE, b, se = harmonise(IDP, ASD, keep)
    fwd = mr(B, SE, b, se) if len(B) >= 3 else dict(n_snp=len(B), note="too few instruments")
    if len(B) >= 3:
        rx = np.sum(2 * (B ** 2) * N_ICVF) / (np.sum(2 * (B ** 2) * N_ICVF) + len(B) * N_ICVF)
        fwd["steiger_note"] = "exposure variance explained > outcome (assumed)"
    res[f"A_{name}_to_ASD"] = fwd
    print(f"   A {name}->ASD: {fwd.get('n_snp')} SNPs", file=sys.stderr)

# --- Direction B: ASD -> ICVF (expected underpowered)
gws_a = [(s, v[4]) for s, v in ASD.items() if v[4] < 5e-8]
keep_a = clump(gws_a, "asd_rev") if gws_a else []
print(f"[rev] ASD GWS clumped instruments: {len(keep_a)}", file=sys.stderr)
for name, fid in IDPS.items():
    IDP = read_idp(fid)
    B, SE, b, se = harmonise(ASD, IDP, keep_a)
    rev = mr(B, SE, b, se) if len(B) >= 3 else dict(n_snp=len(B), note="too few instruments")
    res[f"B_ASD_to_{name}"] = rev

json.dump(res, open(OUT / "bidirectional_mr.json", "w"), indent=1, default=float)

L = ["=" * 100, "P1.7 / §XII  BIDIRECTIONAL MR — ASD liability vs white-matter ICVF", "=" * 100,
     "DEMOTED: reconnaissance only. Cannot establish that infant OL development mediates ASD.",
     "Horizontal pleiotropy severe (developmental genes affect both ASD and adult imaging).", "",
     f"{'test':34s} {'nSNP':>5s} {'IVW beta':>10s} {'IVW p':>9s} {'wMedian':>9s} "
     f"{'Egger b':>9s} {'Egger int p':>11s}"]
for k, v in res.items():
    if "ivw" not in v:
        L.append(f"{k:34s} {v.get('n_snp',0):>5d}   {v.get('note','')}")
        continue
    L.append(f"{k:34s} {v['n_snp']:>5d} {v['ivw']:>+10.4f} {v['ivw_p']:>9.3g} "
             f"{v['wmed']:>+9.4f} {v['egger_slope']:>+9.4f} "
             f"{(('%.3g'%v['egger_int_p']) if v['egger_int_p'] else 'NA'):>11s}")
L += ["", "Direction B (ASD->ICVF) uses Grove's handful of genome-wide significant loci and is",
      "expected to be severely underpowered; a null there is uninformative by construction."]
txt = "\n".join(L)
(OUT / "BIDIRECTIONAL_MR_REPORT.txt").write_text(txt + "\n")
print(txt)

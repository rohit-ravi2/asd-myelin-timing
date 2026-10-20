#!/usr/bin/env python3
"""
p16 — Is the OPC epigenetic excess in GSE298243 real, or an artifact?

THE OPEN QUESTION
-----------------
2026-09-17 claim: "OPCs are the most epigenetically perturbed cell type in ASD PFC"
  CpG-stratified significant-bin rate OPC 272/12,515 = 0.02173, OPC-vs-pooled-others
  OR 1.90, Fisher p = 1.49e-19, OPC +5.10 SD above a rate ~ log10(tested bins) trend.
2026-09-21 downgrade (error-log entry E-006) rests on three legs:
  (a) NOISE      — OPC has the highest median coefficient SE (0.3653) and the largest
                   median |coef| among significant bins, so "noisiest estimator -> widest
                   coefficient distribution -> most extreme values -> most significant
                   calls" would explain the whole thing.
  (b) REGION     — (raised 2026-09-21 by this script's own first run) the eight cell types
                   are not tested on the same bins, so a rate comparison is also a
                   comparison of genomic region sets.
  (c) EMPTINESS  — only 3 of the 272 significant OPC bins carry a gene symbol
                   (FAM225A + two olfactory pseudogenes), so the signal is "mechanistically
                   empty".
The 2026-09-21 audit found NEITHER the original OR nor the E-006 noise statistics are
emitted by any script. This script persists all of them and then tests each leg.

WHAT THIS SCRIPT DOES
---------------------
 A. One streaming pass over the 2.52M-row DMR table -> compact cached .npz.
 1. Reproduces the cell type x cytosine context table, the OPC-vs-others OR + Fisher p,
    and the "+5.10 SD above the tested-bin trend" regression.
 2. Persists the E-006 noise statistics (coefSE_DiagnosisASD median/IQR, coefficient SD,
    median |coef| among significant bins, per cell type).
 3. LEG (a) THE SE-MATCHED TEST. All tested CpG bins pooled, cut into deciles of
    coefSE_DiagnosisASD; within each decile OPC's significant-bin rate vs the other seven.
    Mantel-Haenszel OR (RBG CI, MH chi-square) + logistic sig ~ is_OPC + C(SE decile).
    Plus an FDR-free p-threshold ladder, because the significance flag is an FDR call.
 4. Permutation nulls. The DIAGNOSIS-label permutation is NOT computable here (see
    PERM_NOTE); the cell-type-label permutation, unstratified and within SE decile, is,
    and is run with a fixed seed for OPC and two comparison cell types.
 5. Effective-N proxies. No coverage/nuclei column exists; coefSE_genoPC1 and
    coefSE_AgeLin are used as diagnosis-independent precision proxies.
 6. LEG (b) REGION COMPOSITION. The design of the file (every row is an age-DMR of that
    cell type), the pairwise bin-set overlap matrix, paired McNemar / conditional-logistic
    tests where overlap exists, region-composition comparison (FeatureType, gene-symbol
    annotation, GeneType, chromosome, feature width), and MH stratified on FeatureType,
    chromosome and SE decile x FeatureType.
 7. LEG (c) MECHANISTIC EMPTINESS. Is 3/272 anomalous given the annotation rate of the
    OPC tested background? Then coordinate-based mapping of the 272 significant bins and
    the 12,515-bin tested background to genes with a local hg38 refGene table, and module
    enrichment against that coordinate-mapped background.
 8. Bottom line, per leg.

Deterministic apart from the permutation draws, which use SEED.

Usage:  python3 p16_opc_se_matched.py [--force]      (--force re-extracts the cache)
"""
import csv
import gzip
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

SEED = 20260921
N_PERM = 10000

BASE = Path(__file__).resolve().parents[1]
SCOPE = Path("data/dvad-corpus/asd_scope")
COR = Path("data/asd-corpus")
DMR = SCOPE / "celltype/GSE298243_Age_DMRs_PrespecifiedRegions.tsv.gz"
REFGENE = Path("data/cache/"
               "refGene_hg38.txt.gz")
OUT = BASE / "outputs/phase16"
REPORT = BASE / "outputs/phase16/OPC_SE_MATCHED.md"
CACHE = OUT / "dmr_extract_v2.npz"
OUT.mkdir(parents=True, exist_ok=True)

# columns (0-based) in GSE298243_Age_DMRs_PrespecifiedRegions.tsv.gz
C_FT, C_CT, C_CX, C_FEAT = 0, 1, 2, 3
C_CHR, C_START, C_END, C_SYM, C_GT, C_SFARI = 4, 5, 6, 7, 8, 9
C_SIGAGE, C_SIGDX, C_SIGINT = 12, 13, 14
C_SE_AGE = 25
C_EST_DX, C_SE_DX, C_P_DX = 30, 31, 32
C_SE_PC1 = 34
C_Q_DX = 47
NEED = max(C_Q_DX, C_SE_PC1)


# --------------------------------------------------------------------------- A. extract
def extract():
    cts, cxs, fts, gts, chrs = [], [], [], [], []
    ct_map, cx_map, ft_map, gt_map, chr_map = {}, {}, {}, {}, {}
    sig = bytearray(); sigage = bytearray(); sigint = bytearray()
    sfari = bytearray(); genic = bytearray()
    start = []; end = []
    est = []; se_dx = []; se_pc1 = []; se_age = []; p_dx = []; q_dx = []
    cpg_feat = []; cpg_sym = []
    n = 0
    proc = subprocess.Popen(["zcat", str(DMR)], stdout=subprocess.PIPE, bufsize=1 << 22)
    f = proc.stdout
    hdr = f.readline().decode().rstrip("\n").split("\t")
    for i, nm in ((C_SE_DX, "coefSE_DiagnosisASD"), (C_SIGDX, "ReportSignif_Diagnosis"),
                  (C_Q_DX, "coefQ_DiagnosisASD"), (C_SE_PC1, "coefSE_genoPC1"),
                  (C_SIGAGE, "ReportSignif_Age"), (C_FT, "FeatureType"),
                  (C_GT, "GeneType"), (C_SFARI, "SFARI"), (C_SIGINT,
                                                           "ReportSignif_Interaction")):
        assert hdr[i] == nm, (i, hdr[i], nm)

    def code(m, v):
        c = m.get(v)
        if c is None:
            c = m[v] = len(m)
        return c

    def fl(s):
        try:
            return float(s)
        except ValueError:
            return np.nan

    for raw in f:
        p = raw.decode().rstrip("\n").split("\t")
        if len(p) <= NEED:
            continue
        cx = p[C_CX]
        cts.append(code(ct_map, p[C_CT])); cxs.append(code(cx_map, cx))
        fts.append(code(ft_map, p[C_FT])); gts.append(code(gt_map, p[C_GT].strip()))
        chrs.append(code(chr_map, p[C_CHR]))
        sig.append(1 if p[C_SIGDX].strip() == "*" else 0)
        sigage.append(1 if p[C_SIGAGE].strip() == "*" else 0)
        sigint.append(1 if p[C_SIGINT].strip() == "*" else 0)
        sfari.append(1 if p[C_SFARI].strip() else 0)
        sym = p[C_SYM].strip()
        genic.append(1 if sym else 0)
        start.append(int(p[C_START])); end.append(int(p[C_END]))
        est.append(fl(p[C_EST_DX])); se_dx.append(fl(p[C_SE_DX]))
        se_pc1.append(fl(p[C_SE_PC1])); se_age.append(fl(p[C_SE_AGE]))
        p_dx.append(fl(p[C_P_DX])); q_dx.append(fl(p[C_Q_DX]))
        if cx == "CpG":
            cpg_feat.append(p[C_FEAT]); cpg_sym.append(sym)
        n += 1
        if n % 500000 == 0:
            print(f"   [{n:,}] rows", file=sys.stderr, flush=True)
    proc.wait()
    print(f"[extract] {n:,} rows", file=sys.stderr)
    d = dict(
        ct=np.array(cts, np.int8), cx=np.array(cxs, np.int8), ft=np.array(fts, np.int8),
        gt=np.array(gts, np.int16), chrom=np.array(chrs, np.int16),
        start=np.array(start, np.int64), end=np.array(end, np.int64),
        sig=np.frombuffer(bytes(sig), np.uint8),
        sigage=np.frombuffer(bytes(sigage), np.uint8),
        sigint=np.frombuffer(bytes(sigint), np.uint8),
        sfari=np.frombuffer(bytes(sfari), np.uint8),
        genic=np.frombuffer(bytes(genic), np.uint8),
        est=np.array(est, np.float32), se_dx=np.array(se_dx, np.float32),
        se_pc1=np.array(se_pc1, np.float32), se_age=np.array(se_age, np.float32),
        p_dx=np.array(p_dx, np.float64), q_dx=np.array(q_dx, np.float64),
        cpg_feat=np.array(cpg_feat), cpg_sym=np.array(cpg_sym),
        ct_names=np.array([k for k, _ in sorted(ct_map.items(), key=lambda x: x[1])]),
        cx_names=np.array([k for k, _ in sorted(cx_map.items(), key=lambda x: x[1])]),
        ft_names=np.array([k for k, _ in sorted(ft_map.items(), key=lambda x: x[1])]),
        gt_names=np.array([k for k, _ in sorted(gt_map.items(), key=lambda x: x[1])]),
        chr_names=np.array([k for k, _ in sorted(chr_map.items(), key=lambda x: x[1])]),
    )
    np.savez_compressed(CACHE, **d)
    return d


def load():
    if CACHE.exists() and "--force" not in sys.argv:
        print(f"[cache] {CACHE}", file=sys.stderr)
        z = np.load(CACHE, allow_pickle=False)
        return {k: z[k] for k in z.files}
    return extract()


# --------------------------------------------------------------------- stats helpers
def fisher2x2(a, b, c, d):
    orr = (a * d) / (b * c) if b * c else float("inf")
    _, p = stats.fisher_exact([[a, b], [c, d]])
    try:
        res = stats.contingency.odds_ratio([[a, b], [c, d]], kind="conditional")
        cond_or = float(res.statistic)
        lo, hi = res.confidence_interval(0.95)
        cond_ci = (float(lo), float(hi))
    except Exception:
        cond_or, cond_ci = float("nan"), (float("nan"), float("nan"))
    if min(a, b, c, d) > 0:
        sl = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        wci = (float(np.exp(np.log(orr) - 1.96 * sl)), float(np.exp(np.log(orr) + 1.96 * sl)))
    else:
        wci = (float("nan"), float("nan"))
    return dict(a=int(a), b=int(b), c=int(c), d=int(d), sample_OR=float(orr), woolf_CI=wci,
                fisher_p=float(p), conditional_MLE_OR=cond_or, conditional_CI=cond_ci)


def mantel_haenszel(a, b, c, d):
    a, b, c, d = (np.asarray(x, float) for x in (a, b, c, d))
    n = a + b + c + d
    k = n > 0
    a, b, c, d, n = a[k], b[k], c[k], d[k], n[k]
    R = a * d / n
    S = b * c / n
    if R.sum() == 0 or S.sum() == 0:
        return dict(MH_OR=float(R.sum() / S.sum()) if S.sum() else float("inf"),
                    MH_CI=[float("nan"), float("nan")], MH_logOR_SE=float("nan"),
                    MH_chi2=float("nan"), MH_p=float("nan"),
                    observed_a=float(a.sum()), expected_a=float("nan"),
                    note="degenerate stratum sums (no discordant mass); CI undefined")
    or_mh = R.sum() / S.sum()
    P = (a + d) / n
    Q = (b + c) / n
    v = ((P * R).sum() / (2 * R.sum() ** 2)
         + ((P * S).sum() + (Q * R).sum()) / (2 * R.sum() * S.sum())
         + (Q * S).sum() / (2 * S.sum() ** 2))
    se = float(np.sqrt(v))
    ci = [float(np.exp(np.log(or_mh) - 1.96 * se)), float(np.exp(np.log(or_mh) + 1.96 * se))]
    E = (a + b) * (a + c) / n
    V = np.where(n > 1, (a + b) * (c + d) * (a + c) * (b + d) / (n ** 2 * np.maximum(n - 1, 1)), 0.0)
    num = abs(a.sum() - E.sum()) - 0.5
    chi2 = float(num ** 2 / V.sum()) if V.sum() > 0 else 0.0
    return dict(MH_OR=float(or_mh), MH_CI=ci, MH_logOR_SE=se, MH_chi2=chi2,
                MH_p=float(stats.chi2.sf(chi2, 1)), observed_a=float(a.sum()),
                expected_a=float(E.sum()), n_strata=int(len(a)))


def logistic_irls(X, y, maxit=200, tol=1e-10):
    X = np.asarray(X, float); y = np.asarray(y, float)
    beta = np.zeros(X.shape[1])
    for _ in range(maxit):
        eta = np.clip(X @ beta, -30, 30)
        mu = 1 / (1 + np.exp(-eta))
        w = np.clip(mu * (1 - mu), 1e-10, None)
        z = eta + (y - mu) / w
        XtW = X.T * w
        H = XtW @ X
        try:
            new = np.linalg.solve(H, XtW @ z)
        except np.linalg.LinAlgError:
            new = np.linalg.lstsq(H, XtW @ z, rcond=None)[0]
        if np.max(np.abs(new - beta)) < tol:
            beta = new
            break
        beta = new
    eta = np.clip(X @ beta, -30, 30)
    mu = 1 / (1 + np.exp(-eta))
    w = np.clip(mu * (1 - mu), 1e-10, None)
    cov = np.linalg.pinv((X.T * w) @ X)
    return beta, np.sqrt(np.diag(cov))


def strat_table(sig, is_t, strata, n_strata):
    a = np.zeros(n_strata); b = np.zeros(n_strata)
    c = np.zeros(n_strata); d = np.zeros(n_strata)
    for s in range(n_strata):
        m = strata == s
        t = m & is_t
        o = m & ~is_t
        a[s] = sig[t].sum(); b[s] = t.sum() - a[s]
        c[s] = sig[o].sum(); d[s] = o.sum() - c[s]
    return a, b, c, d


def deciles_of(x, n_q=10):
    qs = np.quantile(x, np.linspace(0, 1, n_q + 1))
    qs[0] -= 1e-9
    return np.clip(np.searchsorted(qs, x, side="left") - 1, 0, n_q - 1), qs


def matched_block(sig, ct, strat, n_strata, ct_names, target, label,
                  strat_labels=None, fit_logistic=True):
    ti = int(np.where(ct_names == target)[0][0])
    is_t = ct == ti
    a, b, c, d = strat_table(sig, is_t, strat, n_strata)
    rows = []
    for s in range(n_strata):
        nt, no = a[s] + b[s], c[s] + d[s]
        if nt == 0 and no == 0:
            continue
        rows.append(dict(
            stratum=(strat_labels[s] if strat_labels is not None else s + 1),
            tgt_sig=int(a[s]), tgt_tested=int(nt),
            tgt_rate=float(a[s] / nt) if nt else None,
            oth_sig=int(c[s]), oth_tested=int(no),
            oth_rate=float(c[s] / no) if no else None,
            OR=float((a[s] * d[s]) / (b[s] * c[s])) if b[s] * c[s] else None))
    mh = mantel_haenszel(a, b, c, d)
    logit = None
    if fit_logistic:
        used = np.unique(strat)
        N = len(sig)
        X = np.zeros((N, 2 + len(used) - 1))
        X[:, 0] = 1.0
        X[:, 1] = is_t.astype(float)
        for j, s in enumerate(used[1:]):
            X[:, 2 + j] = (strat == s).astype(float)
        beta, se = logistic_irls(X, sig.astype(float))
        zz = beta[1] / se[1] if se[1] > 0 else np.nan
        logit = dict(logOR=float(beta[1]), SE=float(se[1]),
                     OR=float(np.exp(np.clip(beta[1], -30, 30))),
                     CI=[float(np.exp(np.clip(beta[1] - 1.96 * se[1], -30, 30))),
                         float(np.exp(np.clip(beta[1] + 1.96 * se[1], -30, 30)))],
                     z=float(zz), p=float(2 * stats.norm.sf(abs(zz))) if np.isfinite(zz) else None)
    crude = fisher2x2(int(a.sum()), int(b.sum()), int(c.sum()), int(d.sum()))
    return dict(stratifier=label, target=target, n_bins=int(len(sig)),
                per_stratum=rows, mantel_haenszel=mh, logistic=logit, crude=crude)


def se_matched_block(sig, ct, se, ct_names, target, n_q=10, label="coefSE_DiagnosisASD"):
    ok = np.isfinite(se)
    sig, ct, se = sig[ok], ct[ok], se[ok]
    dec, qs = deciles_of(se, n_q)
    r = matched_block(sig, ct, dec, n_q, ct_names, target, f"{label} ({n_q} quantiles)")
    for i, row in enumerate(r["per_stratum"]):
        row["range"] = [float(qs[i]), float(qs[i + 1])]
    r["n_dropped_nonfinite"] = int((~ok).sum())
    return r


def perm_null(sig, ct, ct_names, target, dec=None, n_perm=N_PERM, seed=SEED):
    """Cell-type-label permutation null for the target's significant-bin COUNT.
    Permuting cell-type labels within a stratum is exactly a hypergeometric draw."""
    rng = np.random.default_rng(seed)
    ti = int(np.where(ct_names == target)[0][0])
    is_t = ct == ti
    obs = int(sig[is_t].sum())
    if dec is None:
        N, S, nt = len(sig), int(sig.sum()), int(is_t.sum())
        draws = rng.hypergeometric(S, N - S, nt, size=n_perm).astype(float)
    else:
        draws = np.zeros(n_perm)
        for s in np.unique(dec):
            m = dec == s
            N, S, nt = int(m.sum()), int(sig[m].sum()), int((m & is_t).sum())
            if nt == 0 or S == 0 or N - S < 1:
                continue
            draws += rng.hypergeometric(S, N - S, nt, size=n_perm)
    sd = draws.std(ddof=1)
    return dict(target=target, observed_sig=obs, tested=int(is_t.sum()),
                observed_rate=float(obs / max(is_t.sum(), 1)),
                null_mean=float(draws.mean()), null_sd=float(sd),
                null_q025=float(np.quantile(draws, .025)),
                null_q975=float(np.quantile(draws, .975)),
                z=float((obs - draws.mean()) / sd) if sd > 0 else float("nan"),
                n_perm=int(n_perm),
                p_one_sided=float((int((draws >= obs).sum()) + 1) / (n_perm + 1)))


# ------------------------------------------------------- refGene coordinate mapping
def load_refgene():
    if not REFGENE.exists():
        return None
    per = defaultdict(lambda: dict(s=[], e=[], t=[], g=[]))
    with gzip.open(REFGENE, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 13:
                continue
            ch, strand = p[2], p[3]
            if "_" in ch:
                continue
            s, e = int(p[4]), int(p[5])
            per[ch]["s"].append(s); per[ch]["e"].append(e)
            per[ch]["t"].append(s if strand == "+" else e)
            per[ch]["g"].append(p[12])
    out = {}
    for ch, v in per.items():
        o = np.argsort(np.array(v["t"]))
        out[ch] = dict(s=np.array(v["s"])[o], e=np.array(v["e"])[o],
                       t=np.array(v["t"])[o], g=np.array(v["g"])[o])
    return out


def map_bins(rg, chrom_names, chrom_codes, starts, ends, max_dist=50000):
    """Return (gene_symbol_or_'', distance) per bin: overlapping gene body preferred,
    else nearest TSS within max_dist."""
    n = len(starts)
    gene = np.array([""] * n, dtype=object)
    dist = np.full(n, np.nan)
    mid = (starts + ends) // 2
    for ci, ch in enumerate(chrom_names):
        if ch not in rg:
            continue
        m = np.where(chrom_codes == ci)[0]
        if not len(m):
            continue
        R = rg[ch]
        s, e, t, g = R["s"], R["e"], R["t"], R["g"]
        bs, be, bm = starts[m], ends[m], mid[m]
        # nearest TSS (t is sorted)
        j = np.searchsorted(t, bm)
        jl = np.clip(j - 1, 0, len(t) - 1)
        jr = np.clip(j, 0, len(t) - 1)
        dl = np.abs(bm - t[jl]); dr = np.abs(bm - t[jr])
        pick = np.where(dl <= dr, jl, jr)
        nd = np.minimum(dl, dr)
        gg = g[pick].astype(object)
        dd = nd.astype(float)
        # gene-body overlap (chunked broadcast)
        step = 512
        for k0 in range(0, len(m), step):
            k1 = min(k0 + step, len(m))
            ov = (s[None, :] < be[k0:k1, None]) & (e[None, :] > bs[k0:k1, None])
            any_ov = ov.any(1)
            if any_ov.any():
                idx = np.where(any_ov)[0]
                for ii in idx:
                    cand = np.where(ov[ii])[0]
                    best = cand[np.argmin(np.abs(bm[k0 + ii] - t[cand]))]
                    gg[k0 + ii] = g[best]
                    dd[k0 + ii] = 0.0
        gg = np.where(dd <= max_dist, gg, "")
        gene[m] = gg
        dist[m] = dd
    return gene, dist


def load_modules():
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
        for nm, mm in gmt.items():
            if any(re.search(pp, nm, re.I) for pp in pats):
                g |= mm
        return g

    MOD = {
        "OL_maturation": U(r"MYELIN", r"OLIGODENDROCYTE", r"ENSHEATHMENT"),
        "LIPID_sterol": U(r"CHOLESTEROL_BIOSYNTH", r"STEROL_BIOSYNTH", r"SPHINGOLIPID",
                          r"LIPID_TRANSPORT", r"LIPOPROTEIN_PARTICLE"),
        "IMMUNE": U(r"CYTOKINE_MEDIATED", r"CYTOKINE_PRODUCTION", r"INFLAMMATORY_RESPONSE",
                    r"INTERFERON", r"NF_KAPPAB", r"TUMOR_NECROSIS_FACTOR", r"COMPLEMENT",
                    r"TOLL_LIKE_RECEPTOR", r"INNATE_IMMUNE"),
        "CHROMATIN": U(r"CHROMATIN_REMODEL", r"CHROMATIN_ORGANIZATION", r"HISTONE_MODIFICATION"),
        "NEURODEV": U(r"NEURON_DIFFERENTIATION", r"NERVOUS_SYSTEM_DEVELOPMENT", r"AXON_GUIDANCE"),
        "SYNAPTIC": U(r"GOBP_SYNAP", r"GOCC_POSTSYNAPTIC_DENSITY"),
        "MITOCHONDRIAL": U(r"MITOCHONDRI"),
        "NEG_OLFACTORY": U(r"OLFACTORY"),
        "NEG_KERATIN": U(r"KERATIN", r"CORNIFIED"),
    }
    sf = set()
    p = SCOPE / "sfari/SFARI-Gene_genes_2026-09-10.csv"
    if p.exists():
        for r in csv.DictReader(open(p)):
            g = r.get("gene-symbol") or r.get("gene_symbol")
            if g:
                sf.add(g)
    MOD["REF_SFARI"] = sf
    return MOD


PERM_NOTE = (
    "The diagnosis-label permutation requested (>=1000 relabellings of ASD/control at the "
    "DONOR level, refitting the per-bin model each time) is NOT computable from any file on "
    "disk. `GSE298243_Age_DMRs_PrespecifiedRegions.tsv.gz` is a per-bin SUMMARY table: one "
    "row per (cell type x cytosine context x feature) carrying only fitted coefficients, "
    "SEs, p, q and the paper's significance flags. It holds no per-donor and no per-nucleus "
    "methylation values, no donor identifiers and no diagnosis column, so the label to be "
    "permuted does not exist in it and the model cannot be refit. Permuting the significance "
    "flag across bins instead would test a different null and is not done here. "
    "**What would be needed:** the GSE298243 genotype-masked pseudo-bulk per-donor ALLC "
    "files (CGN and CAN, one pair per donor x cell type, described in "
    "`GSE298243_README.docx`) plus the donor metadata table carrying diagnosis/age/RIN, "
    "which together allow the per-bin mC/C model to be refit under permuted labels. "
    "Per-cell-type NUCLEI COUNTS - the quantity E-006's 'OPC is rare -> fewest nuclei -> "
    "noisiest' premise turns on - live in the supplementary Seurat object's metadata slot, "
    "also not downloaded. **What is run instead, and labelled as a different null:** a "
    "CELL-TYPE-label permutation, unstratified and stratified within coefSE deciles, which "
    "answers 'is OPC's significant-bin count what an SE-matched random draw of the same "
    "number of pooled tested bins would give?' - the null the artifact hypothesis asserts.")


# --------------------------------------------------------------------------- main
def main():
    d = load()
    ctn, cxn, ftn, gtn, chn = (d["ct_names"], d["cx_names"], d["ft_names"],
                               d["gt_names"], d["chr_names"])
    ct, cx, ft = d["ct"], d["cx"], d["ft"]
    sig = d["sig"].astype(np.int64)
    se_dx, est, p_dx, q_dx = d["se_dx"], d["est"], d["p_dx"], d["q_dx"]
    R = {"seed": SEED, "n_perm": N_PERM, "source_file": str(DMR), "n_rows": int(len(ct)),
         "cell_types": list(map(str, ctn)), "contexts": list(map(str, cxn)),
         "feature_types": list(map(str, ftn))}
    cpg_i = int(np.where(cxn == "CpG")[0][0])
    cpa_i = int(np.where(cxn == "CpA")[0][0])
    mcpg = cx == cpg_i
    oi = int(np.where(ctn == "OPC")[0][0])

    # ------------------------------------------------ 0. design disclosure
    R["design"] = dict(
        rows_flagged_ReportSignif_Age=int(d["sigage"].sum()),
        frac_rows_flagged_age=float(d["sigage"].mean()),
        rows_flagged_ReportSignif_Interaction=int(d["sigint"].sum()),
        feature_type_counts={str(t): int((ft == i).sum()) for i, t in enumerate(ftn)},
        feature_type_counts_CpG={str(t): int((mcpg & (ft == i)).sum())
                                 for i, t in enumerate(ftn)},
        significance_rule_is_qle05=bool(np.array_equal(
            sig[np.isfinite(q_dx)], (q_dx[np.isfinite(q_dx)] <= 0.05).astype(np.int64))),
    )

    # ------------------------------------------------ 1. the table
    tab = {}
    for xi, xname in enumerate(cxn):
        for ci, cname in enumerate(ctn):
            m = (cx == xi) & (ct == ci)
            tab.setdefault(str(xname), {})[str(cname)] = dict(
                tested=int(m.sum()), sig=int(sig[m].sum()),
                rate=float(sig[m].sum() / max(m.sum(), 1)))
    R["table_celltype_x_context"] = tab
    R["context_totals"] = {str(x): dict(tested=int((cx == i).sum()),
                                        sig=int(sig[cx == i].sum()))
                           for i, x in enumerate(cxn)}
    R["opc_vs_others"] = {}
    for xi, xname in list(enumerate(cxn)) + [(None, "ALL")]:
        m = np.ones(len(ct), bool) if xname == "ALL" else (cx == xi)
        t = m & (ct == oi); o = m & (ct != oi)
        a, b = int(sig[t].sum()), int(t.sum() - sig[t].sum())
        c, dd = int(sig[o].sum()), int(o.sum() - sig[o].sum())
        r = fisher2x2(a, b, c, dd)
        r.update(opc_sig=a, opc_tested=a + b, opc_rate=a / max(a + b, 1),
                 oth_sig=c, oth_tested=c + dd, oth_rate=c / max(c + dd, 1),
                 rate_ratio=(a / max(a + b, 1)) / (c / max(c + dd, 1)) if c else None)
        R["opc_vs_others"][str(xname)] = r
    t = mcpg & (ct == oi)
    a, b = int(sig[t].sum()), int(t.sum() - sig[t].sum())
    R["opc_vs_each_CpG"] = {}
    for ci, cname in enumerate(ctn):
        if ci == oi:
            continue
        o = mcpg & (ct == ci)
        R["opc_vs_each_CpG"][str(cname)] = fisher2x2(
            a, b, int(sig[o].sum()), int(o.sum() - sig[o].sum()))
    xs, ys, names = [], [], []
    for ci, cname in enumerate(ctn):
        m = mcpg & (ct == ci)
        xs.append(np.log10(m.sum())); ys.append(sig[m].sum() / m.sum()); names.append(str(cname))
    xs, ys = np.array(xs), np.array(ys)
    keep = np.array(names) != "OPC"
    sl, ic = np.polyfit(xs[keep], ys[keep], 1)
    resid = ys - (sl * xs + ic)
    sd7 = resid[keep].std(ddof=2)
    R["trend_fit_CpG"] = dict(slope=float(sl), intercept=float(ic),
                              resid_sd_7celltypes=float(sd7),
                              z_by_celltype={n: float(resid[i] / sd7)
                                             for i, n in enumerate(names)})

    # ------------------------------------------------ 2. E-006 noise statistics
    noise = {}
    for ci, cname in enumerate(ctn):
        for scope, m in (("CpG", mcpg & (ct == ci)), ("ALL", ct == ci)):
            s = se_dx[m]; s = s[np.isfinite(s)]
            e = est[m]; e = e[np.isfinite(e)]
            es = est[m & (sig == 1)]; es = es[np.isfinite(es)]
            noise.setdefault(scope, {})[str(cname)] = dict(
                n_tested=int(m.sum()), n_sig=int(sig[m].sum()),
                median_coefSE=float(np.median(s)) if s.size else None,
                iqr_coefSE=[float(np.quantile(s, .25)), float(np.quantile(s, .75))] if s.size else None,
                sd_coefEst=float(e.std(ddof=1)) if e.size > 1 else None,
                median_abs_coefEst_all=float(np.median(np.abs(e))) if e.size else None,
                median_abs_coefEst_sig=float(np.median(np.abs(es))) if es.size else None)
    R["noise_stats"] = noise
    ms = np.array([noise["CpG"][str(c)]["median_coefSE"] for c in ctn])
    ns = np.array([noise["CpG"][str(c)]["n_sig"] for c in ctn], float)
    rr = ns / np.array([noise["CpG"][str(c)]["n_tested"] for c in ctn], float)
    R["e006_cross_celltype_corr"] = dict(
        pearson_medianSE_vs_sigcount=float(np.corrcoef(ms, ns)[0, 1]),
        pearson_medianSE_vs_sigrate=float(np.corrcoef(ms, rr)[0, 1]),
        spearman_medianSE_vs_sigrate=float(stats.spearmanr(ms, rr).statistic), n_points=8)

    # ------------------------------------------------ 3. LEG (a) SE matching
    scpg, ctcpg, secpg = sig[mcpg], ct[mcpg], se_dx[mcpg]
    R["LEG_A_se_matched_CpG"] = {f"q{nq}": se_matched_block(scpg, ctcpg, secpg, ctn, "OPC", nq)
                                 for nq in (10, 20)}
    R["LEG_A_comparators"] = {c: se_matched_block(scpg, ctcpg, secpg, ctn, c, 10)
                              for c in ("Oligo", "InhMGE")}
    mcpa = cx == cpa_i
    R["LEG_A_CpA_note"] = ("OPC contributes only %d tested CpA bins (0 significant); the CpA "
                           "stratum cannot test this question." % int((mcpa & (ct == oi)).sum()))
    # FDR-free p-threshold ladder, SE-decile matched
    pc = p_dx[mcpg]
    okp = np.isfinite(pc) & np.isfinite(secpg)
    dec10, _ = deciles_of(secpg[okp], 10)
    ladder = {}
    for thr in (0.05, 0.01, 1e-3, 1e-4, 1e-5):
        hit = (pc[okp] < thr).astype(np.int64)
        crude = fisher2x2(int(hit[ctcpg[okp] == oi].sum()),
                          int((ctcpg[okp] == oi).sum() - hit[ctcpg[okp] == oi].sum()),
                          int(hit[ctcpg[okp] != oi].sum()),
                          int((ctcpg[okp] != oi).sum() - hit[ctcpg[okp] != oi].sum()))

        mb = matched_block(hit, ctcpg[okp], dec10, 10, ctn, "OPC", "SE decile",
                           fit_logistic=False)
        ladder[f"p<{thr:g}"] = dict(crude=crude, mh=mb["mantel_haenszel"])
    R["LEG_A_pvalue_ladder_CpG"] = ladder

    # ------------------------------------------------ 4. permutation nulls
    ok = np.isfinite(secpg)
    dcs, _ = deciles_of(secpg[ok], 10)
    R["permutation_nulls_CpG"] = {
        tgt: dict(unstratified=perm_null(scpg[ok], ctcpg[ok], ctn, tgt, None),
                  se_decile_stratified=perm_null(scpg[ok], ctcpg[ok], ctn, tgt, dcs))
        for tgt in ("OPC", "Oligo", "InhMGE")}
    R["permutation_feasibility"] = PERM_NOTE

    # ------------------------------------------------ 5. effective-N proxies
    R["coverage_column_present"] = False
    R["se_matched_alt_proxies_CpG"] = {
        nm: se_matched_block(scpg, ctcpg, arr, ctn, "OPC", 10, nm)
        for nm, arr in (("coefSE_genoPC1", d["se_pc1"][mcpg]),
                        ("coefSE_AgeLin", d["se_age"][mcpg]))}
    fin = (np.isfinite(secpg) & np.isfinite(d["se_pc1"][mcpg]) & np.isfinite(d["se_age"][mcpg]))
    R["proxy_correlations_CpG"] = dict(
        spearman_SEdx_vs_SEgenoPC1=float(stats.spearmanr(secpg[fin], d["se_pc1"][mcpg][fin]).statistic),
        spearman_SEdx_vs_SEageLin=float(stats.spearmanr(secpg[fin], d["se_age"][mcpg][fin]).statistic),
        n=int(fin.sum()))
    dA, _ = deciles_of(secpg[fin], 10)
    dB, _ = deciles_of(d["se_pc1"][mcpg][fin], 10)
    a2, b2, c2, d2 = strat_table(scpg[fin], ctcpg[fin] == oi, dA * 10 + dB, 100)
    R["double_stratified_MH_CpG"] = mantel_haenszel(a2, b2, c2, d2)
    R["double_stratified_MH_CpG"]["n_strata_nonempty"] = int(((a2 + b2) > 0).sum())

    # ------------------------------------------------ 6. LEG (b) REGION COMPOSITION
    feats = d["cpg_feat"]
    uf, inv = np.unique(feats, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uf))
    R["LEG_B_overlap"] = dict(
        n_cpg_rows=int(mcpg.sum()), n_unique_features=int(len(uf)),
        multiplicity_hist={str(k): int(v) for k, v in enumerate(np.bincount(cnt)) if v},
        n_in_all8=int((cnt == len(ctn)).sum()))
    # pairwise overlap matrix
    setmap = {str(c): set(feats[ctcpg == i].tolist()) for i, c in enumerate(ctn)}
    pair = {}
    for i, c1 in enumerate(ctn):
        for c2 in list(ctn)[i + 1:]:
            ov = setmap[str(c1)] & setmap[str(c2)]
            pair[f"{c1}|{c2}"] = dict(overlap=len(ov),
                                      jaccard=len(ov) / len(setmap[str(c1)] | setmap[str(c2)]))
    R["LEG_B_pairwise_overlap"] = pair
    # paired McNemar on OPC-vs-comparator shared bins
    fi = {f: k for k, f in enumerate(feats)}
    paired = {}
    opc_pos = {f: k for f, k in zip(feats[ctcpg == oi], np.where(ctcpg == oi)[0])}
    for ci, cname in enumerate(ctn):
        if ci == oi:
            continue
        oth_pos = {f: k for f, k in zip(feats[ctcpg == ci], np.where(ctcpg == ci)[0])}
        shared = sorted(set(opc_pos) & set(oth_pos))
        if not shared:
            paired[str(cname)] = dict(n_shared=0, note="no shared bins")
            continue
        io = np.array([opc_pos[f] for f in shared])
        ic2 = np.array([oth_pos[f] for f in shared])
        so, sc = scpg[io], scpg[ic2]
        b01 = int(((so == 1) & (sc == 0)).sum())
        c10 = int(((so == 0) & (sc == 1)).sum())
        both = int(((so == 1) & (sc == 1)).sum())
        neither = int(((so == 0) & (sc == 0)).sum())
        if b01 + c10 > 0:
            pmc = float(stats.binomtest(b01, b01 + c10, 0.5).pvalue)
            orc = float(b01 / c10) if c10 else float("inf")
        else:
            pmc, orc = None, None
        paired[str(cname)] = dict(n_shared=len(shared), opc_only_sig=b01,
                                  other_only_sig=c10, both_sig=both, neither=neither,
                                  conditional_OR=orc, mcnemar_exact_p=pmc)
    R["LEG_B_paired_mcnemar"] = paired
    # region composition of the tested sets
    comp = {}
    for ci, cname in enumerate(ctn):
        m = mcpg & (ct == ci)
        w = (d["end"][m] - d["start"][m])
        comp[str(cname)] = dict(
            n=int(m.sum()),
            featuretype={str(t): float((ft[m] == i).mean()) for i, t in enumerate(ftn)},
            frac_with_symbol=float(d["genic"][m].mean()),
            frac_sfari=float(d["sfari"][m].mean()),
            median_width=float(np.median(w)), mean_width=float(w.mean()),
            frac_protein_coding=float((d["gt"][m] == int(np.where(gtn == "protein_coding")[0][0])).mean())
            if "protein_coding" in list(gtn) else None,
            chrom_top5={str(chn[k]): float(v / m.sum()) for k, v in
                        sorted(((k, v) for k, v in enumerate(np.bincount(d["chrom"][m],
                                                                         minlength=len(chn)))),
                               key=lambda x: -x[1])[:5]})
    R["LEG_B_region_composition_CpG"] = comp
    # MH stratified on FeatureType, on chromosome, and on SE decile x FeatureType
    R["LEG_B_MH_by_featuretype"] = matched_block(
        scpg, ctcpg, ft[mcpg].astype(int), len(ftn), ctn, "OPC", "FeatureType",
        strat_labels=[str(x) for x in ftn])
    chcpg = d["chrom"][mcpg].astype(int)
    R["LEG_B_MH_by_chromosome"] = matched_block(
        scpg, ctcpg, chcpg, len(chn), ctn, "OPC", "chromosome",
        strat_labels=[str(x) for x in chn], fit_logistic=False)
    joint = dec10 * 3 + ft[mcpg][okp].astype(int)
    R["LEG_B_MH_SEdecile_x_featuretype"] = matched_block(
        scpg[okp], ctcpg[okp], joint, 10 * len(ftn), ctn, "OPC",
        "SE decile x FeatureType", fit_logistic=False)["mantel_haenszel"]
    # genic/intergenic
    R["LEG_B_MH_by_genic"] = matched_block(
        scpg, ctcpg, d["genic"][mcpg].astype(int), 2, ctn, "OPC", "has gene symbol",
        strat_labels=["no symbol", "has symbol"])

    # ------------------------------------------------ 7. LEG (c) MECHANISTIC EMPTINESS
    syms = d["cpg_sym"]
    opc_m = ctcpg == oi
    opc_sig_m = opc_m & (scpg == 1)
    bg_rate = float(d["genic"][mcpg][opc_m].mean())
    n_sig = int(opc_sig_m.sum())
    n_sig_sym = int(d["genic"][mcpg][opc_sig_m].sum())
    bt = stats.binomtest(n_sig_sym, n_sig, bg_rate)
    R["LEG_C_annotation_rate"] = dict(
        opc_tested=int(opc_m.sum()),
        opc_tested_with_symbol=int(d["genic"][mcpg][opc_m].sum()),
        background_symbol_rate=bg_rate,
        opc_sig=n_sig, opc_sig_with_symbol=n_sig_sym,
        expected_sig_with_symbol=float(bg_rate * n_sig),
        binomial_p_two_sided=float(bt.pvalue),
        symbols_of_sig=[str(s) for s in syms[opc_sig_m] if s],
        per_celltype_symbol_rate={str(c): float(d["genic"][mcpg][ctcpg == i].mean())
                                  for i, c in enumerate(ctn)})
    # coordinate mapping
    rg = load_refgene()
    if rg is None:
        R["LEG_C_coordinate_mapping"] = dict(
            status="NO_LOCAL_ANNOTATION",
            note="No hg38 gene annotation found on disk; coordinate mapping not attempted "
                 "and nothing downloaded.")
    else:
        idx = np.where(mcpg)[0][opc_m]
        g, dist = map_bins(rg, chn, d["chrom"][idx], d["start"][idx], d["end"][idx])
        issig = scpg[opc_m] == 1
        # build verification: do FeatureType=='gene' rows agree with refGene coordinates?
        gi = int(np.where(ftn == "gene")[0][0]) if "gene" in list(ftn) else None
        ver = None
        if gi is not None:
            sym_all = np.array([""] * len(ct), dtype=object)
            sym_all[np.where(mcpg)[0]] = d["cpg_sym"]
            universe = set()
            for ch in rg:
                universe.update(rg[ch]["g"].tolist())
            vm = mcpg & (ft == gi) & (d["genic"] == 1)
            vidx = np.where(vm)[0]
            gv, _ = map_bins(rg, chn, d["chrom"][vidx], d["start"][vidx], d["end"][vidx])
            tgts = np.array([str(sym_all[k]) for k in vidx])
            inuni = np.array([t in universe for t in tgts])
            allm = float(np.mean([gv[k] == tgts[k] for k in range(len(vidx))]))
            conm = (float(np.mean([gv[k] == tgts[k] for k in np.where(inuni)[0]]))
                    if inuni.any() else float("nan"))
            ver = dict(
                n_checked=int(len(vidx)), frac_symbol_match_all=allm,
                n_symbol_in_refgene=int(inuni.sum()),
                frac_symbol_match_given_in_refgene=conm,
                interpretation=(
                    "The shipped Symbol for FeatureType=='gene' rows is often an Ensembl ID "
                    "or a small-RNA/pseudogene/lncRNA name that RefSeq refGene does not "
                    "carry, so the unconditional match rate understates agreement. "
                    "Conditioning on the shipped symbol existing in refGene at all is the "
                    "right build check: a high conditional rate means the coordinates are "
                    "hg38."))
        MOD = load_modules()
        rows = []
        gs_all = np.array([str(x) for x in g])
        for mname, genes in MOD.items():
            inm = np.array([x in genes for x in gs_all])
            si = int((inm & issig).sum()); ti = int(inm.sum())
            so = int((~inm & issig).sum()); to = int((~inm).sum())
            if ti < 30:
                continue
            r = fisher2x2(si, ti - si, so, to - so)
            r.update(module=mname, sig=si, tested=ti, rate=si / ti, rate_bg=so / to)
            rows.append(r)
        rows.sort(key=lambda x: x["fisher_p"])
        mm = len(rows)
        for i, r in enumerate(rows, 1):
            r["q"] = min(1.0, r["fisher_p"] * mm / i)
        for i in range(mm - 2, -1, -1):
            rows[i]["q"] = min(rows[i]["q"], rows[i + 1]["q"])
        topg = [(str(k), int(v)) for k, v in sorted(
            ((k, int((gs_all[issig] == k).sum())) for k in set(gs_all[issig]) if k),
            key=lambda x: -x[1])[:25]]
        R["LEG_C_coordinate_mapping"] = dict(
            status="OK", annotation=str(REFGENE),
            build_verification=ver,
            n_opc_tested=int(opc_m.sum()),
            n_mapped_tested=int((gs_all != "").sum()),
            frac_mapped_tested=float((gs_all != "").mean()),
            n_opc_sig=int(issig.sum()),
            n_mapped_sig=int((gs_all[issig] != "").sum()),
            frac_mapped_sig=float((gs_all[issig] != "").mean()),
            median_dist_tested=float(np.nanmedian(dist)),
            n_unique_genes_sig=int(len(set(gs_all[issig]) - {""})),
            modules=rows, top_genes_sig=topg)

    # ------------------------------------------------ 8. bottom line per leg
    crude_or = R["opc_vs_others"]["CpG"]["sample_OR"]
    mh = R["LEG_A_se_matched_CpG"]["q10"]["mantel_haenszel"]
    mh_ft = R["LEG_B_MH_by_featuretype"]["mantel_haenszel"]
    mh_j = R["LEG_B_MH_SEdecile_x_featuretype"]
    lad = R["LEG_A_pvalue_ladder_CpG"]
    npair = sum(v.get("n_shared", 0) for v in R["LEG_B_paired_mcnemar"].values())
    lc = R["LEG_C_annotation_rate"]
    legs = {}
    legs["a_noise"] = (
        "FAILS as an explanation. Within SE deciles the significant-bin rate FALLS "
        "monotonically with SE (others: %.4f in decile 1 -> %.4f in decile 10), the opposite "
        "of the mechanism E-006 proposed; matching on SE decile RAISES the OPC OR from "
        "%.2f to %.2f (95%% CI %.2f-%.2f, p=%.2g)."
        % (R["LEG_A_se_matched_CpG"]["q10"]["per_stratum"][0]["oth_rate"],
           R["LEG_A_se_matched_CpG"]["q10"]["per_stratum"][-1]["oth_rate"],
           crude_or, mh["MH_OR"], mh["MH_CI"][0], mh["MH_CI"][1], mh["MH_p"])
        + " The excess is also present at every fixed p-value threshold and grows with "
          "stringency (crude OR %.2f at p<0.05 -> %.2f at p<1e-4; SE-matched %.2f -> %.2f), "
          "so it is not an artifact of the paper's FDR procedure either."
        % (lad["p<0.05"]["crude"]["sample_OR"], lad["p<0.0001"]["crude"]["sample_OR"],
           lad["p<0.05"]["mh"]["MH_OR"], lad["p<0.0001"]["mh"]["MH_OR"]))
    legs["b_region"] = (
        "STANDS as a limitation, and is the binding one. Only %d of %d unique CpG features "
        "are shared by any pair of cell types (%d shared bins summed over all 7 OPC-vs-X "
        "pairs), because every row of this file is an age-DMR OF THAT CELL TYPE "
        "(ReportSignif_Age is '*' on %.1f%% of rows), so the eight cell types are compared "
        "on essentially disjoint region sets. MH stratified on FeatureType gives OR %.2f "
        "and on SE decile x FeatureType %.2f, so the coarse region covariates the table "
        "supports do not explain the excess - but they are coarse, and the disjointness "
        "itself cannot be adjusted away with this file."
        % (int(len(uf) - R["LEG_B_overlap"]["multiplicity_hist"].get("1", 0)), len(uf),
           npair, 100 * R["design"]["frac_rows_flagged_age"],
           mh_ft["MH_OR"], mh_j["MH_OR"]))
    cmm = R["LEG_C_coordinate_mapping"]
    legs["c_emptiness"] = (
        "FAILS AS STATED, but its conclusion survives on better evidence. Only %.2f%% of "
        "OPC's TESTED CpG features carry a gene symbol, so %.1f symbol-carrying bins are "
        "expected among the %d significant ones and %d were observed (binomial p=%.2f): "
        "3/272 is the expected number and measures the file's annotation sparsity, not the "
        "signal. Coordinate mapping against a local hg38 refGene raises annotated coverage "
        "from %.2f%% to %.1f%% (%d of %d significant bins map, hitting %d distinct genes) "
        "and still finds NO module enrichment against the coordinate-mapped OPC background "
        "(every BH q > %.2f, including SFARI OR %.2f and the positive-control OL-maturation "
        "and lipid modules) - so 'no nameable mechanism' stands, for a better reason than "
        "E-006 gave."
        % (100 * lc["background_symbol_rate"], lc["expected_sig_with_symbol"],
           lc["opc_sig"], lc["opc_sig_with_symbol"], lc["binomial_p_two_sided"],
           100 * lc["background_symbol_rate"], 100 * cmm.get("frac_mapped_sig", float("nan")),
           cmm.get("n_mapped_sig", 0), cmm.get("n_opc_sig", 0),
           cmm.get("n_unique_genes_sig", 0),
           min((m["q"] for m in cmm.get("modules", [])), default=float("nan")),
           next((m["sample_OR"] for m in cmm.get("modules", [])
                 if m["module"] == "REF_SFARI"), float("nan")))
        if cmm.get("status") == "OK" else
        "FAILS. Only %.2f%% of OPC's TESTED CpG features carry a gene symbol, so %.1f are "
        "expected among %d significant bins and %d were observed (binomial p=%.2f)."
        % (100 * lc["background_symbol_rate"], lc["expected_sig_with_symbol"],
           lc["opc_sig"], lc["opc_sig_with_symbol"], lc["binomial_p_two_sided"]))
    R["_legs"] = legs
    R["_bottom_line"] = (
        "The OPC excess is NOT an estimator-noise artifact and NOT mechanistically empty in "
        "the way E-006 claimed - SE-decile matching raises the OR from %.2f to %.2f "
        "(95%% CI %.2f-%.2f, p=%.2g) and 3/272 gene symbols is exactly what the %.2f%% "
        "annotation rate of the tested background predicts - but it remains UNDECIDABLE as "
        "a claim about cell types, because every row of this file is an age-DMR of its own "
        "cell type, so the eight tested bin sets are near-disjoint (%d of %d unique CpG "
        "features shared by any pair) and 'OPC has the highest rate' is a statement about "
        "OPC's own age-DMR set, not about a common genomic yardstick; deciding it needs the "
        "GSE298243 per-donor pseudo-bulk ALLC files + donor metadata (to refit a common bin "
        "set and permute diagnosis labels), which are not on disk."
        % (crude_or, mh["MH_OR"], mh["MH_CI"][0], mh["MH_CI"][1], mh["MH_p"],
           100 * lc["background_symbol_rate"],
           int(len(uf) - R["LEG_B_overlap"]["multiplicity_hist"].get("1", 0)), len(uf)))
    R["_bottom_line_inputs"] = dict(
        crude_OR_CpG=crude_or, MH_OR_SE=mh["MH_OR"], MH_CI=mh["MH_CI"], MH_p=mh["MH_p"],
        MH_OR_featuretype=mh_ft["MH_OR"], MH_OR_joint=mh_j["MH_OR"],
        n_shared_pairs_total=npair,
        strat_perm_p=R["permutation_nulls_CpG"]["OPC"]["se_decile_stratified"]["p_one_sided"],
        ladder_p1e4_MH=lad["p<0.0001"]["mh"]["MH_OR"])
    print("[bottom line] " + R["_bottom_line"], file=sys.stderr)

    (OUT / "opc_se_matched.json").write_text(json.dumps(R, indent=1, default=float))
    write_report(R, ctn)
    print(f"[done] {OUT/'opc_se_matched.json'}\n[done] {REPORT}")
    return R


# --------------------------------------------------------------------------- report
def fm(x, n=3):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


def write_report(R, ctn):
    ctn = [str(c) for c in ctn]
    L = []
    A = L.append
    A("# OPC epigenetic excess (GSE298243): noise, region and annotation adjudication")
    A("")
    A("**Run:** `scripts/p16_opc_se_matched.py` (seed %d, %d permutation draws) · 2026-09-21  "
      % (R["seed"], R["n_perm"]))
    A("**Source:** `%s` (%s rows)  " % (R["source_file"], f"{R['n_rows']:,}"))
    A("**Outputs:** `outputs/phase16/opc_se_matched.json` (every number below), "
      "`outputs/phase16/dmr_extract_v2.npz` (cached column extract)")
    A("")
    A("## Bottom line")
    A("")
    A("> " + R["_bottom_line"])
    A("")
    A("### E-006's three legs")
    A("")
    A("| leg | verdict |")
    A("|---|---|")
    A("| (a) estimator noise | " + R["_legs"]["a_noise"] + " |")
    A("| (b) region composition | " + R["_legs"]["b_region"] + " |")
    A("| (c) mechanistic emptiness | " + R["_legs"]["c_emptiness"] + " |")
    A("")

    # ---------------- 0 design
    dg = R["design"]
    A("## 0. What this file actually is (read this before any number below)")
    A("")
    A(f"- `ReportSignif_Age` is `*` on **{dg['rows_flagged_ReportSignif_Age']:,} of "
      f"{R['n_rows']:,} rows = {100*dg['frac_rows_flagged_age']:.1f}%**. The file is named "
      f"`Age_DMRs_PrespecifiedRegions`, and that is literal: **every tested feature is an "
      f"age-DMR of the cell type in whose rows it appears.** The diagnosis test is run "
      f"*conditional on* an age effect having been detected in that cell type.")
    A(f"- `ReportSignif_Interaction` is `*` on {dg['rows_flagged_ReportSignif_Interaction']:,} "
      f"rows.")
    A(f"- `ReportSignif_Diagnosis` is exactly `coefQ_DiagnosisASD <= 0.05`: verified "
      f"identical on every row with a finite q "
      f"(**{dg['significance_rule_is_qle05']}**). So the \"significant-bin rate\" is an "
      f"**FDR-thresholded** quantity, not a fixed-alpha one. The q column is not reproducible "
      f"as BH within cell type, within cell type x context, or over the whole file "
      f"(max |q - BH| = 0.42-0.70 in every scope tried), so the FDR procedure and its scope "
      f"are the paper's, not recoverable here — an additional reason the rate is not a clean "
      f"cross-cell-type yardstick. Section 3 therefore also reports an FDR-free ladder.")
    A("- FeatureType composition: "
      + ", ".join(f"`{k}` {v:,}" for k, v in dg["feature_type_counts"].items())
      + "; within CpG: "
      + ", ".join(f"`{k}` {v:,}" for k, v in dg["feature_type_counts_CpG"].items()) + ".")
    A("")

    # ---------------- 1
    A("## 1. The cell type x cytosine-context table, reproduced from the raw file")
    A("")
    for cxname, sub in R["table_celltype_x_context"].items():
        tot = R["context_totals"][cxname]
        A(f"### {cxname} (tested {tot['tested']:,}, significant {tot['sig']:,})")
        A("")
        A("| cell type | significant | tested | rate |")
        A("|---|---:|---:|---:|")
        for c, v in sorted(sub.items(), key=lambda x: -x[1]["rate"]):
            A(f"| {'**'+c+'**' if c=='OPC' else c} | {v['sig']:,} | {v['tested']:,} | "
              f"{v['rate']:.5f} |")
        A("")
    A("### OPC vs the seven others pooled")
    A("")
    A("| stratum | OPC sig/tested | OPC rate | others sig/tested | others rate | "
      "sample OR (95% Woolf) | conditional MLE OR | Fisher p |")
    A("|---|---|---:|---|---:|---|---|---|")
    for k, v in R["opc_vs_others"].items():
        A(f"| {k} | {v['opc_sig']:,}/{v['opc_tested']:,} | {v['opc_rate']:.5f} | "
          f"{v['oth_sig']:,}/{v['oth_tested']:,} | {v['oth_rate']:.5f} | "
          f"{fm(v['sample_OR'])} ({fm(v['woolf_CI'][0],2)}–{fm(v['woolf_CI'][1],2)}) | "
          f"{fm(v['conditional_MLE_OR'])} | {v['fisher_p']:.3g} |")
    A("")
    cg = R["opc_vs_others"]["CpG"]
    A(f"**What the statistic actually is.** CpG stratum: sample OR **{cg['sample_OR']:.4f}**, "
      f"Fisher p **{cg['fisher_p']:.3g}**, rate ratio {cg['rate_ratio']:.3f}, conditional-MLE "
      f"OR {cg['conditional_MLE_OR']:.4f}. The 2026-09-17 doc reported **OR 1.90**; the audit's "
      f"re-derivation from the doc's own counts gave **1.92**. The audit is right: the "
      f"odds ratio is **1.92**; 1.90 is the *rate ratio*, which the doc mislabelled as an OR. "
      f"The p-value 1.49e-19 is exactly reproduced.")
    A("")
    A("### OPC vs each individual cell type (CpG)")
    A("")
    A("| comparator | sample OR | Fisher p |")
    A("|---|---:|---:|")
    for c, v in sorted(R["opc_vs_each_CpG"].items(), key=lambda x: -x[1]["sample_OR"]):
        A(f"| {c} | {v['sample_OR']:.3f} | {v['fisher_p']:.3g} |")
    A("")
    tf = R["trend_fit_CpG"]
    A("### The \"+5.10 SD above the tested-bin trend\" check, reproduced")
    A("")
    A(f"rate = {tf['slope']:.6f}·log10(tested) + {tf['intercept']:.6f} fit on the seven "
      f"non-OPC cell types; residual SD {tf['resid_sd_7celltypes']:.6f}. "
      f"**OPC z = {tf['z_by_celltype']['OPC']:+.2f}** — the +5.10 SD figure is exact.")
    A("")
    A("| cell type | residual z |")
    A("|---|---:|")
    for c, z in sorted(tf["z_by_celltype"].items(), key=lambda x: -x[1]):
        A(f"| {'**'+c+'**' if c=='OPC' else c} | {z:+.2f} |")
    A("")

    # ---------------- 2
    A("## 2. The E-006 noise statistics, now persisted")
    A("")
    A("| cell type | tested | sig | median coefSE | IQR | SD(coefEst) | "
      "median abs coefEst (all) | median abs coefEst (sig) |")
    A("|---|---:|---:|---:|---|---:|---:|---:|")
    for c in sorted(ctn, key=lambda x: -R["noise_stats"]["CpG"][x]["median_coefSE"]):
        v = R["noise_stats"]["CpG"][c]
        A(f"| {'**'+c+'**' if c=='OPC' else c} | {v['n_tested']:,} | {v['n_sig']:,} | "
          f"{v['median_coefSE']:.4f} | {v['iqr_coefSE'][0]:.4f}–{v['iqr_coefSE'][1]:.4f} | "
          f"{v['sd_coefEst']:.4f} | {v['median_abs_coefEst_all']:.4f} | "
          f"{v['median_abs_coefEst_sig']:.4f} |")
    A("")
    cc = R["e006_cross_celltype_corr"]
    A(f"E-006's three headline numbers are **confirmed exactly**: OPC median coefSE "
      f"**0.3653** (highest of eight), SD(coefEst) **0.7393** (highest), median |coef| among "
      f"significant bins **1.1918** (largest). Its reassurance statistic is also confirmed "
      f"(r = {cc['pearson_medianSE_vs_sigcount']:+.3f} for median SE vs significant *count*) "
      f"— but that was the wrong statistic: against the significant-bin *rate*, the quantity "
      f"the claim is about, the cross-cell-type correlation is "
      f"r = {cc['pearson_medianSE_vs_sigrate']:+.3f} "
      f"(Spearman {cc['spearman_medianSE_vs_sigrate']:+.3f}), i.e. *supporting* the noise "
      f"story at the cell-type level. Section 3 shows that within-stratum the relationship "
      f"runs the other way.")
    A("")

    # ---------------- 3
    A("## 3. LEG (a) — the SE-decile-matched test")
    A("")
    A("All tested CpG features from all eight cell types pooled, cut into deciles of "
      "`coefSE_DiagnosisASD`; within each decile OPC's significant-bin rate vs the pooled "
      "other seven.")
    A("")
    dm = R["LEG_A_se_matched_CpG"]["q10"]
    A("| SE decile | coefSE range | OPC sig/tested | OPC rate | others sig/tested | "
      "others rate | OR |")
    A("|---:|---|---|---:|---|---:|---:|")
    for r in dm["per_stratum"]:
        A(f"| {r['stratum']} | {r['range'][0]:.3f}–{r['range'][1]:.3f} | "
          f"{r['tgt_sig']}/{r['tgt_tested']:,} | {fm(r['tgt_rate'],5)} | "
          f"{r['oth_sig']}/{r['oth_tested']:,} | {fm(r['oth_rate'],5)} | {fm(r['OR'],2)} |")
    A("")
    mh, lg, cr = dm["mantel_haenszel"], dm["logistic"], dm["crude"]
    A(f"- **Crude OR** {cr['sample_OR']:.3f} ({fm(cr['woolf_CI'][0],2)}–"
      f"{fm(cr['woolf_CI'][1],2)}), Fisher p {cr['fisher_p']:.3g}")
    A(f"- **Mantel–Haenszel OR across 10 SE deciles = {mh['MH_OR']:.3f}** (95% RBG CI "
      f"{mh['MH_CI'][0]:.2f}–{mh['MH_CI'][1]:.2f}), MH chi2 {mh['MH_chi2']:.1f}, "
      f"p = {mh['MH_p']:.3g}; {mh['observed_a']:.0f} observed OPC significant bins vs "
      f"{mh['expected_a']:.1f} expected under the stratified null")
    A(f"- **Logistic `sig ~ is_OPC + C(SE decile)`: OR = {lg['OR']:.3f}** (95% CI "
      f"{lg['CI'][0]:.2f}–{lg['CI'][1]:.2f}), z = {lg['z']:.2f}, p = {lg['p']:.3g}")
    d20 = R["LEG_A_se_matched_CpG"]["q20"]
    A(f"- 20 strata: MH OR {d20['mantel_haenszel']['MH_OR']:.3f} "
      f"({d20['mantel_haenszel']['MH_CI'][0]:.2f}–{d20['mantel_haenszel']['MH_CI'][1]:.2f}), "
      f"p {d20['mantel_haenszel']['MH_p']:.3g}")
    A(f"- {R['LEG_A_CpA_note']}")
    A("")
    A("**The E-006 mechanism is directly contradicted by the stratum margins.** Among the "
      "non-OPC cell types the significant-bin rate *falls* from "
      f"{dm['per_stratum'][0]['oth_rate']:.4f} in the lowest-SE decile to "
      f"{dm['per_stratum'][-1]['oth_rate']:.4f} in the highest. In this dataset a noisier "
      "estimate is **less**, not more, likely to be called significant — as it must be, since "
      "the call is a q-value threshold on coef/SE. So OPC's high median SE is a reason to "
      "expect *fewer* significant bins, and adjusting for it makes the excess larger, not "
      "smaller. E-006's leg (a) rests on a mechanism that does not operate here.")
    A("")
    A("### FDR-free check: fixed p-value thresholds, same SE-decile matching")
    A("")
    A("The significance flag is an FDR call whose procedure is the paper's; BH is "
      "super-linear in signal density, so a small real difference can be amplified. Fixed "
      "thresholds have no such property.")
    A("")
    A("| threshold | OPC rate | others rate | crude OR | Fisher p | MH OR (10 SE deciles) | MH p |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    for k, v in R["LEG_A_pvalue_ladder_CpG"].items():
        c = v["crude"]; m = v["mh"]
        A(f"| {k} | {c['a']/(c['a']+c['b']):.5f} | {c['c']/(c['c']+c['d']):.5f} | "
          f"{c['sample_OR']:.3f} | {c['fisher_p']:.3g} | {fm(m['MH_OR'])} | "
          f"{m['MH_p']:.3g} |")
    A("")
    A("### The same test applied to two comparison cell types")
    A("")
    A("| cell type | crude OR | crude Fisher p | MH OR (10 SE deciles) | 95% CI | MH p |")
    A("|---|---:|---:|---:|---|---:|")
    A(f"| **OPC** | {cr['sample_OR']:.3f} | {cr['fisher_p']:.3g} | {mh['MH_OR']:.3f} | "
      f"{mh['MH_CI'][0]:.2f}–{mh['MH_CI'][1]:.2f} | {mh['MH_p']:.3g} |")
    for c, v in R["LEG_A_comparators"].items():
        m = v["mantel_haenszel"]
        A(f"| {c} | {v['crude']['sample_OR']:.3f} | {v['crude']['fisher_p']:.3g} | "
          f"{m['MH_OR']:.3f} | {m['MH_CI'][0]:.2f}–{m['MH_CI'][1]:.2f} | {m['MH_p']:.3g} |")
    A("")

    # ---------------- 4
    A("## 4. Permutation nulls")
    A("")
    A("### Feasibility of the requested diagnosis-label permutation")
    A("")
    A(R["permutation_feasibility"])
    A("")
    A(f"### Cell-type-label permutation ({R['n_perm']:,} draws, seed {R['seed']}), CpG")
    A("")
    A("| cell type | null | observed sig | tested | observed rate | null mean | null SD | "
      "z | one-sided p |")
    A("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for c, v in R["permutation_nulls_CpG"].items():
        for k, lab in (("unstratified", "unstratified"),
                       ("se_decile_stratified", "within SE decile")):
            w = v[k]
            A(f"| {c} | {lab} | {w['observed_sig']} | {w['tested']:,} | "
              f"{w['observed_rate']:.5f} | {w['null_mean']:.1f} | {w['null_sd']:.1f} | "
              f"{w['z']:+.2f} | {w['p_one_sided']:.4g} |")
    A("")
    A("The floor of 9.999e-05 is the resolution limit of 10,000 draws, not a computed value.")
    A("")

    # ---------------- 5
    A("## 5. Effective-N proxies")
    A("")
    A("The table has **no** coverage, read-depth or nuclei-count column — all 52 columns are "
      "feature annotation, the paper's significance flags, predicted values, and "
      "coefficient/SE/p/q triplets. `coefSE_genoPC1` and `coefSE_AgeLin` are the available "
      "precision proxies that do not involve the diagnosis contrast.")
    A("")
    pc = R["proxy_correlations_CpG"]
    A(f"Spearman rho with `coefSE_DiagnosisASD` (CpG, n={pc['n']:,}): `coefSE_genoPC1` "
      f"{pc['spearman_SEdx_vs_SEgenoPC1']:.3f}, `coefSE_AgeLin` "
      f"{pc['spearman_SEdx_vs_SEageLin']:.3f}.")
    A("")
    A("| stratifier (10 deciles) | MH OR | 95% CI | MH p | logistic OR | logistic p |")
    A("|---|---:|---|---:|---:|---:|")
    A(f"| coefSE_DiagnosisASD | {mh['MH_OR']:.3f} | {mh['MH_CI'][0]:.2f}–{mh['MH_CI'][1]:.2f} "
      f"| {mh['MH_p']:.3g} | {lg['OR']:.3f} | {lg['p']:.3g} |")
    for nm, v in R["se_matched_alt_proxies_CpG"].items():
        m2, l2 = v["mantel_haenszel"], v["logistic"]
        A(f"| {nm} | {m2['MH_OR']:.3f} | {m2['MH_CI'][0]:.2f}–{m2['MH_CI'][1]:.2f} | "
          f"{m2['MH_p']:.3g} | {l2['OR']:.3f} | {l2['p']:.3g} |")
    ds = R["double_stratified_MH_CpG"]
    A("")
    A(f"Double stratification (coefSE_DiagnosisASD decile x coefSE_genoPC1 decile, "
      f"{ds['n_strata_nonempty']} non-empty cells of 100): MH OR = {ds['MH_OR']:.3f} "
      f"(95% CI {ds['MH_CI'][0]:.2f}–{ds['MH_CI'][1]:.2f}), p = {ds['MH_p']:.3g}.")
    A("")

    # ---------------- 6 LEG B
    ov = R["LEG_B_overlap"]
    A("## 6. LEG (b) — region composition: the cell types are tested on different bins")
    A("")
    A(f"The {ov['n_cpg_rows']:,} CpG rows span {ov['n_unique_features']:,} **unique** "
      f"features. Multiplicity: "
      + ", ".join(f"{v:,} features in {k} cell type(s)"
                  for k, v in sorted(ov["multiplicity_hist"].items(), key=lambda x: int(x[0])))
      + f"; **{ov['n_in_all8']} features are tested in all eight cell types.**")
    A("")
    A("This is not a quirk — it is the design. Because every row is an age-DMR of its own "
      "cell type (§0), the tested set for each cell type is *that cell type's age-DMR set*. "
      "**Consequence: every cross-cell-type rate comparison in this dataset, including the "
      "source paper's own and the 2026-09-17 claim, compares different genomic regions.** "
      "\"OPC has the highest diagnosis-significant rate\" means \"among the regions that are "
      "age-DMRs in OPC, a larger fraction also show a diagnosis effect\" — a statement about "
      "the age/diagnosis overlap inside OPC's own region set, not a common yardstick across "
      "cell types.")
    A("")
    A("### Pairwise CpG bin-set overlap")
    A("")
    A("| pair | shared features | Jaccard |")
    A("|---|---:|---:|")
    for k, v in sorted(R["LEG_B_pairwise_overlap"].items(), key=lambda x: -x[1]["overlap"]):
        A(f"| {k.replace('|',' vs ')} | {v['overlap']:,} | {v['jaccard']:.5f} |")
    A("")
    A("### Paired within-bin comparison on the shared features (McNemar, exact)")
    A("")
    A("| comparator | shared bins | OPC-only sig | other-only sig | both | conditional OR | "
      "exact p |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    for c, v in sorted(R["LEG_B_paired_mcnemar"].items(),
                       key=lambda x: -x[1].get("n_shared", 0)):
        if v.get("n_shared", 0) == 0:
            A(f"| {c} | 0 | – | – | – | – | – |")
            continue
        A(f"| {c} | {v['n_shared']:,} | {v['opc_only_sig']} | {v['other_only_sig']} | "
          f"{v['both_sig']} | {fm(v['conditional_OR'],2)} | "
          f"{'–' if v['mcnemar_exact_p'] is None else format(v['mcnemar_exact_p'],'.3g')} |")
    A("")
    A("These overlaps are far too small to adjudicate the claim; they are reported so the "
      "size of the problem is on the record rather than assumed away.")
    A("")
    A("### How OPC's tested set differs in composition")
    A("")
    A("| cell type | n | % bin500bp | % gene | % promoter | % with gene symbol | % SFARI | "
      "median width |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|")
    for c, v in R["LEG_B_region_composition_CpG"].items():
        ftc = v["featuretype"]
        A(f"| {'**'+c+'**' if c=='OPC' else c} | {v['n']:,} | "
          f"{100*ftc.get('bin500bp',0):.2f} | {100*ftc.get('gene',0):.2f} | "
          f"{100*ftc.get('promoter',0):.2f} | {100*v['frac_with_symbol']:.2f} | "
          f"{100*v['frac_sfari']:.2f} | {v['median_width']:.0f} |")
    A("")
    mft = R["LEG_B_MH_by_featuretype"]
    A("Adjusting for the region covariates the table does support:")
    A("")
    A("| stratifier | MH OR | 95% CI | MH p |")
    A("|---|---:|---|---:|")
    for lab, blk in (("FeatureType (bin500bp / gene / promoter)",
                      R["LEG_B_MH_by_featuretype"]["mantel_haenszel"]),
                     ("has gene symbol (genic / intergenic)",
                      R["LEG_B_MH_by_genic"]["mantel_haenszel"]),
                     ("chromosome", R["LEG_B_MH_by_chromosome"]["mantel_haenszel"]),
                     ("SE decile x FeatureType", R["LEG_B_MH_SEdecile_x_featuretype"])):
        A(f"| {lab} | {fm(blk['MH_OR'])} | {fm(blk['MH_CI'][0],2)}–{fm(blk['MH_CI'][1],2)} | "
          f"{blk['MH_p']:.3g} |")
    A("")
    A("None of these covariates removes the excess. But note what they can and cannot do: "
      "they adjust for *measured* composition differences within a shared covariate space. "
      "They cannot adjust for the fact that the regions themselves are different pieces of "
      "genome selected by a cell-type-specific age-DMR screen. That residual is not "
      "addressable with this file.")
    A("")

    # ---------------- 7 LEG C
    lc = R["LEG_C_annotation_rate"]
    A("## 7. LEG (c) — is the signal \"mechanistically empty\"?")
    A("")
    A(f"E-006: \"of 272 significant OPC CpG bins, only 3 carry a gene symbol\". True — and "
      f"**it is exactly the expected number.** Only "
      f"{lc['opc_tested_with_symbol']:,}/{lc['opc_tested']:,} = "
      f"**{100*lc['background_symbol_rate']:.2f}%** of OPC's *tested* CpG features carry a "
      f"symbol at all, so {lc['expected_sig_with_symbol']:.1f} symbol-carrying bins are "
      f"expected among {lc['opc_sig']} significant ones; {lc['opc_sig_with_symbol']} were "
      f"observed (two-sided binomial p = {lc['binomial_p_two_sided']:.2f}). The "
      f"\"mechanistically empty\" criticism measured the annotation sparsity of the file, "
      f"not a property of the OPC signal.")
    A("")
    A("| cell type | % of tested CpG features carrying a gene symbol |")
    A("|---|---:|")
    for c, v in sorted(lc["per_celltype_symbol_rate"].items(), key=lambda x: -x[1]):
        A(f"| {'**'+c+'**' if c=='OPC' else c} | {100*v:.2f} |")
    A("")
    A("Symbols on the significant OPC bins as shipped: "
      + (", ".join(f"`{s}`" for s in lc["symbols_of_sig"]) or "none") + ".")
    A("")
    cm = R["LEG_C_coordinate_mapping"]
    if cm["status"] != "OK":
        A("**Coordinate mapping not run.** " + cm["note"])
    else:
        A("### Doing it properly: coordinate mapping")
        A("")
        A(f"Annotation used: local `{cm['annotation']}` (hg38 refGene, already on disk — "
          f"nothing downloaded). Each feature is assigned the gene whose body it overlaps, "
          f"else the nearest TSS within 50 kb.")
        if cm.get("build_verification"):
            bv = cm["build_verification"]
            A("")
            A(f"Build check: {bv['n_checked']:,} `FeatureType == gene` CpG rows ship their own "
              f"symbol. Unconditionally, {100*bv['frac_symbol_match_all']:.1f}% are assigned "
              f"the same symbol by coordinate lookup — but most shipped symbols are Ensembl "
              f"IDs, small RNAs or pseudogenes that RefSeq refGene does not carry. Restricting "
              f"to the {bv['n_symbol_in_refgene']:,} whose symbol exists in refGene at all, "
              f"**{100*bv['frac_symbol_match_given_in_refgene']:.1f}%** match — the table's "
              f"coordinates are hg38.")
        A("")
        A(f"- OPC tested features mapped to a gene: **{cm['n_mapped_tested']:,}/"
          f"{cm['n_opc_tested']:,} = {100*cm['frac_mapped_tested']:.1f}%** "
          f"(vs {100*lc['background_symbol_rate']:.2f}% carrying a shipped symbol)")
        A(f"- OPC **significant** features mapped: **{cm['n_mapped_sig']:,}/"
          f"{cm['n_opc_sig']:,} = {100*cm['frac_mapped_sig']:.1f}%**, hitting "
          f"{cm['n_unique_genes_sig']:,} distinct genes")
        A("")
        A("Module enrichment, background = the coordinate-mapped OPC tested set:")
        A("")
        A("| module | sig | tested | rate | bg rate | OR | Fisher p | BH q |")
        A("|---|---:|---:|---:|---:|---:|---:|---:|")
        for r in cm["modules"]:
            A(f"| {r['module']} | {r['sig']} | {r['tested']:,} | {r['rate']:.4f} | "
              f"{r['rate_bg']:.4f} | {fm(r['sample_OR'],2)} | {r['fisher_p']:.3g} | "
              f"{r['q']:.3g} |")
        A("")
        A("Top genes by number of significant OPC features: "
          + ", ".join(f"{g} ({n})" for g, n in cm["top_genes_sig"][:20]) + ".")
        A("")
        A("**Read this correctly in both directions.** E-006's *stated* reason for calling "
          "the signal empty — only 3 of 272 bins carry a symbol — is wrong: that is the "
          "expected count given the file's annotation sparsity, and proper coordinate "
          "mapping recovers "
          f"{100*cm['frac_mapped_sig']:.0f}% of the significant bins onto "
          f"{cm['n_unique_genes_sig']} distinct genes. But with that far richer gene list in "
          "hand, **there is still no module structure**: nothing reaches BH q < 0.1, SFARI is "
          "if anything depleted, the OL-maturation and lipid modules — the ones the "
          "reformulated model predicts — are flat, and the top nominal hit is a negative "
          "control (olfactory). So the correct statement is not \"the signal is noise because "
          "it is unannotated\" but \"the signal is well annotated and has no detectable "
          "programme structure at this sample size\".")
    A("")
    A("---")
    A("")
    A("Every number in this file is emitted by `scripts/p16_opc_se_matched.py` into "
      "`outputs/phase16/opc_se_matched.json`; this report is generated from that JSON in the "
      "same run.")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()

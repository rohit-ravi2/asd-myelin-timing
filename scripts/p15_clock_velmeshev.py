#!/usr/bin/env python3
"""
p15 STEP 2 -- apply the atlas-2 (GSE155488) OL-lineage maturation CLOCK to Velmeshev 2019
ASD cortex snRNA-seq, and test for a DONOR-LEVEL stage shift in ASD vs control.

GENERATIVE / EXPLORATORY. The deliverable is a map + an honest power statement.
No result is called positive unless the control battery passes.

------------------------------------------------------------------------------------
THE CLOCK (definition, reused verbatim from scripts/p6_atlas2_replication.py)
------------------------------------------------------------------------------------
For a set of cells C and five stage marker panels OPC/COP/NFOL/MFOL/MOL:
  1. X = log1p( count / cell_total_UMI * 1e4 )                       (log-CPM10k)
  2. v_s = mean over the panel's genes present in the matrix          (s = 0..4)
  3. z_s = (v_s - mean_C v_s) / sd_C v_s                              z-scored ACROSS CELLS
  4. w   = softmax over s of z                                        (p6: exp(z-max)/sum)
  5. stage(cell) = sum_s w_s * s                                      soft stage in [0,4]
p6 z-scores across the pooled atlas; p10 z-scores per sample. We use the p6 pooling
(all nuclei of one lineage pooled) because per-donor z-scoring would by construction
destroy the between-donor contrast this analysis is about.

Internal validation of the clock in atlas 2 (p6, outputs/phase6): mean assigned stage
rises with donor age without being told the ages -- fetal 0.51 -> 2yr 1.91 -> 13yr 2.35
-> adult ~2.2. NOTE the dynamic range is almost entirely fetal->2yr; the postnatal
slope 2yr->adult is ~0.3-0.45 stage units and NON-monotonic at the adult end. Velmeshev
donors are 4-22 yr, i.e. entirely inside the flat part of the clock's demonstrated
range. This is a sensitivity limitation and is reported as such.

------------------------------------------------------------------------------------
SEED 20260921. Re-runnable: python3 scripts/p15_extract_velmeshev.py && python3 this.
------------------------------------------------------------------------------------
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy import stats

BASE = Path(__file__).resolve().parents[1]
VEL = Path("data/dvad-corpus/asd_scope/velmeshev2019")
EXT = BASE / "outputs/phase15/extract"
OUT = BASE / "outputs/phase15"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20260921
N_PERM = 10000
N_NULL_SETS = 500
N_SUBSAMPLE = 200
rng = np.random.default_rng(SEED)

sys.path.insert(0, str(BASE / "scripts"))
from p15_extract_velmeshev import STAGE_P6, STAGE_P10, ASTRO_PANEL, NEURON_PANEL  # noqa: E402

ORDER5 = ["OPC", "COP", "NFOL", "MFOL", "MOL"]
LINEAGES = {
    "OL":  ["OPC", "Oligodendrocytes"],
    "AST": ["AST-PP", "AST-FB"],
    "EXC": ["L2/3", "L4", "L5/6", "L5/6-CC"],
}


# ====================================================================== stats helpers
def ols(y, cols, names):
    """OLS with t-distribution p-values. cols: list of 1-D covariate arrays."""
    X = np.column_stack([np.ones(len(y))] + list(cols))
    ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    X, yy = X[ok], y[ok]
    b, *_ = np.linalg.lstsq(X, yy, rcond=None)
    r = yy - X @ b
    df = len(yy) - X.shape[1]
    s2 = (r @ r) / df
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = b / se
    p = 2 * stats.t.sf(np.abs(t), df)
    return dict(names=["const"] + names, beta=b.tolist(), se=se.tolist(),
                t=t.tolist(), p=p.tolist(), df=int(df), n=int(len(yy)))


def ols_t_first(y, cols):
    """Fast path: t-stat of the FIRST covariate only (for permutation loops)."""
    X = np.column_stack([np.ones(len(y))] + list(cols))
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    df = len(y) - X.shape[1]
    if df <= 0:
        return np.nan
    s2 = (r @ r) / df
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(cov[1, 1])
    return b[1] / se if se > 0 else np.nan


def cohen_d(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    return float((a.mean() - b.mean()) / sp) if sp > 0 else np.nan


def mde_two_sample(n1, n2, alpha=0.05, power=0.80, extra_df_loss=0):
    """Minimum detectable Cohen's d by exact noncentral-t search."""
    df = n1 + n2 - 2 - extra_df_loss
    crit = stats.t.ppf(1 - alpha / 2, df)
    lo, hi = 0.01, 5.0
    for _ in range(200):
        d = (lo + hi) / 2
        nc = d * np.sqrt(1 / (1 / n1 + 1 / n2))
        pw = stats.nct.sf(crit, df, nc) + stats.nct.cdf(-crit, df, nc)
        if pw < power:
            lo = d
        else:
            hi = d
    return float((lo + hi) / 2)


# ====================================================================== load
def load():
    t0 = time.time()
    z = np.load(EXT / "velmeshev_subset.npz", allow_pickle=True)
    sym = z["symbols"].astype(str)
    gene_total = z["gene_total"]
    cell_total = z["cell_total"].astype(np.float64)
    cell_nnz = z["cell_nnz"].astype(np.float64)
    eg, ec, ev = z["ent_gene"], z["ent_cell"], z["ent_val"]

    meta_rows = []
    with open(VEL / "meta.tsv") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            meta_rows.append(line.rstrip("\n").split("\t"))
    H = {k: i for i, k in enumerate(hdr)}
    M = np.array(meta_rows, dtype=object)
    meta = dict(
        cell=M[:, H["cell"]].astype(str), cluster=M[:, H["cluster"]].astype(str),
        sample=M[:, H["sample"]].astype(str), individual=M[:, H["individual"]].astype(str),
        region=M[:, H["region"]].astype(str), age=M[:, H["age"]].astype(float),
        sex=M[:, H["sex"]].astype(str), diagnosis=M[:, H["diagnosis"]].astype(str),
        capbatch=M[:, H["Capbatch"]].astype(str), seqbatch=M[:, H["Seqbatch"]].astype(str),
        pmi=M[:, H["post-mortem interval (hours)"]].astype(float),
        rin=M[:, H["RNA Integrity Number"]].astype(float),
        genes=M[:, H["genes"]].astype(float), umis=M[:, H["UMIs"]].astype(float),
        mito=M[:, H["RNA mitochondr. percent"]].astype(float),
        ribo=M[:, H["RNA ribosomal percent"]].astype(float),
    )
    # barcode order check
    bc = np.array([l.strip() for l in open(VEL / "raw/barcodes.tsv")])
    assert len(bc) == len(meta["cell"]), (len(bc), len(meta["cell"]))
    n_match = int((bc == meta["cell"]).sum())
    print(f"[align] barcodes vs meta cell ids identical for {n_match}/{len(bc)}",
          file=sys.stderr)
    assert n_match == len(bc), "barcode/meta order mismatch -- refusing to continue"
    # UMI total check
    rho = stats.pearsonr(cell_total, meta["umis"])[0]
    print(f"[check] mtx per-cell total vs meta UMIs: r={rho:.6f}, "
          f"max|diff|={np.abs(cell_total-meta['umis']).max():.0f}", file=sys.stderr)
    print(f"[load] {time.time()-t0:.0f}s", file=sys.stderr)
    return sym, gene_total, cell_total, cell_nnz, eg, ec, ev, meta


def build_logcpm(sym, gene_total, cell_total, eg, ec, ev, wanted):
    """Dense genes x cells log-CPM10k for `wanted` symbols. Duplicated symbols -> the
    matrix row with the highest total count."""
    ncell = len(cell_total)
    rows, present = {}, []
    for g in wanted:
        idx = np.flatnonzero(sym == g)
        if idx.size == 0:
            continue
        r = int(idx[np.argmax(gene_total[idx])]) + 1   # 1-based mtx row
        rows[g] = r
        present.append(g)
    row2k = {r: k for k, r in enumerate(rows[g] for g in present)}
    X = np.zeros((len(present), ncell), dtype=np.float32)
    sel = np.isin(eg, np.array(list(row2k.keys()), dtype=eg.dtype))
    gg, cc, vv = eg[sel], ec[sel], ev[sel].astype(np.float32)
    kk = np.array([row2k[int(r)] for r in gg], dtype=np.int64)
    X[kk, cc - 1] = vv
    tot = cell_total.copy()
    tot[tot <= 0] = 1
    X = np.log1p(X / tot[None, :] * 1e4).astype(np.float32)
    return X, present


def soft_stage(X, ci, panels, order):
    """Steps 2-5 of the clock on the columns of X (already restricted to one lineage)."""
    sc = []
    for s in order:
        ix = [ci[m] for m in panels[s] if m in ci]
        v = X[ix].mean(axis=0) if ix else np.zeros(X.shape[1], dtype=np.float32)
        sc.append((v - v.mean()) / (v.std() + 1e-9))
    S = np.vstack(sc)
    w = np.exp(S - S.max(axis=0))
    w /= w.sum(axis=0)
    stage = (w * np.arange(len(order))[:, None]).sum(axis=0)
    return stage.astype(np.float64), np.argmax(S, axis=0)


# ====================================================================== donor tables
def donor_table(meta, mask):
    """One row per individual over the cells selected by `mask`."""
    ind = meta["individual"][mask]
    uniq = sorted(set(ind.tolist()))
    idx = {d: np.flatnonzero(ind == d) for d in uniq}
    D = dict(
        donor=np.array(uniq),
        n=np.array([len(idx[d]) for d in uniq], float),
        asd=np.array([1.0 if meta["diagnosis"][mask][idx[d]][0] == "ASD" else 0.0
                      for d in uniq]),
        age=np.array([meta["age"][mask][idx[d]][0] for d in uniq]),
        female=np.array([1.0 if meta["sex"][mask][idx[d]][0] == "F" else 0.0 for d in uniq]),
        pmi=np.array([meta["pmi"][mask][idx[d]][0] for d in uniq]),
        rin=np.array([meta["rin"][mask][idx[d]][0] for d in uniq]),
        mlogumi=np.array([np.log10(meta["umis"][mask][idx[d]]).mean() for d in uniq]),
        mloggene=np.array([np.log10(meta["genes"][mask][idx[d]]).mean() for d in uniq]),
        mmito=np.array([meta["mito"][mask][idx[d]].mean() for d in uniq]),
        fracPFC=np.array([(meta["region"][mask][idx[d]] == "PFC").mean() for d in uniq]),
    )
    return D, idx


def donor_means(vals, idx, donors):
    return np.array([vals[idx[d]].mean() for d in donors])


def full_test(y, D, label, extra=None, extra_names=None):
    """Welch + Mann-Whitney + OLS(dx, age, sex) + OLS(+depth) + donor-level permutation."""
    a, c = y[D["asd"] == 1], y[D["asd"] == 0]
    cols = [D["asd"], D["age"], D["female"]]
    names = ["ASD", "age", "female"]
    if extra is not None:
        cols = cols + list(extra)
        names = names + list(extra_names)
    m = ols(y, cols, names)
    mdepth = ols(y, [D["asd"], D["age"], D["female"], D["mlogumi"], D["mloggene"],
                     np.log10(D["n"])],
                 ["ASD", "age", "female", "mean_log10UMI", "mean_log10genes", "log10_ncells"])
    # donor-level permutation of the diagnosis label, refitting the adjusted model
    t_obs = m["t"][1]
    base = [D["age"], D["female"]]
    perm = np.empty(N_PERM)
    lab = D["asd"].copy()
    r = np.random.default_rng(SEED + 7)
    for i in range(N_PERM):
        perm[i] = ols_t_first(y, [r.permutation(lab)] + base)
    p_perm = float((np.abs(perm) >= abs(t_obs) - 1e-12).sum() + 1) / (N_PERM + 1)
    return dict(
        label=label,
        n_asd=int(len(a)), n_ctrl=int(len(c)),
        mean_asd=float(a.mean()), mean_ctrl=float(c.mean()),
        sd_asd=float(a.std(ddof=1)), sd_ctrl=float(c.std(ddof=1)),
        delta=float(a.mean() - c.mean()), cohen_d=cohen_d(a, c),
        welch_t=float(stats.ttest_ind(a, c, equal_var=False).statistic),
        welch_p=float(stats.ttest_ind(a, c, equal_var=False).pvalue),
        mwu_p=float(stats.mannwhitneyu(a, c, alternative="two-sided").pvalue),
        ols_adj=m, ols_adj_depth=mdepth,
        perm_p_adj=p_perm, perm_n=N_PERM,
        perm_t_null_sd=float(np.nanstd(perm)),
    )


# ====================================================================== main
def main():
    t0 = time.time()
    sym, gene_total, cell_total, cell_nnz, eg, ec, ev, meta = load()
    R = {"seed": SEED, "n_perm": N_PERM, "n_null_sets": N_NULL_SETS}

    # ---------------- dataset description
    ucl, ccl = np.unique(meta["cluster"], return_counts=True)
    R["dataset"] = dict(
        path=str(VEL), n_nuclei=int(len(meta["cell"])), n_genes_matrix=int(len(sym)),
        clusters={str(a): int(b) for a, b in zip(ucl, ccl)},
        n_donors=int(len(set(meta["individual"].tolist()))),
        diagnosis={k: int((meta["diagnosis"] == k).sum()) for k in ["ASD", "Control"]},
    )

    # ---------------- build marker matrix once (all panels)
    wanted = sorted(set([g for v in STAGE_P6.values() for g in v] +
                        [g for v in STAGE_P10.values() for g in v] +
                        [g for v in ASTRO_PANEL.values() for g in v] +
                        [g for v in NEURON_PANEL.values() for g in v]))
    Xm, present = build_logcpm(sym, gene_total, cell_total, eg, ec, ev, wanted)
    ci = {g: i for i, g in enumerate(present)}
    R["markers"] = dict(requested=len(wanted), present=len(present),
                        missing=sorted(set(wanted) - set(present)),
                        p6_panel_coverage={s: [g for g in STAGE_P6[s] if g in ci]
                                           for s in ORDER5})
    print(f"[markers] {len(present)}/{len(wanted)} present", file=sys.stderr)

    # ================================================================ OL lineage
    ol = np.isin(meta["cluster"], LINEAGES["OL"])
    stage_ol, assign_ol = soft_stage(Xm[:, ol], ci, STAGE_P6, ORDER5)
    R["OL_cells"] = int(ol.sum())

    # --- clock internal validation inside Velmeshev
    cl_ol = meta["cluster"][ol]
    R["clock_validation"] = dict(
        mean_stage_by_cluster={c: float(stage_ol[cl_ol == c].mean()) for c in LINEAGES["OL"]},
        n_by_cluster={c: int((cl_ol == c).sum()) for c in LINEAGES["OL"]},
        argmax_stage_by_cluster={c: {ORDER5[k]: int((assign_ol[cl_ol == c] == k).sum())
                                     for k in range(5)} for c in LINEAGES["OL"]},
        separation_cohen_d=cohen_d(stage_ol[cl_ol == "Oligodendrocytes"],
                                   stage_ol[cl_ol == "OPC"]),
    )

    D, idx = donor_table(meta, ol)
    y_ol = donor_means(stage_ol, idx, D["donor"])
    ctrl = D["asd"] == 0
    R["clock_validation"]["age_corr_controls_spearman"] = dict(
        rho=float(stats.spearmanr(D["age"][ctrl], y_ol[ctrl]).statistic),
        p=float(stats.spearmanr(D["age"][ctrl], y_ol[ctrl]).pvalue), n=int(ctrl.sum()))
    R["clock_validation"]["age_corr_all_spearman"] = dict(
        rho=float(stats.spearmanr(D["age"], y_ol).statistic),
        p=float(stats.spearmanr(D["age"], y_ol).pvalue), n=int(len(y_ol)))

    R["donors"] = dict(
        n=int(len(D["donor"])), n_asd=int(D["asd"].sum()), n_ctrl=int((1 - D["asd"]).sum()),
        age_range_asd=[float(D["age"][D["asd"] == 1].min()), float(D["age"][D["asd"] == 1].max())],
        age_range_ctrl=[float(D["age"][ctrl].min()), float(D["age"][ctrl].max())],
        age_mean_asd=float(D["age"][D["asd"] == 1].mean()),
        age_mean_ctrl=float(D["age"][ctrl].mean()),
        age_welch_p=float(stats.ttest_ind(D["age"][D["asd"] == 1], D["age"][ctrl],
                                          equal_var=False).pvalue),
        female_asd=int(D["female"][D["asd"] == 1].sum()), female_ctrl=int(D["female"][ctrl].sum()),
        ol_cells_per_donor=dict(min=int(D["n"].min()), median=float(np.median(D["n"])),
                                max=int(D["n"].max())),
        ol_cells_asd=float(D["n"][D["asd"] == 1].mean()), ol_cells_ctrl=float(D["n"][ctrl].mean()),
        ncell_welch_p=float(stats.ttest_ind(D["n"][D["asd"] == 1], D["n"][ctrl],
                                            equal_var=False).pvalue),
    )

    # ---------------- PRIMARY
    R["PRIMARY_OL_stage"] = full_test(y_ol, D, "OL-lineage donor-mean soft stage (p6 panels)")

    # ---------------- composition vs within-class
    opc_m = cl_ol == "OPC"
    olig_m = cl_ol == "Oligodendrocytes"
    ind_ol = meta["individual"][ol]
    frac_opc = np.array([(opc_m[idx[d]]).mean() for d in D["donor"]])
    R["OPC_fraction"] = full_test(frac_opc, D, "per-donor OPC fraction of OL lineage")
    y_within_opc = np.array([stage_ol[idx[d]][opc_m[idx[d]]].mean()
                             if opc_m[idx[d]].sum() >= 20 else np.nan for d in D["donor"]])
    y_within_ol = np.array([stage_ol[idx[d]][olig_m[idx[d]]].mean()
                            if olig_m[idx[d]].sum() >= 20 else np.nan for d in D["donor"]])
    for nm, yy in [("WITHIN_OPCcluster_stage", y_within_opc),
                   ("WITHIN_Oligocluster_stage", y_within_ol)]:
        ok = np.isfinite(yy)
        Dk = {k: (v[ok] if isinstance(v, np.ndarray) else v) for k, v in D.items()}
        R[nm] = full_test(yy[ok], Dk, nm)
    # composition-adjusted primary
    R["PRIMARY_OL_stage_compadj"] = full_test(
        y_ol, D, "OL stage adjusted for OPC fraction", extra=[frac_opc],
        extra_names=["frac_OPC"])

    # ---------------- sensitivity: p10 panels
    stage_p10, _ = soft_stage(Xm[:, ol], ci, STAGE_P10, ORDER5)
    R["SENS_p10_panels"] = full_test(donor_means(stage_p10, idx, D["donor"]), D,
                                     "OL stage, p10 marker panels")

    # ---------------- sensitivity: per-nucleus depth residualisation
    dep = np.column_stack([np.ones(ol.sum()), np.log10(meta["umis"][ol]),
                           np.log10(meta["genes"][ol]), meta["mito"][ol], meta["ribo"][ol]])
    bb, *_ = np.linalg.lstsq(dep, stage_ol, rcond=None)
    stage_res = stage_ol - dep @ bb
    R["SENS_depth_residualised"] = full_test(donor_means(stage_res, idx, D["donor"]), D,
                                             "OL stage residualised on per-nucleus depth")
    R["SENS_depth_corr"] = dict(
        cell_stage_vs_logUMI=float(stats.pearsonr(stage_ol, np.log10(meta["umis"][ol]))[0]),
        cell_stage_vs_loggenes=float(stats.pearsonr(stage_ol, np.log10(meta["genes"][ol]))[0]),
        donor_stage_vs_mean_logUMI=float(stats.pearsonr(y_ol, D["mlogumi"])[0]),
        donor_stage_vs_log_ncells=float(stats.pearsonr(y_ol, np.log10(D["n"]))[0]),
        donor_stage_vs_pmi=float(stats.pearsonr(y_ol, D["pmi"])[0]),
        donor_stage_vs_rin=float(stats.pearsonr(y_ol, D["rin"])[0]),
        donor_stage_vs_fracPFC=float(stats.pearsonr(y_ol, D["fracPFC"])[0]),
    )

    # ---------------- sensitivity: equal cells per donor
    nmin = int(D["n"].min())
    r2 = np.random.default_rng(SEED + 11)
    ts, ds = [], []
    for _ in range(N_SUBSAMPLE):
        ysub = np.array([stage_ol[r2.choice(idx[d], nmin, replace=False)].mean()
                         for d in D["donor"]])
        ts.append(ols_t_first(ysub, [D["asd"], D["age"], D["female"]]))
        ds.append(cohen_d(ysub[D["asd"] == 1], ysub[ctrl]))
    R["SENS_equal_cells"] = dict(n_cells_per_donor=nmin, reps=N_SUBSAMPLE,
                                 t_mean=float(np.mean(ts)), t_sd=float(np.std(ts)),
                                 t_q=[float(q) for q in np.percentile(ts, [2.5, 50, 97.5])],
                                 d_mean=float(np.mean(ds)))

    # ---------------- sensitivity: region
    for reg in ["PFC", "ACC"]:
        m = ol & (meta["region"] == reg)
        Dr, idxr = donor_table(meta, m)
        st, _ = soft_stage(Xm[:, m], ci, STAGE_P6, ORDER5)
        yr = donor_means(st, idxr, Dr["donor"])
        R[f"SENS_region_{reg}"] = full_test(yr, Dr, f"OL stage, {reg} only")

    # ================================================================ CONTROL LINEAGES
    R["CONTROLS"] = {}
    for lin in ["AST", "EXC"]:
        m = np.isin(meta["cluster"], LINEAGES[lin])
        Dl, idxl = donor_table(meta, m)
        # (a) identical OL stage machinery in a non-OL lineage
        st, _ = soft_stage(Xm[:, m], ci, STAGE_P6, ORDER5)
        R["CONTROLS"][f"{lin}_OLpanel_stage"] = full_test(
            donor_means(st, idxl, Dl["donor"]), Dl,
            f"{lin}: OL stage panels applied verbatim (global-artefact control)")
        R["CONTROLS"][f"{lin}_n_cells"] = int(m.sum())
        # (b) lineage-analogous maturation score
        pan = ASTRO_PANEL if lin == "AST" else NEURON_PANEL
        order = sorted(pan.keys())
        st2, _ = soft_stage(Xm[:, m], ci, pan, order)
        R["CONTROLS"][f"{lin}_analogous_maturation"] = full_test(
            donor_means(st2, idxl, Dl["donor"]), Dl,
            f"{lin}: analogous {len(order)}-state maturation score")

    # ================================================================ AMBIENT-RNA INDEX
    # Myelin transcripts leaking into non-OL nuclei are a standard snRNA-seq artefact and
    # would move a marker-mean score in EVERY lineage at once. Index = donor-mean log-CPM of
    # the MOL panel measured in Microglia + Endothelial nuclei, where it should be ~absent.
    amb_m = np.isin(meta["cluster"], ["Microglia", "Endothelial"])
    mol_ix = [ci[g] for g in STAGE_P6["MOL"] if g in ci]
    amb_cell = Xm[mol_ix][:, amb_m].mean(axis=0)
    Da, idxa = donor_table(meta, amb_m)
    amb = donor_means(amb_cell, idxa, Da["donor"])
    order_a = {d: i for i, d in enumerate(Da["donor"])}
    amb_aligned = np.array([amb[order_a[d]] if d in order_a else np.nan for d in D["donor"]])
    R["AMBIENT"] = dict(
        n_cells=int(amb_m.sum()),
        test=full_test(amb_aligned[np.isfinite(amb_aligned)],
                       {k: (v[np.isfinite(amb_aligned)] if isinstance(v, np.ndarray) else v)
                        for k, v in D.items()},
                       "ambient MOL-panel signal in Microglia+Endothelial nuclei"),
        corr_with_OL_stage=float(stats.pearsonr(y_ol[np.isfinite(amb_aligned)],
                                                amb_aligned[np.isfinite(amb_aligned)])[0]),
    )
    ok = np.isfinite(amb_aligned)
    Dk = {k: (v[ok] if isinstance(v, np.ndarray) else v) for k, v in D.items()}
    R["PRIMARY_OL_stage_ambientadj"] = full_test(
        y_ol[ok], Dk, "OL stage adjusted for the ambient index",
        extra=[amb_aligned[ok]], extra_names=["ambient_MOL"])

    # ================================================================ AGE x DIAGNOSIS
    # POST-HOC / EXPLORATORY. Prompted by the clock validating on age in CONTROL donors
    # (rho=0.77) but not across all donors. Reported as a hypothesis, not a test result.
    def interaction(y, Dx, label):
        ag = (Dx["age"] - Dx["age"].mean())
        m = ols(y, [Dx["asd"] * ag, Dx["asd"], ag, Dx["female"]],
                ["ASDxage", "ASD", "age_c", "female"])
        r = np.random.default_rng(SEED + 23)
        perm = np.array([ols_t_first(y, [(lp := r.permutation(Dx["asd"])) * ag, lp, ag,
                                         Dx["female"]]) for _ in range(N_PERM)])
        return dict(label=label, ols=m,
                    perm_p=float((np.abs(perm) >= abs(m["t"][1]) - 1e-12).sum() + 1) /
                           (N_PERM + 1),
                    rho_ctrl=float(stats.spearmanr(Dx["age"][Dx["asd"] == 0],
                                                   y[Dx["asd"] == 0]).statistic),
                    rho_asd=float(stats.spearmanr(Dx["age"][Dx["asd"] == 1],
                                                  y[Dx["asd"] == 1]).statistic))
    R["POSTHOC_age_x_dx"] = {"OL_stage": interaction(y_ol, D, "OL stage")}
    for lin in ["AST", "EXC"]:
        m = np.isin(meta["cluster"], LINEAGES[lin])
        Dl, idxl = donor_table(meta, m)
        st, _ = soft_stage(Xm[:, m], ci, STAGE_P6, ORDER5)
        R["POSTHOC_age_x_dx"][f"{lin}_OLpanel"] = interaction(
            donor_means(st, idxl, Dl["donor"]), Dl, f"{lin} OL-panel score")
    R["POSTHOC_age_x_dx"]["OPC_fraction"] = interaction(frac_opc, D, "OPC fraction")
    R["POSTHOC_age_x_dx"]["within_Oligo"] = interaction(
        y_within_ol[np.isfinite(y_within_ol)],
        {k: (v[np.isfinite(y_within_ol)] if isinstance(v, np.ndarray) else v)
         for k, v in D.items()}, "within-Oligodendrocyte stage")
    R["POSTHOC_age_x_dx"]["ambient_index"] = interaction(amb_aligned[ok], Dk,
                                                         "ambient MOL index (control)")

    # interaction after adjusting for the ambient index and donor depth
    def interaction_adj(y, Dx, extra, extra_names, label):
        ag = Dx["age"] - Dx["age"].mean()
        return dict(label=label, ols=ols(
            y, [Dx["asd"] * ag, Dx["asd"], ag, Dx["female"]] + list(extra),
            ["ASDxage", "ASD", "age_c", "female"] + list(extra_names)))
    oko = np.isfinite(y_within_ol) & ok
    Do = {k: (v[oko] if isinstance(v, np.ndarray) else v) for k, v in D.items()}
    R["POSTHOC_age_x_dx"]["ADJ_OL_ambient_depth"] = interaction_adj(
        y_ol[ok], Dk, [amb_aligned[ok], Dk["mlogumi"], np.log10(Dk["n"]), Dk["rin"], Dk["pmi"]],
        ["ambient", "mean_log10UMI", "log10_ncells", "RIN", "PMI"],
        "OL stage x age, + ambient + depth + RIN + PMI")
    R["POSTHOC_age_x_dx"]["ADJ_withinOligo_ambient_depth"] = interaction_adj(
        y_within_ol[oko], Do,
        [amb_aligned[oko], Do["mlogumi"], np.log10(Do["n"]), Do["rin"], Do["pmi"]],
        ["ambient", "mean_log10UMI", "log10_ncells", "RIN", "PMI"],
        "within-Oligo stage x age, + ambient + depth + RIN + PMI")

    # leave-one-donor-out stability of the interaction
    def loo(y, Dx):
        ag = Dx["age"] - Dx["age"].mean()
        out = []
        for i in range(len(y)):
            m = np.ones(len(y), bool); m[i] = False
            agm = Dx["age"][m] - Dx["age"][m].mean()
            out.append((str(Dx["donor"][i]),
                        float(ols_t_first(y[m], [Dx["asd"][m] * agm, Dx["asd"][m], agm,
                                                 Dx["female"][m]]))))
        return out
    lo = loo(y_ol, D)
    lw = loo(y_within_ol[np.isfinite(y_within_ol)],
             {k: (v[np.isfinite(y_within_ol)] if isinstance(v, np.ndarray) else v)
              for k, v in D.items()})
    R["POSTHOC_age_x_dx"]["LOO"] = dict(
        OL_stage=dict(t=[x[1] for x in lo], worst_donor=min(lo, key=lambda x: abs(x[1]))[0],
                      t_range=[float(min(x[1] for x in lo)), float(max(x[1] for x in lo))],
                      n_sign_flips=int(sum(1 for x in lo if x[1] > 0))),
        within_Oligo=dict(t=[x[1] for x in lw],
                          worst_donor=min(lw, key=lambda x: abs(x[1]))[0],
                          t_range=[float(min(x[1] for x in lw)), float(max(x[1] for x in lw))],
                          n_sign_flips=int(sum(1 for x in lw if x[1] > 0))))

    # ================================================================ MATCHED-RANDOM NULL
    # Are arbitrary 5-panel softmax scores, built from expression-matched random genes,
    # just as "shifted" in ASD as the OL clock is?
    pool_rows = np.load(EXT / "velmeshev_subset.npz")["keep_rows"]
    gt = gene_total
    # expression-matched sampling: bucket by log10(total count) decile of the real markers
    real_rows = []
    for s in ORDER5:
        for g in STAGE_P6[s]:
            i = np.flatnonzero(sym == g)
            if i.size:
                real_rows.append(int(i[np.argmax(gt[i])]) + 1)
    real_rows = np.array(sorted(set(real_rows)))
    cand = np.setdiff1d(pool_rows, real_rows)
    cand = cand[gt[cand - 1] > 0]
    lc = np.log10(gt[cand - 1] + 1)
    sizes = [len([g for g in STAGE_P6[s] if g in ci]) for s in ORDER5]
    targets = []
    for s in ORDER5:
        tg = [np.log10(gt[np.flatnonzero(sym == g)].max() + 1)
              for g in STAGE_P6[s] if (sym == g).any()]
        targets.append(np.array(tg))
    r3 = np.random.default_rng(SEED + 13)
    null_t, null_d, null_ti = [], [], []
    age_c = D["age"] - D["age"].mean()
    Xol_cols = np.flatnonzero(ol)
    for j in range(N_NULL_SETS):
        chosen = []
        for k, s in enumerate(ORDER5):
            pick = []
            for tval in targets[k]:
                near = np.flatnonzero(np.abs(lc - tval) < 0.25)
                if near.size == 0:
                    near = np.array([int(np.argmin(np.abs(lc - tval)))])
                pick.append(int(cand[r3.choice(near)]))
            chosen.append(pick)
        Xr, pres = build_logcpm_rows(chosen, cell_total, eg, ec, ev, Xol_cols)
        if Xr is None:
            continue
        cir = {i: i for i in range(Xr.shape[0])}
        pan = {ORDER5[k]: pres[k] for k in range(5)}
        st, _ = soft_stage(Xr, cir, pan, ORDER5)
        yy = donor_means(st, idx, D["donor"])
        null_t.append(ols_t_first(yy, [D["asd"], D["age"], D["female"]]))
        null_d.append(cohen_d(yy[D["asd"] == 1], yy[ctrl]))
        null_ti.append(ols_t_first(yy, [D["asd"] * age_c, D["asd"], age_c, D["female"]]))
    null_t = np.array([x for x in null_t if np.isfinite(x)])
    null_ti = np.array([x for x in null_ti if np.isfinite(x)])
    ti_real = R["POSTHOC_age_x_dx"]["OL_stage"]["ols"]["t"][1]
    R["POSTHOC_age_x_dx"]["null_calibration"] = dict(
        n_sets=int(len(null_ti)), null_abs_t_q=[float(q) for q in
                                                np.percentile(np.abs(null_ti), [50, 90, 95, 99])],
        real_t=float(ti_real),
        p_vs_null=float((np.abs(null_ti) >= abs(ti_real)).sum() + 1) / (len(null_ti) + 1))
    t_real = R["PRIMARY_OL_stage"]["ols_adj"]["t"][1]
    R["MATCHED_RANDOM_NULL"] = dict(
        n_sets=int(len(null_t)),
        null_t_mean=float(null_t.mean()), null_t_sd=float(null_t.std()),
        null_abs_t_q=[float(q) for q in np.percentile(np.abs(null_t), [50, 90, 95, 99])],
        real_t=float(t_real),
        p_vs_null=float((np.abs(null_t) >= abs(t_real)).sum() + 1) / (len(null_t) + 1),
        null_d_mean=float(np.mean(null_d)), null_d_sd=float(np.std(null_d)),
    )

    # ================================================================ POWER
    n1 = R["donors"]["n_asd"]; n2 = R["donors"]["n_ctrl"]
    sd_pool = float(np.sqrt(((n1 - 1) * y_ol[D["asd"] == 1].var(ddof=1) +
                             (n2 - 1) * y_ol[ctrl].var(ddof=1)) / (n1 + n2 - 2)))
    R["POWER"] = dict(
        n_asd=n1, n_ctrl=n2,
        mde_d_80pct_unadjusted=mde_two_sample(n1, n2),
        mde_d_80pct_adjusted_age_sex=mde_two_sample(n1, n2, extra_df_loss=2),
        mde_d_50pct=mde_two_sample(n1, n2, power=0.50),
        observed_d=R["PRIMARY_OL_stage"]["cohen_d"],
        donor_sd_stage_units=sd_pool,
        mde_stage_units_80pct=mde_two_sample(n1, n2, extra_df_loss=2) * sd_pool,
        clock_dynamic_range_atlas2_stage_units=dict(
            fetal=0.51, yr2=1.91, yr13=2.35, adult=2.2,
            postnatal_2yr_to_13yr=0.44, note="Velmeshev ages 4-22 sit inside the flat "
                                             "postnatal segment of the clock"),
        power_for_atlas2_postnatal_slope=None,
    )
    # what power would we have for an effect the size of the clock's own 2yr->13yr step?
    d_ref = 0.44 / sd_pool
    df = n1 + n2 - 4
    crit = stats.t.ppf(0.975, df)
    nc = d_ref * np.sqrt(1 / (1 / n1 + 1 / n2))
    R["POWER"]["power_for_atlas2_postnatal_slope"] = dict(
        effect_stage_units=0.44, effect_as_d=float(d_ref),
        power=float(stats.nct.sf(crit, df, nc) + stats.nct.cdf(-crit, df, nc)))

    # ================================================================ persist
    (OUT / "p15_results.json").write_text(json.dumps(R, indent=2, default=float))
    # donor-level table, every number in the report traceable to this file
    hdr = ["donor", "diagnosis", "age", "sex", "n_OL_cells", "frac_OPC", "stage_OL",
           "stage_withinOPC", "stage_withinOligo", "stage_p10", "stage_depthres",
           "mean_log10UMI", "mean_log10genes", "pmi", "rin", "fracPFC"]
    y_p10 = donor_means(stage_p10, idx, D["donor"])
    y_res = donor_means(stage_res, idx, D["donor"])
    with open(OUT / "p15_donor_table.tsv", "w") as f:
        f.write("\t".join(hdr) + "\n")
        for i, d in enumerate(D["donor"]):
            f.write("\t".join(map(str, [
                d, "ASD" if D["asd"][i] else "Control", D["age"][i],
                "F" if D["female"][i] else "M", int(D["n"][i]), round(frac_opc[i], 4),
                round(y_ol[i], 4), round(y_within_opc[i], 4), round(y_within_ol[i], 4),
                round(y_p10[i], 4), round(y_res[i], 4), round(D["mlogumi"][i], 4),
                round(D["mloggene"][i], 4), D["pmi"][i], D["rin"][i],
                round(D["fracPFC"][i], 3)])) + "\n")
    np.save(OUT / "p15_null_t.npy", null_t)
    print(f"[done] {time.time()-t0:.0f}s -> {OUT}", file=sys.stderr)
    print(json.dumps({k: R[k] for k in ["PRIMARY_OL_stage", "MATCHED_RANDOM_NULL", "POWER"]},
                     indent=2, default=float)[:4000])


def build_logcpm_rows(chosen, cell_total, eg, ec, ev, cols):
    """Dense panel matrix for a matched-random gene set, restricted to `cols` cells."""
    flat = sorted({r for grp in chosen for r in grp})
    row2k = {r: k for k, r in enumerate(flat)}
    colmap = -np.ones(len(cell_total), dtype=np.int64)
    colmap[cols] = np.arange(len(cols))
    X = np.zeros((len(flat), len(cols)), dtype=np.float32)
    sel = np.isin(eg, np.array(flat, dtype=eg.dtype))
    gg, cc, vv = eg[sel], ec[sel], ev[sel].astype(np.float32)
    cj = colmap[cc - 1]
    keep = cj >= 0
    if keep.sum() == 0:
        return None, None
    kk = np.array([row2k[int(r)] for r in gg[keep]], dtype=np.int64)
    X[kk, cj[keep]] = vv[keep]
    tot = cell_total[cols].copy()
    tot[tot <= 0] = 1
    X = np.log1p(X / tot[None, :] * 1e4).astype(np.float32)
    pres = [[row2k[r] for r in grp] for grp in chosen]
    return X, pres


if __name__ == "__main__":
    main()

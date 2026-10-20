#!/usr/bin/env python3
"""
p15 STEP 3 -- assemble the traceable number table for the write-up.

Reads outputs/phase15/p15_results.json + p15_donor_table.tsv, writes:
  outputs/phase15/p15_summary_table.tsv   one row per donor-level test
  outputs/phase15/p15_extra.json          BH correction over the PRE-SPECIFIED family,
                                          confound correlations, power curve
No new modelling: every value here is read from, or a BH/arithmetic transform of, step 2.
"""
import json
from pathlib import Path
import numpy as np
from scipy import stats

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "outputs/phase15"
R = json.loads((OUT / "p15_results.json").read_text())

# ---- pre-specified family (committed in the task brief: primary + the named controls)
PRESPEC = [
    ("PRIMARY_OL_stage", R["PRIMARY_OL_stage"]),
    ("CONTROL_AST_OLpanel", R["CONTROLS"]["AST_OLpanel_stage"]),
    ("CONTROL_EXC_OLpanel", R["CONTROLS"]["EXC_OLpanel_stage"]),
    ("CONTROL_AST_analogous", R["CONTROLS"]["AST_analogous_maturation"]),
    ("CONTROL_EXC_analogous", R["CONTROLS"]["EXC_analogous_maturation"]),
]
EXTRA = [
    ("DECOMP_OPC_fraction", R["OPC_fraction"]),
    ("DECOMP_within_OPCcluster", R["WITHIN_OPCcluster_stage"]),
    ("DECOMP_within_Oligocluster", R["WITHIN_Oligocluster_stage"]),
    ("ADJ_composition", R["PRIMARY_OL_stage_compadj"]),
    ("ADJ_ambient", R["PRIMARY_OL_stage_ambientadj"]),
    ("SENS_p10_panels", R["SENS_p10_panels"]),
    ("SENS_depth_residualised", R["SENS_depth_residualised"]),
    ("SENS_region_PFC", R["SENS_region_PFC"]),
    ("SENS_region_ACC", R["SENS_region_ACC"]),
    ("CONTROL_ambient_index", R["AMBIENT"]["test"]),
]

rows = []
for tag, d in PRESPEC + EXTRA:
    o = d["ols_adj"]
    i = o["names"].index("ASD")
    rows.append(dict(test=tag, family="prespecified" if (tag, d) in PRESPEC else "secondary",
                     n_asd=d["n_asd"], n_ctrl=d["n_ctrl"],
                     mean_asd=d["mean_asd"], mean_ctrl=d["mean_ctrl"], delta=d["delta"],
                     cohen_d=d["cohen_d"], welch_p=d["welch_p"], mwu_p=d["mwu_p"],
                     beta_adj=o["beta"][i], t_adj=o["t"][i], p_adj=o["p"][i], df=o["df"],
                     perm_p=d["perm_p_adj"]))


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    prev = 1.0
    for k in range(m - 1, -1, -1):
        prev = min(prev, p[o[k]] * m / (k + 1))
        q[o[k]] = prev
    return q


pre = [r for r in rows if r["family"] == "prespecified"]
qs = bh([r["perm_p"] for r in pre])
for r, q in zip(pre, qs):
    r["BH_q_prespec"] = float(q)
for r in rows:
    r.setdefault("BH_q_prespec", "")

cols = ["test", "family", "n_asd", "n_ctrl", "mean_asd", "mean_ctrl", "delta", "cohen_d",
        "welch_p", "mwu_p", "beta_adj", "t_adj", "p_adj", "df", "perm_p", "BH_q_prespec"]
with open(OUT / "p15_summary_table.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(f"{r[c]:.5g}" if isinstance(r[c], float) else str(r[c])
                          for c in cols) + "\n")

# ---- donor-level confound correlations straight off the donor table
T = [l.rstrip("\n").split("\t") for l in open(OUT / "p15_donor_table.tsv")]
H = {k: i for i, k in enumerate(T[0])}
D = T[1:]


def col(name, f=float):
    return np.array([f(r[H[name]]) if r[H[name]] != "nan" else np.nan for r in D])


asd = np.array([r[H["diagnosis"]] == "ASD" for r in D])
age, n, fopc = col("age"), col("n_OL_cells"), col("frac_OPC")
stg, rin, pmi = col("stage_OL"), col("rin"), col("pmi")


def pr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return dict(r=float(stats.pearsonr(a[m], b[m])[0]), p=float(stats.pearsonr(a[m], b[m])[1]),
                n=int(m.sum()))


extra = dict(
    confounds=dict(
        fracOPC_vs_log_ncells=pr(fopc, np.log10(n)),
        fracOPC_vs_age_controls=pr(fopc[~asd], age[~asd]),
        fracOPC_vs_age_asd=pr(fopc[asd], age[asd]),
        stage_vs_log_ncells=pr(stg, np.log10(n)),
        stage_vs_RIN=pr(stg, rin), stage_vs_PMI=pr(stg, pmi),
        ncells_vs_age=pr(n, age),
    ),
    bh_prespecified={r["test"]: r["BH_q_prespec"] for r in pre},
)

# ---- power curve for the donor-level design actually in hand
n1, n2 = R["donors"]["n_asd"], R["donors"]["n_ctrl"]
df = n1 + n2 - 4
crit = stats.t.ppf(0.975, df)
sd = R["POWER"]["donor_sd_stage_units"]
curve = []
for d in [0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5]:
    nc = d * np.sqrt(1 / (1 / n1 + 1 / n2))
    curve.append(dict(cohen_d=d, stage_units=round(d * sd, 3),
                      power=float(stats.nct.sf(crit, df, nc) + stats.nct.cdf(-crit, df, nc))))
extra["power_curve_alpha05_adj"] = curve
extra["power_note"] = (
    f"n_ASD={n1}, n_ctrl={n2}, model stage ~ dx + age + sex (df={df}). "
    f"Donor-level pooled SD of the OL stage score = {sd:.3f} stage units.")
(OUT / "p15_extra.json").write_text(json.dumps(extra, indent=2))
print(json.dumps(extra, indent=2))

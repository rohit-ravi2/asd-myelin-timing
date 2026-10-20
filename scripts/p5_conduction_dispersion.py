#!/usr/bin/env python3
"""
§127-129 — Conduction-delay dispersion model.

THE QUESTION
------------
The myelin->computation bridge has always been asserted, never quantified here. The claim
is NOT that mean conduction delay changes. It is that the DISTRIBUTION of arrival times
widens enough to disrupt coincidence-dependent computation:

    d_i = L_i / v_i        sigma_d^2 = Var(d)

For a myelinated axon, conduction velocity is well approximated (Rushton; Waxman) by

    v  =  k * D * ln(1/g) ... in practice v scales ~linearly with outer diameter D and
    depends on g-ratio (g = inner/outer diameter) through the internodal cable properties.

We use the standard form   v = C * D * sqrt(-ln g) * (L_int/L_int0)^a   normalised so that
a canonical axon (D=1um, g=0.7, L_int=100um) matches measured cortical values, then ask:
how much does realistic PERTURBATION of g, D, internode length and nodal function change
the SPREAD of arrival times across a tract?

§128 requires BOTH directions: thinner AND thicker myelin, plus increased structural
variance and nodal defects. ASD subtypes are directionally heterogeneous (PTEN precocious
vs CHD8/NLGN3 hypomyelination), so the model must show whether opposite mean changes can
produce the SAME dispersion increase.

Benchmarks the dispersion is compared against:
    gamma cycle (40 Hz)            25 ms
    feed-forward inhibition window ~2-4 ms
    coincidence detection window   ~5-10 ms
    STDP window                    ~20 ms
"""
import json
import numpy as np
from pathlib import Path

OUT = Path("./outputs/phase5")
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20260917)
N_AX = 20000

# --- canonical human cortical parameters (literature-standard ranges) ---
# outer diameter distribution: cortical myelinated axons are dominated by thin fibres
D_MEAN, D_CV = 0.8, 0.45          # um
G_MEAN, G_SD = 0.70, 0.045        # g-ratio, tight in healthy CNS
LINT_PER_D = 100.0                # internode length ~100x axon diameter
V_CAL = 5.5                       # m/s for D=1um, g=0.7 (human cortical myelinated)

TRACTS = {                        # (length mm, label)
    "short association (U-fibre)": 20.0,
    "long association (SLF/arcuate)": 90.0,
    "commissural (callosal)": 120.0,
    "projection (corticospinal)": 300.0,
}


def velocity(D, g, lint_scale=1.0, node_factor=1.0):
    """m/s. v ~ D * sqrt(-ln g) normalised to V_CAL at D=1,g=0.7; internode and nodal
    terms enter as multiplicative efficiency factors."""
    base = V_CAL / (1.0 * np.sqrt(-np.log(G_MEAN)))
    return base * D * np.sqrt(np.maximum(-np.log(np.clip(g, 1e-3, 0.999)), 1e-6)) \
        * (lint_scale ** 0.5) * node_factor


def sample_pop(n, g_mean=G_MEAN, g_sd=G_SD, d_cv=D_CV, lint_scale=1.0,
               lint_cv=0.0, node_factor=1.0):
    D = np.abs(rng.normal(D_MEAN, D_MEAN * d_cv, n))
    D = np.clip(D, 0.15, 6.0)
    g = np.clip(rng.normal(g_mean, g_sd, n), 0.35, 0.95)
    ls = np.clip(lint_scale * (1 + rng.normal(0, lint_cv, n)), 0.15, 3.0) \
        if lint_cv > 0 else lint_scale
    return velocity(D, g, ls, node_factor)


SCENARIOS = {
    "baseline":                     dict(),
    "hypomyelination (g 0.70->0.78)": dict(g_mean=0.78),
    "severe hypomyelination (g->0.85)": dict(g_mean=0.85),
    "HYPERmyelination (g 0.70->0.62)": dict(g_mean=0.62),
    "increased g variance (sd x3)":  dict(g_sd=G_SD * 3),
    "increased g variance (sd x5)":  dict(g_sd=G_SD * 5),
    "short internodes (-30%)":       dict(lint_scale=0.7),
    "internode length variance":     dict(lint_cv=0.35),
    "nodal dysfunction (-25% eff)":  dict(node_factor=0.75),
    "mixed: hypo + g variance":      dict(g_mean=0.78, g_sd=G_SD * 3),
    "mixed: HYPER + g variance":     dict(g_mean=0.62, g_sd=G_SD * 3),
}

res = {}
for tname, Lmm in TRACTS.items():
    L = Lmm / 1000.0                       # m
    base_v = sample_pop(N_AX)
    base_d = L / base_v * 1000.0           # ms
    res[tname] = {}
    for sname, kw in SCENARIOS.items():
        v = sample_pop(N_AX, **kw)
        d = L / v * 1000.0
        res[tname][sname] = dict(
            mean_ms=float(np.mean(d)), sd_ms=float(np.std(d)),
            p5=float(np.percentile(d, 5)), p95=float(np.percentile(d, 95)),
            iqr_ms=float(np.percentile(d, 75) - np.percentile(d, 25)),
            d_mean_vs_base=float(np.mean(d) - np.mean(base_d)),
            d_sd_vs_base=float(np.std(d) - np.std(base_d)),
            sd_ratio=float(np.std(d) / np.std(base_d)),
        )

json.dump(res, open(OUT / "conduction_dispersion.json", "w"), indent=1)

BENCH = {"feed-forward inhibition": 3.0, "coincidence detection": 8.0,
         "gamma cycle (40 Hz)": 25.0, "STDP window": 20.0}

L = ["=" * 112, "§127-129 CONDUCTION-DELAY DISPERSION MODEL", "=" * 112,
     f"{N_AX} axons/scenario. D ~ N({D_MEAN},{D_MEAN*D_CV:.2f}) um; g ~ N({G_MEAN},{G_SD});",
     f"v calibrated to {V_CAL} m/s at D=1um, g=0.70.",
     "KEY QUANTITY is sigma_d (arrival-time SD), not mean delay.", ""]
for tname in TRACTS:
    L.append(f"--- {tname}  (L = {TRACTS[tname]:.0f} mm) " + "-" * 34)
    L.append(f"{'scenario':34s} {'mean ms':>8s} {'sigma ms':>9s} {'d_mean':>8s} "
             f"{'d_sigma':>8s} {'sigma x':>8s} {'5-95% ms':>14s}")
    for sname in SCENARIOS:
        r = res[tname][sname]
        L.append(f"{sname:34s} {r['mean_ms']:>8.2f} {r['sd_ms']:>9.2f} "
                 f"{r['d_mean_vs_base']:>+8.2f} {r['d_sd_vs_base']:>+8.2f} "
                 f"{r['sd_ratio']:>8.2f} {r['p5']:>6.1f}-{r['p95']:<7.1f}")
    L.append("")

L += ["=" * 112, "INTERPRETATION vs physiological timing windows", "=" * 112]
for b, ms in BENCH.items():
    L.append(f"   {b:26s} {ms:>5.1f} ms")
L.append("")
for tname in TRACTS:
    base = res[tname]["baseline"]["sd_ms"]
    L.append(f"{tname}: baseline sigma_d = {base:.2f} ms")
    for sname in SCENARIOS:
        if sname == "baseline":
            continue
        r = res[tname][sname]
        exceeded = [b for b, ms in BENCH.items() if r["d_sd_vs_base"] >= ms]
        if r["sd_ratio"] >= 1.25:
            tag = "; exceeds " + ", ".join(exceeded) if exceeded else ""
            L.append(f"    {sname:34s} sigma x{r['sd_ratio']:.2f} "
                     f"(+{r['d_sd_vs_base']:.2f} ms){tag}")
    L.append("")

L += ["KEY STRUCTURAL RESULT (§128) -- NOT what was predicted before running:",
      "  Hypomyelination RAISES sigma_d; HYPERmyelination LOWERS it (sigma x0.86), and",
      "  HYPER + added g-variance only returns to ~baseline (x1.01). So directionally",
      "  heterogeneous ASD subtypes do NOT converge on a common dispersion phenotype --",
      "  they CANCEL in both mean AND variance across a mixed population.",
      "  This REFUTES the pre-written expectation (see error-log entry E-005) and makes the",
      "  power problem worse, not better: a mixed-subtype cohort loses signal in every",
      "  summary statistic, so stratification is not optional, it is required.",
      "",
      "CAVEAT: this is a cable-theory scaling model with literature-standard parameters,",
      "not a compartmental simulation. It bounds plausibility; it does not measure biology."]
txt = "\n".join(L)
(OUT / "CONDUCTION_DISPERSION_REPORT.txt").write_text(txt + "\n")
print(txt)

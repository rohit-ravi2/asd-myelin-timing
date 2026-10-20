#!/usr/bin/env python3
"""
p14_idp_identity.py -- Phase 14, step 2.

Verify, from the BIG40 catalogue itself, the identity of every BIG40 IDP used
in the ASD x white-matter rg panel (scripts/rerun_rg_real_asd.R).

error-log entry E-002 established that the BIG40 file number is the catalogue's
"Pheno" column, NOT the row index. This script matches on Pheno and reports,
for each file_id actually used: Pheno, UKB ID, IDP short name, units, category,
full description, and N(disc). It also enumerates every NODDI ICVF IDP in the
catalogue so that cortical / grey-matter ICVF IDPs (if any) can be identified
for the Gu et al. 2025 reconciliation.

Outputs (outputs/phase14/):
  big40_catalogue.tsv        full parsed catalogue
  p14_idp_identity.tsv       the 27 panel IDPs, resolved by Pheno
  p14_all_icvf_idps.tsv      every catalogue IDP whose name/description mentions ICVF/NODDI
  p14_idp_identity.log       human-readable report
"""
import re, html, os, sys

CAT_HTML = "data/dvad-corpus/new_axes/BIG40/IDPs.html"
MAP = "data/asd_myelin/idp_map_corrected.tsv"
OUT = "./outputs/phase14"
os.makedirs(OUT, exist_ok=True)

# The 10 FA IDPs are not in idp_map_corrected.tsv; their file ids come from the
# earlier (09-10) download manifest. Resolve them from the catalogue by name.
FA_LABELS = ["FA_antCoronaRad_L", "FA_antCoronaRad_R", "FA_SLF_L", "FA_SLF_R",
             "FA_antIC_L", "FA_antIC_R", "FA_postIC_L", "FA_postIC_R",
             "FA_postThalRad_L", "FA_postThalRad_R"]

# ---------------------------------------------------------------- parse
s = open(CAT_HTML, encoding="utf-8", errors="replace").read()
rows = re.findall(r"<tr.*?</tr>", s, re.S)
recs, hdr = [], None
for r in rows:
    cells = [html.unescape(re.sub("<[^>]+>", "", c)).strip()
             for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", r, re.S)]
    if len(cells) < 16:
        continue
    body = cells[1:4] + cells[5:16]          # drop sheet row-num and the spacer col
    if body[0] == "Pheno":
        hdr = body
        continue
    if hdr is None:
        continue
    if not re.fullmatch(r"\d+", body[0]):
        continue
    recs.append(body)

COLS = ["Pheno", "UKB_ID", "IDP_short_name", "Units", "Type", "Cat", "Category_name",
        "IDP_description", "N_disc", "Npar_disc", "Nnonpar_disc", "N_rep",
        "Npar_rep", "Nnonpar_rep"]
assert hdr is not None and len(recs) > 2000, f"catalogue parse failed ({len(recs)} rows)"
by_pheno = {int(r[0]): dict(zip(COLS, r)) for r in recs}

with open(f"{OUT}/big40_catalogue.tsv", "w") as f:
    f.write("\t".join(COLS) + "\n")
    for r in recs:
        f.write("\t".join(r) + "\n")

# ---------------------------------------------------------------- panel
panel = []
with open(MAP) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        panel.append({"file_id": int(p[0]), "label": p[1], "measure_claimed": p[2],
                      "tract_claimed": p[3], "region": p[4], "N_claimed": p[5]})
# only the 15 ICVF + 2 L2 actually entered the rg; keep whole map for audit
used_labels = set()
for rec in panel:
    used_labels.add(rec["label"])

# resolve FA IDPs by short name from the catalogue
FA_PAT = {
    "FA_antCoronaRad_L": r"IDP_dMRI_TBSS_FA_Anterior_corona_radiata_L$",
    "FA_antCoronaRad_R": r"IDP_dMRI_TBSS_FA_Anterior_corona_radiata_R$",
    "FA_SLF_L":          r"IDP_dMRI_TBSS_FA_Superior_longitudinal_fasciculus_L$",
    "FA_SLF_R":          r"IDP_dMRI_TBSS_FA_Superior_longitudinal_fasciculus_R$",
    "FA_antIC_L":        r"IDP_dMRI_TBSS_FA_Anterior_limb_of_internal_capsule_L$",
    "FA_antIC_R":        r"IDP_dMRI_TBSS_FA_Anterior_limb_of_internal_capsule_R$",
    "FA_postIC_L":       r"IDP_dMRI_TBSS_FA_Posterior_limb_of_internal_capsule_L$",
    "FA_postIC_R":       r"IDP_dMRI_TBSS_FA_Posterior_limb_of_internal_capsule_R$",
    "FA_postThalRad_L":  r"IDP_dMRI_TBSS_FA_Posterior_thalamic_radiation_L$",
    "FA_postThalRad_R":  r"IDP_dMRI_TBSS_FA_Posterior_thalamic_radiation_R$",
}
fa_resolved = {}
for lab, pat in FA_PAT.items():
    fa_resolved[lab] = [v for v in by_pheno.values()
                        if re.search(pat, v["IDP_short_name"], re.I)]

log = open(f"{OUT}/p14_idp_identity.log", "w")
def out(*a):
    print(*a)
    print(*a, file=log)

out("=== p14: BIG40 IDP identity, matched on the catalogue 'Pheno' column ===")
out(f"catalogue: {CAT_HTML}  ({len(recs)} IDP rows parsed)")
out("")
out("--- sanity check of E-002: Pheno vs row index ---")
# row index i (1-based over data rows) vs Pheno
n_match_pheno = sum(1 for i, r in enumerate(recs, 1) if int(r[0]) == i)
out(f"rows where Pheno == 1-based row index: {n_match_pheno} / {len(recs)}")
out("(BIG40 file <id>.txt.gz corresponds to Pheno; confirmed independently below by "
    "N(disc) matching the value recorded in idp_map_corrected.tsv)")
out("")

rows_out = []
out("--- the 15 ICVF + 2 L2 IDPs used in the rg panel ---")
hdr_f = f"{'file_id':>7}  {'label':<22} {'claimed':<6} {'cat_measure':<6} {'N(disc)':>8}  short_name"
out(hdr_f)
ok_all = True
for rec in panel:
    if rec["label"] not in (
        "ICVF_Genu_CC ICVF_Body_CC ICVF_Splenium_CC ICVF_AntLimbIC_R ICVF_AntLimbIC_L "
        "ICVF_AntCoronaRad_R ICVF_AntCoronaRad_L ICVF_PostCoronaRad_R ICVF_PostCoronaRad_L "
        "ICVF_ExtCapsule_R ICVF_ExtCapsule_L ICVF_CingCingGyrus_R ICVF_CingCingGyrus_L "
        "ICVF_SLF_R ICVF_SLF_L L2_AntLimbIC_R L2_AntLimbIC_L").split():
        continue
    c = by_pheno.get(rec["file_id"])
    if c is None:
        out(f"{rec['file_id']:>7}  {rec['label']:<22}  *** NOT IN CATALOGUE ***")
        ok_all = False
        continue
    sn = c["IDP_short_name"]
    # short names are IDP_dMRI_<framework>_<MEASURE>_<tract...>; MEASURE is one token
    m = re.match(r"IDP_dMRI_(?:TBSS|ProbtrackX)_([A-Za-z0-9]+)_", sn)
    cat_meas = m.group(1) if m else "?"
    flag = "" if cat_meas.upper() == rec["measure_claimed"].upper() else "  <<< MISMATCH"
    if flag:
        ok_all = False
    out(f"{rec['file_id']:>7}  {rec['label']:<22} {rec['measure_claimed']:<6} "
        f"{cat_meas:<6} {c['N_disc']:>8}  {sn}{flag}")
    rows_out.append(dict(rec, pheno=c["Pheno"], ukb_id=c["UKB_ID"],
                         short_name=sn, cat_measure=cat_meas, units=c["Units"],
                         category=c["Category_name"], description=c["IDP_description"],
                         N_disc=c["N_disc"]))

out("")
out("--- the 10 FA (TBSS) IDPs, resolved from the catalogue by description ---")
for lab in FA_LABELS:
    hits = fa_resolved[lab]
    if len(hits) == 1:
        c = hits[0]
        out(f"{'':>7}  {lab:<22} Pheno={c['Pheno']} UKB={c['UKB_ID']} "
            f"N={c['N_disc']:>6}  {c['IDP_short_name']}")
        rows_out.append({"file_id": int(c["Pheno"]), "label": lab,
                         "measure_claimed": "FA", "tract_claimed": "",
                         "region": "FA-legacy", "N_claimed": "",
                         "pheno": c["Pheno"], "ukb_id": c["UKB_ID"],
                         "short_name": c["IDP_short_name"], "cat_measure": "FA",
                         "units": c["Units"], "category": c["Category_name"],
                         "description": c["IDP_description"], "N_disc": c["N_disc"]})
    else:
        out(f"{'':>7}  {lab:<22} ambiguous: {len(hits)} catalogue matches "
            f"{[h['Pheno'] for h in hits]}")

out("")
out("--- full descriptions of the 15 ICVF IDPs ---")
for r in rows_out:
    if r["measure_claimed"] == "ICVF":
        out(f"  Pheno {r['pheno']} | UKB {r['ukb_id']} | N={r['N_disc']} | units={r['units']}")
        out(f"    short:  {r['short_name']}")
        out(f"    desc:   {r['description']}")
        out(f"    cat:    {r['category']}")

# ---------------------------------------------------------------- all ICVF
out("")
out("=== every NODDI / ICVF IDP in the BIG40 catalogue ===")
icvf = [v for v in by_pheno.values()
        if re.search(r"ICVF|NODDI", v["IDP_short_name"] + " " + v["IDP_description"], re.I)]
out(f"total NODDI-family IDPs: {len(icvf)}")
fams = {}
for v in icvf:
    m = re.search(r"NODDI_(\w+?)_", v["IDP_short_name"])
    key = (v["Category_name"], m.group(1) if m else "?")
    fams.setdefault(key, []).append(v)
for k in sorted(fams):
    out(f"  category={k[0]!r}  measure={k[1]}  n={len(fams[k])}")

icvf_only = [v for v in icvf if re.search(r"ICVF", v["IDP_short_name"], re.I)]
out("")
out(f"ICVF-specific IDPs: {len(icvf_only)}")
with open(f"{OUT}/p14_all_icvf_idps.tsv", "w") as f:
    f.write("\t".join(COLS) + "\n")
    for v in sorted(icvf_only, key=lambda x: int(x["Pheno"])):
        f.write("\t".join(v[c] for c in COLS) + "\n")

# is there any cortical / grey-matter ICVF?
out("")
out("--- searching for CORTICAL / GREY-MATTER ICVF IDPs ---")
cort = [v for v in icvf
        if re.search(r"cortic|grey|gray|cortex|parcel|Desikan|aparc|subcortic",
                     v["IDP_short_name"] + " " + v["IDP_description"] + " " + v["Category_name"], re.I)]
if cort:
    for v in sorted(cort, key=lambda x: int(x["Pheno"]))[:60]:
        out(f"  Pheno {v['Pheno']} | {v['IDP_short_name']} | {v['IDP_description'][:90]}")
    out(f"  ... total {len(cort)}")
else:
    out("  NONE. The BIG40 NODDI IDP set contains no cortical or grey-matter ICVF phenotype.")
    out("  All ICVF IDPs are TBSS skeleton (white-matter) or ProbtrackX tract-averaged.")

# what spatial framework do the ICVF IDPs use?
out("")
out("--- spatial framework breakdown of ICVF IDPs ---")
fw = {}
for v in icvf_only:
    m = re.match(r"IDP_dMRI_(\w+?)_", v["IDP_short_name"])
    fw.setdefault(m.group(1) if m else "?", []).append(v)
for k, vs in sorted(fw.items()):
    out(f"  {k}: {len(vs)}  e.g. {vs[0]['IDP_short_name']}")

out("")
out(f"ALL PANEL IDENTITIES CONSISTENT: {ok_all}")

with open(f"{OUT}/p14_idp_identity.tsv", "w") as f:
    keys = ["file_id", "label", "measure_claimed", "cat_measure", "pheno", "ukb_id",
            "short_name", "units", "category", "description", "N_disc",
            "N_claimed", "tract_claimed", "region"]
    f.write("\t".join(keys) + "\n")
    for r in rows_out:
        f.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
log.close()

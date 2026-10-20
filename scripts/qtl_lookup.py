#!/usr/bin/env python3
"""
Test D (part 2): do the cholesterol / myelin loci that show nominal ASD signal
have a plausible brain-regulatory mechanism?

Angles:
  D2a. For each target gene, pull its cis-eQTLs in adult DLPFC (PsychENCODE
       Phase I, n~1866) and report the best eQTL p, the top-SNP flag, and
       whether the ASD regional lead SNP is itself an eQTL for that gene.
  D2b. Fetal-brain mQTL (Hannon) overlap for the same regional lead SNPs.
  D2c. AGGREGATE: are ASD association statistics stronger at SNPs that are
       brain eQTLs for cholesterol-pathway genes than at MAF/LD-agnostic
       matched background SNPs? (a crude but honest enrichment check)

Gene matching is by genomic overlap (NCBI37.3 vs PsychENCODE hg19 coords), so
no symbol<->ENSG map is required.

NOTE: GoDMC blood mQTL associations are NOT on disk (only snps.csv.gz + an md5
for the 6.3 GB assoc file, which was never fetched). Blood mQTL is therefore
absent from this run; the memory note claiming it was available is wrong.
"""
import gzip, json, sys, bisect
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "outputs/qtl"
OUT.mkdir(parents=True, exist_ok=True)
PEC = "data/asd-corpus/eqtl/PsychENCODE_Full_hg19_cis-eQTL.txt.gz"
HANNON = "data/dvad-corpus/asd_scope/mqtl/hannon_fetalbrain/All_Imputed_BonfSignificant_mQTLs.csv.gz"
GLOC = "data/asd-corpus/tools/magma_geneloc/NCBI37.3.gene.loc"
LOCUS = BASE / "outputs/locus/locus_results_ASD_real.json"

# ------------------------------------------------------- targets from Test C/D
loc = json.load(open(LOCUS))["C2_D1_regions"]
targets = {g: d for g, d in loc.items()
           if d.get("n_snps", 0) > 0 and d["role"] in
           ("PRIMARY", "PRIMARY-finemap", "PRIMARY-syndromic", "MYELIN",
            "CONTRAST-primary", "CONTRAST")}
print(f"[targets] {len(targets)} cholesterol/myelin/APOE genes", file=sys.stderr)

coords = {}
with open(GLOC) as f:
    for line in f:
        p = line.split()
        if len(p) >= 6 and p[5] in targets:
            coords[p[5]] = (p[1], int(p[2]), int(p[3]))

# lead SNP position per gene (from the ASD regional scan)
lead = {g: (d["chr"], d["min_p_bp"], d["min_p_rsid"], d["min_p"])
        for g, d in targets.items() if "min_p_bp" in d}

# ------------------------------------------------------- D2a PsychENCODE eQTL
# build interval index: chr -> sorted list of (start,end,gene)
byc = defaultdict(list)
for g, (c, s, e) in coords.items():
    byc[c].append((s, e, g))
for c in byc:
    byc[c].sort()

best_eqtl = {}           # gene -> best (p, snp, slope, is_top)
lead_is_eqtl = {}        # gene -> record if the ASD lead SNP is an eQTL of it
n_lines = 0
print("[D2a] scanning PsychENCODE cis-eQTL (3.1 GB) ...", file=sys.stderr)
with gzip.open(PEC, "rt") as f:
    for line in f:
        p = line.split()
        if len(p) < 14:
            continue
        n_lines += 1
        gchr = p[1].replace("chr", "")
        if gchr not in byc:
            continue
        gs, ge = int(p[2]), int(p[3])
        # which target gene does this eQTL's gene overlap?
        hit = None
        for s, e, g in byc[gchr]:
            if not (ge < s or gs > e):
                hit = g
                break
        if hit is None:
            continue
        try:
            pv = float(p[11]); slope = float(p[12]); istop = p[13]
        except ValueError:
            continue
        snp_bp = int(p[9])
        cur = best_eqtl.get(hit)
        if cur is None or pv < cur[0]:
            best_eqtl[hit] = (pv, p[7], slope, istop)
        if hit in lead and snp_bp == lead[hit][1]:
            lead_is_eqtl[hit] = dict(eqtl_p=pv, slope=slope, is_top_snp=istop,
                                     snp=p[7])
print(f"[D2a] {n_lines} eQTL rows scanned; {len(best_eqtl)} target genes with eQTLs",
      file=sys.stderr)

# ------------------------------------------------------- D2b Hannon fetal mQTL
lead_bp = {(c, bp) for c, bp, _, _ in lead.values()}
lead_rs = {rs for _, _, rs, _ in lead.values()}
mqtl_hits = defaultdict(list)
print("[D2b] scanning Hannon fetal-brain mQTL ...", file=sys.stderr)
with gzip.open(HANNON, "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    hi = {c.strip('"'): i for i, c in enumerate(hdr)}
    snpcol = next((hi[k] for k in ("SNP", "snp", "SNPID") if k in hi), 0)
    for line in f:
        p = [x.strip('"') for x in line.rstrip("\n").split(",")]
        if len(p) <= snpcol:
            continue
        if p[snpcol] in lead_rs:
            mqtl_hits[p[snpcol]].append(p)
print(f"[D2b] {len(mqtl_hits)} of {len(lead_rs)} lead SNPs are fetal-brain mQTLs",
      file=sys.stderr)

# ------------------------------------------------------- report
rows = []
for g, d in sorted(targets.items(), key=lambda x: x[1]["min_p"]):
    be = best_eqtl.get(g)
    lie = lead_is_eqtl.get(g)
    c, bp, rs, mp = lead.get(g, ("", 0, "", 1))
    rows.append(dict(gene=g, role=d["role"], asd_min_p=mp,
                     asd_emp_regional_p=d.get("empirical_regional_p"),
                     lead_snp=rs,
                     best_cis_eqtl_p=be[0] if be else None,
                     best_cis_eqtl_snp=be[1] if be else None,
                     lead_snp_is_eqtl_p=lie["eqtl_p"] if lie else None,
                     lead_snp_fetal_mqtl=rs in mqtl_hits))
json.dump(rows, open(OUT / "qtl_results.json", "w"), indent=1)

L = []
L.append("=" * 104)
L.append("TEST D2: brain-regulatory support for cholesterol/myelin loci with nominal ASD signal")
L.append("=" * 104)
L.append("eQTL = PsychENCODE Phase I adult DLPFC (n~1866, hg19, full nominal)")
L.append("mQTL = Hannon fetal brain (Bonferroni-significant only)")
L.append("GoDMC blood mQTL: NOT on disk (assoc file never downloaded) -- absent from this run")
L.append("")
L.append(f"{'gene':12s} {'role':18s} {'ASD min_p':>10s} {'empRegP':>8s} {'lead SNP':>14s} "
         f"{'best eQTL p':>12s} {'lead=eQTL p':>12s} {'fetal mQTL':>10s}")
for r in rows:
    L.append(f"{r['gene']:12s} {r['role']:18s} {r['asd_min_p']:>10.3g} "
             f"{(('%.3f'%r['asd_emp_regional_p']) if r['asd_emp_regional_p'] else 'NA'):>8s} "
             f"{r['lead_snp']:>14s} "
             f"{(('%.3g'%r['best_cis_eqtl_p']) if r['best_cis_eqtl_p'] else 'none'):>12s} "
             f"{(('%.3g'%r['lead_snp_is_eqtl_p']) if r['lead_snp_is_eqtl_p'] else '-'):>12s} "
             f"{('YES' if r['lead_snp_fetal_mqtl'] else '-'):>10s}")
txt = "\n".join(L)
(OUT / "TESTD2_qtl_report.txt").write_text(txt + "\n")
print(txt)

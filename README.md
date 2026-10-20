# asd-myelin-timing

Analysis code and result tables for testing whether autism genetics points at myelin, on public
GWAS, single-cell, methylation and imaging-genetics data. The starting hypothesis (weak lipid
supply leaves late-myelinating frontal cortex under-insulated) did not hold up, but it was
replaced by a sharper one about developmental timing. The main results: autism liability is
genetically correlated with lower neurite density in all 15 white-matter tracts tested (mean
rg = -0.10) with no frontal gradient, the common-variant signal localises to neurons rather than
the oligodendrocyte lineage once gene constraint is controlled for, and a conduction model shows
that under- and over-myelination push arrival-time spread in opposite directions.

The write-up is at https://www.rohitravi.com/projects/asd-myelin-timing.

## Results

Every number here is in a file under `outputs/`.

**Autism x white matter (LDSC rg).** Grove et al. 2019 autism GWAS against the 15 BIG40 UK
Biobank ICVF (neurite density) tracts: every rg is negative, mean -0.102. Anterior and posterior
tracts differ by 0.011, so there is no frontal gradient. Best tract is the left cingulum
(rg = -0.142, SE 0.053, p = 0.0073); BH over the 15 tracts gives a best q of 0.064, so nothing is
significant tract by tract. An allergic-disease GWAS run through the same pipeline is the
negative control. `outputs/phase14/p14_signfix_results.csv` (27 IDPs incl. FA/L2; the `q` column
there is BH over all 27), `outputs/phase14/p14_headtohead.csv` compares the tract estimates with
Gu et al. 2025.

**Where on the oligodendrocyte maturation axis does autism signal sit?** Each gene gets an
expression-weighted lineage stage in a human infant PFC atlas (GSE294786, 3,021 OL-lineage cells,
13,644 genes) and MAGMA gene Z is regressed on it.

| model | t | p |
|---|---|---|
| naive | +2.88 | 0.0040 |
| + LOEUF constraint, gene length, expression breadth | +0.47 | 0.638 |
| stage-label permutation (1,000) | | 0.70 |
| neuron specificity, same conditioning (positive control) | +3.36 | 0.00079 |

The naive effect is a gene-property confound; the same test finds the neuronal signal. A second
atlas (GSE155488, 10,391 cells, fetal to 13 yr, includes mature OL) is null on all three
maturation scores after conditioning (p = 0.53 / 0.81 / 0.87, permutation p = 0.81).
`outputs/phase1/P1_PSEUDOTIME_REPORT.txt`, `outputs/phase6/ATLAS2_REPLICATION_REPORT.txt`.

**Rare de novo variants, carrier level** (denovo-db SSC, 2,140 probands vs 1,614 unaffected
siblings, each module with a synonymous negative control and a length/LOEUF-matched null):

- neurogenic-to-gliogenic switch set, LoF: OR 17.35, p = 3.3e-5. Its synonymous control is
  OR 6.03, so I do not count it.
- chromatin regulators, LoF: OR 2.41, p = 1.6e-7, q = 3.1e-6, clean controls.
- cholesterol, lipid-transport, myelin and APOE modules: null.

`outputs/phase1/P1_RARE_BURDEN_V2_REPORT.txt`.

**Cytokine receptors in infant PFC.** Mean log-CPM per cell type: IL6R 0.811 in microglia vs
0.015 in OPCs (about 54x); IL1R1 is near zero everywhere; complement receptors are absent from
OPCs. The adenosine receptor ADORA1 goes the other way (0.137 OPC vs 0.015 microglia) and peaks
at the committed-precursor stage. `outputs/phase4/CYTOKINE_RECEPTOR_REPORT.txt`,
`outputs/phase10/STAGE_TRAJECTORY_REPORT.txt`.

**Conduction-delay spread.** 20,000 simulated axons per scenario with realistic diameter and
g-ratio distributions. Moderate hypomyelination (g 0.70 to 0.78) adds 5.33 ms of arrival-time
SD on a 90 mm association tract and 1.02 ms on a 20 mm U-fibre. Hypermyelination reduces the
spread (-2.84 ms at 90 mm), so a cohort mixing both directions cancels in variance as well as
mean. `outputs/phase5/CONDUCTION_DISPERSION_REPORT.txt`.

**Cell-type methylation in autistic PFC** (GSE298243). OPCs have the highest rate of
diagnosis-associated DMRs of the eight cell types (272 of 12,586 tested bins, 2.16%). Matching
on estimator standard error raises the excess from OR 1.92 to 2.86 (95% CI 2.50-3.29), so it is
not a noise artifact. It can't be settled as a cell-type comparison because each cell type is
tested on a nearly disjoint set of regions (2,103 of 124,166 features shared by any pair).
`outputs/phase7/PFC_DMR_REPORT.txt`, `outputs/phase16/opc_se_matched.json`, `outputs/p16.log`.

**Maturation clock on autistic tissue.** The OL-lineage stage score built on GSE155488, applied
to 25,494 OL-lineage nuclei from Velmeshev et al. 2019 (31 donors, 15 autism / 16 control).
In controls alone, mean stage tracks age at rho = 0.77. Autism vs control, adjusted for age and
sex: p = 0.158, and the study is underpowered for the effect sizes in play.
`outputs/phase15/p15_results.json`, `outputs/phase15/p15_summary_table.tsv`.

## Repo layout

```
scripts/        analysis code, one script per test (Python, R, bash)
annot/          gene-set definitions (MAGMA .set files + manifest with sources)
outputs/        small result tables, reports and JSON behind the numbers above
  phase1/       pseudotime test, rare-variant burden
  phase4-5/     cytokine receptors, gate competition, conduction model
  phase6-7/     second-atlas replication, PFC methylation
  phase10-16/   stage trajectories, MR, WM factor, rg sign check, clock, OPC DMR checks
  magma/        MAGMA gene-set analysis summaries
  sldsc/results partitioned-heritability results per gene set
```

Script prefixes follow the order the tests were run in (`p1_` ... `p16_`). The unprefixed
scripts are the earlier gene-set, MAGMA and apoptosis/oligodendrocyte tests.

## How to run

Python 3.12 with numpy, pandas, scipy and openpyxl; R with GenomicSEM for the rg scripts;
LDSC (python 2 env, set `LDSC_PYTHON`), MAGMA 1.10 and PLINK 1.9 on the path. Data go under
`data/` in the repo root (or point `DATA` at another directory) and scripts are run from the
repo root:

```
python scripts/build_genesets.py           # gene sets -> annot/
bash   scripts/run_magma.sh                # MAGMA gene + gene-set analysis
Rscript scripts/p14_sign_fix.R             # ASD x ICVF rg with the BIG40 allele fix
python scripts/p1_pseudotime.py            # maturation-axis test
python scripts/p1_rare_burden_v2.py        # de novo carrier burden
python scripts/p5_conduction_dispersion.py # conduction model, no external data
python scripts/p15_extract_velmeshev.py && python scripts/p15_clock_velmeshev.py
```

The subdirectory names under `data/` are the ones the scripts expect (`asd-corpus/`,
`dvad-corpus/asd_scope/`, `ldsc_ref/`); check the path constants at the top of each script.

## Data

None of the input data is in this repo. All of it is public:

- Autism GWAS: Grove et al. 2019 (iPSYCH-PGC), from the PGC download page; Matoba et al. 2020
  (SPARK + iPSYCH + PGC).
- White-matter GWAS: BIG40 UK Biobank IDP summary statistics, https://open.oxcin.ox.ac.uk/ukbiobank/big40/
  (note: a2 is the effect allele; `scripts/p14_sign_fix.R` handles this).
- LDSC references (eur_w_ld_chr, baselineLD v2.2, 1000G EUR Phase 3) from the LDSC resources;
  MAGMA g1000_eur reference and NCBI37.3 gene locations from the MAGMA site.
- Single-cell: GEO GSE294786 (infant PFC), GSE155488 (OL lineage, fetal to adolescent),
  Velmeshev et al. 2019 autism cortex snRNA-seq.
- Methylation: GEO GSE298243 (PFC single-nucleus methylation).
- Rare variants: denovo-db, SSC samples.
- Also used: BrainSpan, MitoCarta 3.0, PsychENCODE cis-eQTL, Hannon et al. fetal-brain mQTL, Ensembl 115 GFF3,
  Gu et al. 2025 (Mol Psychiatry) supplementary tables.

## Notes

- Several outputs were re-run after I found the BIG40 IDP sumstats had been munged with the
  reference allele as A1, which flips the sign of every rg. The corrected run is
  `p14_sign_fix.R` / `p14_signfix_results.csv`; earlier rg tables are not included here.
- The phase-11 bidirectional MR is included for completeness but pleiotropy makes it weak
  evidence either way.

## License

MIT

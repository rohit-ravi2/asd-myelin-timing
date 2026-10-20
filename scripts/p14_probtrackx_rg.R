## p14_probtrackx_rg.R -- Phase 14, step 5.
##
## THE DISCRIMINATING TEST.
##
## Gu et al. 2025 (PMID 39994426) supplementary table S4 reports rg(autism, ICVF)
## for 15 AutoPtx white-matter tracts: 14/15 NEGATIVE, mean ~ -0.074, none
## surviving FDR. This project reports 15/15 POSITIVE, mean +0.102, on 15 TBSS
## FA-skeleton JHU-label tracts. Same nominal quantity, opposite sign.
##
## Three things differ between the two analyses:
##   (i)   the autism GWAS        -- ours iPSYCH-PGC Grove 2019 (18,381/27,969);
##                                   theirs iPSYCH-only (19,870/39,078);
##   (ii)  the ICVF GWAS          -- ours BIG40 (Smith 2021, UKB, fastGWA-free
##                                   pipeline); theirs their own fastGWA on
##                                   UKB+ABCD meta (36,663);
##   (iii) the tract definition   -- ours TBSS FA-skeleton JHU ROIs;
##                                   theirs AutoPtx probabilistic tractography.
##
## This script removes (iii) as an explanation. BIG40 also publishes ProbtrackX
## ICVF IDPs computed on exactly Gu's AutoPtx tracts (Pheno 1950-1976). We run
## OUR ASD GWAS against BIG40 ProbtrackX slf_l/r, cgc_l/r, atr_l -- Gu's own
## tract definitions -- and compare tract-for-tract.
##
##   If ProbtrackX comes out NEGATIVE  -> the tract definition (TBSS skeleton vs
##       tractography) drives the discrepancy, and this project's positive is an
##       artifact of the TBSS skeleton ROI family.
##   If ProbtrackX comes out POSITIVE  -> tract definition is NOT the cause; the
##       discrepancy lives in the autism GWAS or the ICVF GWAS pipeline, and the
##       two results are estimating the same thing with opposite signs.

suppressMessages(library(GenomicSEM))
SRC <- "data/dvad-corpus/asd_myelin/big40_probtrackx"
RG  <- "./outputs/rg"
OUT <- "./outputs/phase14"
LD  <- "data/ldsc_ref/eur_w_ld_chr/"
hm3 <- "data/ldsc_ref/eur_w_ld_chr/w_hm3.snplist"
dir.create(file.path(OUT, "munge"), recursive = TRUE, showWarnings = FALSE)
setwd(OUT)

MAP <- data.frame(
  id    = c(1971, 1972, 1954, 1955, 1952),
  label = c("pICVF_slf_l", "pICVF_slf_r", "pICVF_cgc_l", "pICVF_cgc_r", "pICVF_atr_l"),
  gu_tract = c("slf", "slf", "cgc", "cgc", "atr"),
  stringsAsFactors = FALSE)
## Gu et al. 2025 supplementary S4, white-matter tract block (bilateral tracts)
GU <- data.frame(tract = c("slf", "cgc", "atr"),
                 gu_rg = c(-0.0999, -0.0639, -0.0885),
                 gu_se = c( 0.0564,  0.0591,  0.0531),
                 stringsAsFactors = FALSE)
## this project's TBSS-skeleton estimates for the homologous tracts
TBSS <- data.frame(label = c("pICVF_slf_l", "pICVF_slf_r", "pICVF_cgc_l", "pICVF_cgc_r"),
                   tbss_label = c("ICVF_SLF_L", "ICVF_SLF_R",
                                  "ICVF_CingCingGyrus_L", "ICVF_CingCingGyrus_R"),
                   stringsAsFactors = FALSE)

sink(file.path(OUT, "p14_probtrackx_rg.log"), split = TRUE)
cat("=== p14 step 5: ASD x BIG40 ProbtrackX (AutoPtx) ICVF -- Gu's tract definitions ===\n")
cat("date:", format(Sys.time()), "\n\n")

## ---- prep: BIG40 cols are chr rsid pos a1 a2 beta se pval(-log10) ----------
## ALLELE ORIENTATION (error-log entry E-007): BIG40 documents a1 as the REFERENCE
## allele and a2 as the EFFECT allele -- "the sign of regression beta relates to
## the a2 allele count". munge/LDSC attribute BETA to the declared A1. So the
## munge input must be A1 = source a2 ($5), A2 = source a1 ($4). The original
## prep_idps.sh used A1 = $4 and therefore inverted the sign of every BIG40 IDP.
for (i in seq_len(nrow(MAP))) {
  f  <- file.path(SRC, paste0(MAP$id[i], ".txt.gz"))
  mu <- file.path(OUT, "munge", paste0(MAP$label[i], ".txt"))
  stopifnot(file.exists(f))
  if (!file.exists(mu)) {
    cmd <- sprintf(
      "zcat %s | awk 'BEGIN{OFS=\"\\t\"; print \"SNP\",\"A1\",\"A2\",\"BETA\",\"SE\",\"P\"} NR>1 && $2 ~ /^rs/ && $7+0>0 { lp=$8+0; if(lp>300) lp=300; print $2, toupper($5), toupper($4), $6, $7, 10^(-lp) }' > %s",
      f, mu)
    system(cmd)
  }
  cat(sprintf("prepped %s (%s rows)\n", MAP$label[i],
              format(as.numeric(system(sprintf("wc -l < %s", mu), intern = TRUE)), big.mark = ",")))
}

## ---- munge (same N as the TBSS panel: BIG40 stats33k) ---------------------
todo <- MAP$label[!file.exists(file.path(OUT, paste0(MAP$label, ".sumstats.gz")))]
if (length(todo)) {
  munge(files = file.path(OUT, "munge", paste0(todo, ".txt")), hm3 = hm3,
        trait.names = todo, N = rep(31355, length(todo)),
        info.filter = 0.9, maf.filter = 0.01)
}

## ---- LDSC against the same ASD sumstats -----------------------------------
file.copy(file.path(RG, "ASDreal.sumstats.gz"), file.path(OUT, "ASDreal.sumstats.gz"),
          overwrite = FALSE)
tn <- c("ASDreal", MAP$label)
L <- ldsc(traits = paste0(tn, ".sumstats.gz"),
          sample.prev = c(18381/46350, rep(NA, nrow(MAP))),
          population.prev = c(0.012, rep(NA, nrow(MAP))),
          ld = LD, wld = LD, trait.names = tn, stand = TRUE)
saveRDS(L, file.path(OUT, "p14_probtrackx_ldsc.rds"))

S <- L$S_Stand; colnames(S) <- rownames(S) <- tn
k <- length(tn)
idx <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
sev <- sqrt(diag(L$V_Stand)); se <- matrix(NA, k, k)
for (i in seq_len(nrow(idx))) {
  se[idx[i,1], idx[i,2]] <- sev[i]; se[idx[i,2], idx[i,1]] <- sev[i]
}

res <- data.frame(label = MAP$label, gu_tract = MAP$gu_tract,
                  rg = S["ASDreal", MAP$label],
                  se = se[1, match(MAP$label, tn)], stringsAsFactors = FALSE)
res$z <- res$rg / res$se
res$p <- 2 * pnorm(-abs(res$z))
res <- merge(res, GU, by.x = "gu_tract", by.y = "tract", all.x = TRUE)
old <- read.csv(file.path(RG, "RG_REAL_ASD_results.csv"), stringsAsFactors = FALSE)
res <- merge(res, TBSS, by = "label", all.x = TRUE)
res$tbss_rg <- old$rg[match(res$tbss_label, old$label)]
res$tbss_se <- old$se[match(res$tbss_label, old$label)]

cat("\n--- rg(ASD Grove2019, BIG40 ProbtrackX ICVF) vs Gu S4 vs our TBSS ---\n")
print(res[, c("label", "rg", "se", "z", "p", "gu_rg", "gu_se", "tbss_label",
              "tbss_rg", "tbss_se")], row.names = FALSE, digits = 3)

cat(sprintf("\nProbtrackX mean rg = %+.4f  (%d/%d positive)\n",
            mean(res$rg), sum(res$rg > 0), nrow(res)))
cat(sprintf("matched TBSS mean rg = %+.4f\n", mean(res$tbss_rg, na.rm = TRUE)))
cat(sprintf("Gu S4 mean rg for these tracts = %+.4f\n", mean(unique(res[, c("gu_tract","gu_rg")])$gu_rg)))

## per-tract difference vs Gu, treating the two as independent (conservative:
## both ICVF GWAS draw on UK Biobank, so errors are positively correlated and
## the true |z| is LARGER than reported here)
res$d_vs_gu <- res$rg - res$gu_rg
res$z_vs_gu <- res$d_vs_gu / sqrt(res$se^2 + res$gu_se^2)
res$p_vs_gu <- 2 * pnorm(-abs(res$z_vs_gu))
cat("\n--- difference from Gu (conservative, assumes independence) ---\n")
print(res[, c("label", "rg", "gu_rg", "d_vs_gu", "z_vs_gu", "p_vs_gu")],
      row.names = FALSE, digits = 3)

if (!is.null(L$I)) {
  I <- L$I; colnames(I) <- rownames(I) <- tn
  cat("\nunivariate intercepts (expect ~1):\n"); print(round(diag(I), 4))
  cat("cross-trait intercepts vs ASD (expect ~0):\n"); print(round(I["ASDreal", MAP$label], 4))
}

cat("\nVERDICT INPUT:\n")
if (mean(res$rg) < 0) {
  cat("  ProbtrackX rg is NEGATIVE -> tract definition (TBSS skeleton vs AutoPtx\n")
  cat("  tractography) explains the disagreement; this project's TBSS positive does\n")
  cat("  not generalise to the tract family Gu used.\n")
} else {
  cat("  ProbtrackX rg is POSITIVE -> tract definition does NOT explain the\n")
  cat("  disagreement. The discrepancy must live in the autism GWAS (Grove2019\n")
  cat("  iPSYCH-PGC vs iPSYCH-2015) or in the ICVF GWAS pipeline (BIG40 vs the\n")
  cat("  authors' own fastGWA UKB+ABCD meta-analysis).\n")
}
write.csv(res, file.path(OUT, "p14_probtrackx_results.csv"), row.names = FALSE)
sink()
cat("\ndone\n")

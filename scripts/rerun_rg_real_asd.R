## Test 1 + 1b RE-RUN on the REAL ASD GWAS.
##
## The 2026-09-10 and 09-14 runs used data/dvad-corpus/asd_myelin/ASD_Grove2019.gz,
## which is NOT autism -- it is an allergic-disease/asthma GWAS (top loci TSLP, IL33,
## IL1RL1, ORMDL3/GSDMB, SMAD3, TLR1, EMSY/LRRC32, MHC-II, STAT6; N=242,569; Grove's
## five real ASD leads all null in it). Every rg, h2 and gene-window result from those
## runs is VOID.
##
## This re-runs both against the real iPSYCH-PGC_ASD_Nov2017 (validated: rs910805
## p=2.041e-09 matches publication; 9,112,387 rows matches GCST007556 snpCount).
## N = 46,350 (18,381 ca / 27,969 co), sample.prev 0.397, pop.prev 0.012.

suppressMessages(library(GenomicSEM))
OLD <- "data/asd_myelin"
NEW <- "./outputs/rg"
dir.create(NEW, recursive = TRUE, showWarnings = FALSE)
setwd(NEW)
hm3 <- "data/ldsc_ref/eur_w_ld_chr/w_hm3.snplist"
LD  <- "data/ldsc_ref/eur_w_ld_chr/"

## ---- munge the real ASD GWAS -------------------------------------------------
if (!file.exists("munge_ASDreal.txt")) {
  cat("[prep] reading real ASD sumstats\n")
  d <- data.table::fread("data/asd-corpus/genetics/ASD_Grove2019_REAL.gz")
  d <- d[grepl("^rs", SNP) & P > 0 & is.finite(OR) & OR > 0]
  out <- data.frame(SNP = d$SNP, A1 = toupper(d$A1), A2 = toupper(d$A2),
                    BETA = log(d$OR), SE = d$SE, P = d$P)
  data.table::fwrite(out, "munge_ASDreal.txt", sep = "\t")
  cat("[prep] wrote", nrow(out), "SNPs\n")
}
if (!file.exists("ASDreal.sumstats.gz")) {
  munge(files = "munge_ASDreal.txt", hm3 = hm3, trait.names = "ASDreal",
        N = 46350, info.filter = 0.9, maf.filter = 0.01)
}

## ---- assemble the IDP panel (already munged, reuse) --------------------------
map <- read.delim(file.path(OLD, "idp_map_corrected.tsv"), stringsAsFactors = FALSE)
fa <- c("FA_antCoronaRad_L","FA_antCoronaRad_R","FA_SLF_L","FA_SLF_R","FA_antIC_L",
        "FA_antIC_R","FA_postIC_L","FA_postIC_R","FA_postThalRad_L","FA_postThalRad_R")
idps <- c(map$label, fa)
src  <- file.path(OLD, paste0(idps, ".sumstats.gz"))
keep <- file.exists(src)
idps <- idps[keep]; src <- src[keep]
for (i in seq_along(src)) if (!file.exists(basename(src[i]))) file.symlink(src[i], basename(src[i]))
cat("[panel]", length(idps), "IDPs\n")

## ---- LDSC --------------------------------------------------------------------
tn <- c("ASDreal", idps)
L <- ldsc(traits = paste0(tn, ".sumstats.gz"),
          sample.prev = c(18381/46350, rep(NA, length(idps))),
          population.prev = c(0.012, rep(NA, length(idps))),
          ld = LD, wld = LD, trait.names = tn, stand = TRUE)
saveRDS(L, "rg_real_asd_ldsc.rds")

rg <- L$S_Stand; colnames(rg) <- rownames(rg) <- tn
k <- length(tn)
idx <- which(lower.tri(rg, diag = TRUE), arr.ind = TRUE)
sev <- sqrt(diag(L$V_Stand)); semat <- matrix(NA, k, k)
for (i in seq_len(nrow(idx))) {
  semat[idx[i,1], idx[i,2]] <- sev[i]; semat[idx[i,2], idx[i,1]] <- sev[i]
}
meas <- ifelse(idps %in% fa, "FA", map$measure[match(idps, map$label)])
reg  <- ifelse(idps %in% fa, "FA-legacy", map$region_class[match(idps, map$label)])
res <- data.frame(label = idps, measure = meas, region = reg,
                  rg = rg["ASDreal", idps], se = semat[1, match(idps, tn)])
res$z <- res$rg / res$se
res$p <- 2 * pnorm(-abs(res$z))
res$p_bonf <- pmin(1, res$p * nrow(res))
res <- res[order(res$p), ]
write.csv(res, "RG_REAL_ASD_results.csv", row.names = FALSE)

sink("RG_REAL_ASD_REPORT.txt")
cat("=== rg(REAL ASD Grove2019, myelin-sensitive + FA dMRI IDPs) ===\n")
cat("date:", format(Sys.time()), "\n")
cat("SUPERSEDES Test 1 (09-10) and Test 1b (09-14), both of which used an\n")
cat("allergic-disease GWAS mislabelled as ASD. Those results are VOID.\n\n")
h2 <- L$S[1,1]; h2se <- sqrt(L$V[1,1])
cat(sprintf("ASD h2 (liability, prev 0.012) = %.4f (SE %.4f, Z = %.2f)\n",
            h2, h2se, h2/h2se))
cat("  [published Grove2019 SNP-h2 ~ 0.118; the old run's 0.348 was the allergy trait]\n\n")
print(res, row.names = FALSE, digits = 3)
icvf <- res[res$measure == "ICVF", ]
a <- icvf[icvf$region == "anterior", ]; p_ <- icvf[icvf$region == "posterior", ]
if (nrow(a) && nrow(p_)) {
  cat(sprintf("\nICVF anterior mean rg  = %+.4f (n=%d)\n", mean(a$rg), nrow(a)))
  cat(sprintf("ICVF posterior mean rg = %+.4f (n=%d)\n", mean(p_$rg), nrow(p_)))
  cat(sprintf("anterior - posterior   = %+.4f  (descriptive; tracts correlated)\n",
              mean(a$rg) - mean(p_$rg)))
}
cat(sprintf("\nMedian SE = %.4f -> 80%% power for |rg| >= %.3f at alpha=0.05\n",
            median(res$se), 2.8 * median(res$se)))
cat(sprintf("Bonferroni (%d tests) 80%% power for |rg| >= %.3f\n", nrow(res),
            (qnorm(1 - 0.05/(2*nrow(res))) + qnorm(0.8)) * median(res$se)))
sink()
cat("\nwrote RG_REAL_ASD_REPORT.txt\n")
print(res, row.names = FALSE, digits = 3)

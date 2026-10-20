## p14_sign_fix.R -- Phase 14, the correction.
##
## error-log entry E-007. BIG40 documents its allele convention verbatim:
##
##   "The alleles are defined such that a1 is the reference allele, and a2 is the
##    alternative allele ... The effect allele (in the GWAS linear model) is a2,
##    meaning that the sign of regression beta relates to the a2 allele count."
##
## LDSC / GenomicSEM munge attribute the signed statistic to the DECLARED A1.
## A faithful munge input therefore needs A1 = source a2, A2 = source a1.
##
## `prep_idps.sh` (2026-09-14) wrote `print $2, toupper($4), toupper($5), $6, ...`
## i.e. A1 = source a1 = the REFERENCE allele. Every BIG40 IDP sumstats file used
## by this project therefore has its effect direction INVERTED, and every
## cross-trait rg against ASD has the WRONG SIGN. (IDP-vs-IDP rg is unaffected:
## both sides flip, so the product is unchanged -- which is exactly why the
## internal FA/ICVF/L2 sign control in p14_rerun_rg.R passed.)
##
## The sister project proved this empirically on BIG40 IDP 1437 (WMH): the a1
## orientation gives rg(WMH, small-vessel stroke) = -0.395, the a2 orientation
## +0.395, and only the latter is biologically possible.
## See the BIG40 documentation on effect alleles
##
## This script re-preps all 27 BIG40 IDPs with A1 = source a2, re-munges, and
## re-runs LDSC against the same validated ASD GWAS.

suppressMessages(library(GenomicSEM))
OLD  <- "data/asd_myelin"
OUT  <- "./outputs/phase14"
WORK <- "data/dvad-corpus/asd_myelin/p14_signfix"
LD   <- "data/ldsc_ref/eur_w_ld_chr/"
hm3  <- "data/ldsc_ref/eur_w_ld_chr/w_hm3.snplist"
dir.create(file.path(WORK, "munge"), recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SRC_EXTRA <- "data/dvad-corpus/asd_myelin/big40_extra"
SRC_FA    <- "data/dvad-corpus/new_axes/BIG40"

map <- read.delim(file.path(OLD, "idp_map_corrected.tsv"), stringsAsFactors = FALSE)
map <- map[map$measure %in% c("ICVF", "L2") &
           map$label %in% c(grep("^ICVF_", map$label, value = TRUE),
                            "L2_AntLimbIC_L", "L2_AntLimbIC_R"), ]
fa <- data.frame(
  file_id = c(1475, 1474, 1493, 1492, 1469, 1468, 1471, 1470, 1481, 1480),
  label = c("FA_antCoronaRad_L", "FA_antCoronaRad_R", "FA_SLF_L", "FA_SLF_R",
            "FA_antIC_L", "FA_antIC_R", "FA_postIC_L", "FA_postIC_R",
            "FA_postThalRad_L", "FA_postThalRad_R"),
  measure = "FA", tract = "", region_class = "FA-legacy", N = 31356,
  stringsAsFactors = FALSE)
all_idp <- rbind(map, fa)
all_idp$src <- ifelse(all_idp$file_id %in% fa$file_id,
                      file.path(SRC_FA,    paste0(all_idp$file_id, ".txt.gz")),
                      file.path(SRC_EXTRA, paste0(all_idp$file_id, ".txt.gz")))

setwd(WORK)
sink(file.path(OUT, "p14_sign_fix.log"), split = TRUE)
cat("=== p14 SIGN FIX: BIG40 A1 = source a2 (the documented effect allele) ===\n")
cat("date:", format(Sys.time()), "\n\n")

ok <- file.exists(all_idp$src)
cat(sprintf("IDPs with raw source on disk: %d / %d\n", sum(ok), nrow(all_idp)))
if (any(!ok)) print(all_idp[!ok, c("file_id", "label", "src")])
all_idp <- all_idp[ok, ]

## ---- re-prep with the CORRECT orientation ---------------------------------
for (i in seq_len(nrow(all_idp))) {
  mu <- file.path(WORK, "munge", paste0(all_idp$label[i], ".txt"))
  if (file.exists(mu) && file.info(mu)$size > 1e8) next
  cmd <- sprintf(
    "zcat %s | awk 'BEGIN{OFS=\"\\t\"; print \"SNP\",\"A1\",\"A2\",\"BETA\",\"SE\",\"P\"} NR>1 && $2 ~ /^rs/ && $7+0>0 { lp=$8+0; if(lp>300) lp=300; print $2, toupper($5), toupper($4), $6, $7, 10^(-lp) }' > %s",
    all_idp$src[i], mu)
  system(cmd)
  cat(sprintf("[%2d/%2d] prepped %-22s (A1 = source a2)\n", i, nrow(all_idp), all_idp$label[i]))
  flush.console()
}

## ---- munge -----------------------------------------------------------------
todo <- all_idp$label[!file.exists(file.path(WORK, paste0(all_idp$label, ".sumstats.gz")))]
if (length(todo)) {
  sub <- all_idp[all_idp$label %in% todo, ]
  munge(files = file.path(WORK, "munge", paste0(sub$label, ".txt")), hm3 = hm3,
        trait.names = sub$label, N = sub$N, info.filter = 0.9, maf.filter = 0.01)
}
file.copy("./outputs/rg/ASDreal.sumstats.gz",
          file.path(WORK, "ASDreal.sumstats.gz"), overwrite = FALSE)

## ---- LDSC ------------------------------------------------------------------
tn <- c("ASDreal", all_idp$label)
L <- ldsc(traits = paste0(tn, ".sumstats.gz"),
          sample.prev = c(18381/46350, rep(NA, nrow(all_idp))),
          population.prev = c(0.012, rep(NA, nrow(all_idp))),
          ld = LD, wld = LD, trait.names = tn, stand = TRUE)
saveRDS(L, file.path(OUT, "p14_signfix_ldsc.rds"))

S <- L$S_Stand; colnames(S) <- rownames(S) <- tn
k <- length(tn)
idx <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
sev <- sqrt(diag(L$V_Stand)); se <- matrix(NA, k, k)
for (i in seq_len(nrow(idx))) {
  se[idx[i,1], idx[i,2]] <- sev[i]; se[idx[i,2], idx[i,1]] <- sev[i]
}

res <- data.frame(label = all_idp$label, measure = all_idp$measure,
                  region = all_idp$region_class,
                  rg = S["ASDreal", all_idp$label],
                  se = se[1, match(all_idp$label, tn)], stringsAsFactors = FALSE)
res$z <- res$rg / res$se
res$p <- 2 * pnorm(-abs(res$z))
old <- read.csv("./outputs/rg/RG_REAL_ASD_results.csv",
                stringsAsFactors = FALSE)
res$rg_uncorrected <- old$rg[match(res$label, old$label)]
res$sum_check <- res$rg + res$rg_uncorrected     # should be ~0 if a pure sign flip
res <- res[order(res$p), ]

cat("\n--- CORRECTED rg(ASD, BIG40 IDPs), A1 = documented effect allele a2 ---\n")
print(res, row.names = FALSE, digits = 3)

icvf <- res[res$measure == "ICVF", ]
cat(sprintf("\nICVF: n=%d, mean rg = %+.4f, negative = %d/%d\n",
            nrow(icvf), mean(icvf$rg), sum(icvf$rg < 0), nrow(icvf)))
cat(sprintf("FA  : mean rg = %+.4f\n", mean(res$rg[res$measure == "FA"])))
cat(sprintf("L2  : mean rg = %+.4f\n", mean(res$rg[res$measure == "L2"])))
cat(sprintf("\nmax |rg_corrected + rg_uncorrected| = %.4f  (0 => exact sign flip)\n",
            max(abs(res$sum_check))))

cat("\n--- comparison with Gu et al. 2025 S4 (WM tracts: mean -0.0733, 14/15 neg) ---\n")
cat(sprintf("ours corrected ICVF mean  = %+.4f (%d/%d negative)\n",
            mean(icvf$rg), sum(icvf$rg < 0), nrow(icvf)))
cat("Gu  AutoPtx WM tract mean = -0.0733 (14/15 negative)\n")
cat("Gu  global CORTICAL ICVF  = -0.2189 (se 0.0702, p 1.8e-3)\n")

## BH FDR
res$q <- p.adjust(res$p, method = "BH")
icvf <- res[res$measure == "ICVF", ]
cat(sprintf("\nminimum BH q over the 15 ICVF tracts = %.4f\n",
            min(p.adjust(icvf$p, method = "BH"))))

write.csv(res, file.path(OUT, "p14_signfix_results.csv"), row.names = FALSE)

if (!is.null(L$I)) {
  I <- L$I; colnames(I) <- rownames(I) <- tn
  cat("\nASD univariate intercept:", round(I[1,1], 4), "\n")
  cat("cross-trait intercepts vs ASD (expect ~0): range [",
      round(min(I["ASDreal", -1]), 4), ",", round(max(I["ASDreal", -1]), 4), "]\n")
}
h2 <- diag(L$S)
cat(sprintf("ASD liability h2 = %.4f (unchanged by orientation, as expected)\n", h2[1]))

sink()

## ---- re-run the common-factor model on the corrected matrix ---------------
sink(file.path(OUT, "p14_signfix_factor.log"), split = TRUE)
cat("=== common white-matter factor, corrected orientation ===\n")
idps <- grep("^ICVF_", tn, value = TRUE)
m <- tryCatch(usermodel(covstruc = L, estimation = "DWLS",
                        model = paste0("F =~ ", paste(idps, collapse = " + "),
                                       "\nASDreal ~ F\nF ~~ 1*F\n"),
                        CFIcalc = TRUE, std.lv = TRUE),
              error = function(e) e)
if (inherits(m, "error")) {
  cat("model failed:", conditionMessage(m), "\n")
} else {
  r <- m$results
  print(r[r$op == "~", ])
  cat("fit: CFI =", m$modelfit$CFI, " SRMR =", m$modelfit$SRMR, "\n")
  cat("\n2026-09-17 uncorrected value was STD_Genotype = +0.1136, p = 0.0287.\n")
  saveRDS(m, file.path(OUT, "p14_signfix_factor.rds"))
}
sink()
cat("\nP14_SIGNFIX_DONE\n")

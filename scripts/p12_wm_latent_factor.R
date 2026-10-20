## §140 / P5.2 — white-matter genetic LATENT FACTOR, and ASD's relation to it.
##
## The 15 ICVF IDPs are near-collinear genetically (effective dimensions Li&Ji 5.00 /
## Galwey 4.84 of 15). Treating them as 15 tests is wrong in both directions: Bonferroni
## over-corrects, and reporting the best single tract over-states. The right object is a
## common genetic factor.
##
## Two models:
##   M1  common factor over the 15 ICVF IDPs; ASD regressed on the factor.
##   M2  same, but ALLERGY regressed on the factor -- so the ASD loading has a comparison
##       trait analysed identically (allergy localises to microglia, not OL, so this is a
##       meaningful contrast rather than a formality).

suppressMessages(library(GenomicSEM))
setwd("./outputs/rg")
L <- readRDS("rg_asd_vs_allergy_ldsc.rds")
tn <- colnames(L$S)
if (is.null(tn)) tn <- L$trait.names
idps <- setdiff(tn, c("ASDreal", "ALLERGY"))
cat("IDPs in factor:", length(idps), "\n")

out <- list()

## ---- common factor over the IDPs only -------------------------------------
keep <- which(tn %in% idps)
Sx <- L$S[keep, keep]; Vx <- L$V_Stand
## GenomicSEM needs a matched S/V; easiest is to refit commonfactor on the full object
cf <- tryCatch(commonfactor(covstruc = L, estimation = "DWLS"), error = function(e) e)
if (inherits(cf, "error")) {
  cat("commonfactor failed:", conditionMessage(cf), "\n")
} else {
  out$commonfactor <- cf$results
  cat("\n=== common factor over ALL traits (incl ASD/ALLERGY) ===\n")
  print(cf$results)
  cat("\nmodel fit:\n"); print(cf$modelfit)
}

## ---- explicit SEM: F =~ IDPs ; ASD ~ F  ------------------------------------
mk <- function(outcome) {
  paste0("F =~ ", paste(idps, collapse = " + "), "\n",
         outcome, " ~ F\n",
         "F ~~ 1*F\n")
}
for (o in c("ASDreal", "ALLERGY")) {
  m <- tryCatch(usermodel(covstruc = L, estimation = "DWLS", model = mk(o),
                          CFIcalc = TRUE, std.lv = TRUE),
                error = function(e) e)
  if (inherits(m, "error")) {
    cat("\n[", o, "] model failed:", conditionMessage(m), "\n")
    next
  }
  cat("\n=== F =~ 15 ICVF IDPs ;", o, "~ F ===\n")
  r <- m$results
  print(r[r$lhs == o | r$op == "~", ])
  cat("fit: CFI =", m$modelfit$CFI, " SRMR =", m$modelfit$SRMR,
      " chisq =", m$modelfit$chisq, " df =", m$modelfit$df, "\n")
  out[[o]] <- list(results = m$results, fit = m$modelfit)
}

saveRDS(out, "wm_latent_factor.rds")

sink("WM_LATENT_FACTOR_REPORT.txt")
cat("=== §140 white-matter genetic latent factor ===\n")
cat("date:", format(Sys.time()), "\n")
cat("IDPs:", length(idps), "| effective dimensions previously Li&Ji 5.00 / Galwey 4.84\n\n")
for (o in c("ASDreal", "ALLERGY")) {
  if (is.null(out[[o]])) { cat(o, ": model did not converge\n\n"); next }
  r <- out[[o]]$results
  reg <- r[r$op == "~", ]
  cat("---", o, "regressed on the common white-matter factor ---\n")
  print(reg)
  f <- out[[o]]$fit
  cat("fit: CFI =", f$CFI, " SRMR =", f$SRMR, " chisq =", f$chisq, " df =", f$df, "\n\n")
}
cat("Interpretation: a significant ASD loading on the common factor means the association is\n")
cat("with GLOBAL white-matter genetic architecture, not with any specific tract. A null\n")
cat("factor loading alongside nominal single-tract effects would instead indicate the\n")
cat("single-tract results are noise.\n")
sink()
cat("\ndone\n")

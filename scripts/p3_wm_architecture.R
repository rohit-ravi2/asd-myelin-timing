## Section XI — white-matter genetic architecture, done properly.
## P5.1 full LDSC table   P5.2 effective dimensionality   P5.3 tract class/length
## P5.4 ASD vs ALLERGY on the same IDPs (does allergy explain the ASD pattern?)
##
## The allergy GWAS here is the file formerly mislabelled ASD_Grove2019.gz (E-001).
## It is a real, well-powered allergic-disease GWAS (N=242,569, 81 GWS loci: TSLP, IL33,
## IL1RL1, ORMDL3/GSDMB, SMAD3, TLR1, EMSY/LRRC32, MHC-II) and is used here deliberately
## as the comparison trait, not as ASD.

suppressMessages(library(GenomicSEM))
OLD <- "data/asd_myelin"
NEW <- "./outputs/rg"
setwd(NEW)
hm3 <- "data/ldsc_ref/eur_w_ld_chr/w_hm3.snplist"
LD  <- "data/ldsc_ref/eur_w_ld_chr/"

## ---- munge the allergy GWAS -------------------------------------------------
if (!file.exists("ALLERGY.sumstats.gz")) {
  cat("[prep] allergy sumstats\n")
  d <- data.table::fread("data/dvad-corpus/asd_myelin/ASD_Grove2019.gz")
  d <- d[grepl("^rs", RS_ID) & PVALUE > 0 & is.finite(BETA)]
  data.table::fwrite(data.frame(SNP = d$RS_ID, A1 = toupper(d$EFFECT_ALLELE),
                                A2 = toupper(d$OTHER_ALLELE), BETA = d$BETA,
                                SE = d$SE, P = d$PVALUE),
                     "munge_ALLERGY.txt", sep = "\t")
  munge(files = "munge_ALLERGY.txt", hm3 = hm3, trait.names = "ALLERGY",
        N = 242569, info.filter = 0.9, maf.filter = 0.01)
}

## ---- ICVF panel only (keeps the pair count tractable) -----------------------
map <- read.delim(file.path(OLD, "idp_map_corrected.tsv"), stringsAsFactors = FALSE)
idps <- map$label[map$measure == "ICVF"]
idps <- idps[file.exists(paste0(idps, ".sumstats.gz"))]
tn <- c("ASDreal", "ALLERGY", idps)
cat("[panel]", length(tn), "traits\n")

L <- ldsc(traits = paste0(tn, ".sumstats.gz"),
          sample.prev = c(18381/46350, NA, rep(NA, length(idps))),
          population.prev = c(0.012, NA, rep(NA, length(idps))),
          ld = LD, wld = LD, trait.names = tn, stand = TRUE)
saveRDS(L, "rg_asd_vs_allergy_ldsc.rds")

rg <- L$S_Stand; colnames(rg) <- rownames(rg) <- tn
k <- length(tn)
idx <- which(lower.tri(rg, diag = TRUE), arr.ind = TRUE)
sev <- sqrt(diag(L$V_Stand)); semat <- matrix(NA, k, k)
for (i in seq_len(nrow(idx))) {
  semat[idx[i,1], idx[i,2]] <- sev[i]; semat[idx[i,2], idx[i,1]] <- sev[i]
}

## tract class + approximate length rank (P5.3)
cls <- c(ICVF_Genu_CC="commissural", ICVF_Body_CC="commissural", ICVF_Splenium_CC="commissural",
         ICVF_AntLimbIC_R="projection", ICVF_AntLimbIC_L="projection",
         ICVF_AntCoronaRad_R="projection", ICVF_AntCoronaRad_L="projection",
         ICVF_PostCoronaRad_R="projection", ICVF_PostCoronaRad_L="projection",
         ICVF_ExtCapsule_R="short_assoc", ICVF_ExtCapsule_L="short_assoc",
         ICVF_CingCingGyrus_R="long_assoc", ICVF_CingCingGyrus_L="long_assoc",
         ICVF_SLF_R="long_assoc", ICVF_SLF_L="long_assoc")

res <- data.frame(idp = idps,
                  tract_class = cls[idps],
                  rg_ASD = rg["ASDreal", idps], se_ASD = semat[1, match(idps, tn)],
                  rg_ALLERGY = rg["ALLERGY", idps], se_ALL = semat[2, match(idps, tn)])
res$z_ASD <- res$rg_ASD / res$se_ASD
res$p_ASD <- 2 * pnorm(-abs(res$z_ASD))
res$fdr_ASD <- p.adjust(res$p_ASD, "BH")
res$z_ALL <- res$rg_ALLERGY / res$se_ALL
res$p_ALL <- 2 * pnorm(-abs(res$z_ALL))
res$fdr_ALL <- p.adjust(res$p_ALL, "BH")
write.csv(res, "WM_ARCHITECTURE_asd_vs_allergy.csv", row.names = FALSE)

## P5.2 effective dimensionality of the IDP panel (Li & Ji on the IDP rg matrix)
R <- rg[idps, idps]; R[!is.finite(R)] <- 0; diag(R) <- 1
ev <- eigen((R + t(R)) / 2, only.values = TRUE)$values
ev <- pmax(ev, 0)
Meff_LiJi <- sum(ifelse(ev >= 1, 1, 0) + (ev - floor(ev)))
Meff_Galwey <- (sum(sqrt(ev))^2) / sum(ev)

sink("WM_ARCHITECTURE_REPORT.txt")
cat("=== Section XI: white-matter genetic architecture ===\n")
cat("date:", format(Sys.time()), "\n\n")
cat(sprintf("ASD h2 (liab)     = %.4f (SE %.4f, Z %.2f)\n",
            L$S[1,1], sqrt(L$V[1,1]), L$S[1,1]/sqrt(L$V[1,1])))
cat(sprintf("ALLERGY h2 (obs)  = %.4f (SE %.4f, Z %.2f)\n",
            L$S[2,2], sqrt(L$V[2,2]), L$S[2,2]/sqrt(L$V[2,2])))
cat(sprintf("rg(ASD, ALLERGY)  = %+.4f (SE %.4f, p %.3g)\n\n",
            rg["ASDreal","ALLERGY"], semat[1,2],
            2*pnorm(-abs(rg["ASDreal","ALLERGY"]/semat[1,2]))))
cat(sprintf("P5.2 effective number of IDP dimensions: Li&Ji=%.2f  Galwey=%.2f (of %d IDPs)\n",
            Meff_LiJi, Meff_Galwey, length(idps)))
cat(sprintf("     -> effective alpha for 0.05 family-wise = %.4f\n\n", 0.05/Meff_LiJi))
print(res[order(res$p_ASD), c("idp","tract_class","rg_ASD","se_ASD","p_ASD","fdr_ASD",
                              "rg_ALLERGY","p_ALL")], row.names = FALSE, digits = 3)
cat("\nP5.3 by tract class (ASD):\n")
for (cc in unique(res$tract_class)) {
  s <- res[res$tract_class == cc, ]
  cat(sprintf("   %-12s n=%d  mean rg=%+.4f  (range %+.3f .. %+.3f)\n",
              cc, nrow(s), mean(s$rg_ASD), min(s$rg_ASD), max(s$rg_ASD)))
}
cat("\nP5.4 ASD vs ALLERGY on the same IDPs:\n")
cat(sprintf("   mean rg ASD     = %+.4f\n", mean(res$rg_ASD)))
cat(sprintf("   mean rg ALLERGY = %+.4f\n", mean(res$rg_ALLERGY)))
ct <- cor.test(res$rg_ASD, res$rg_ALLERGY)
cat(sprintf("   correlation of the two rg profiles across IDPs: r=%.3f p=%.3g\n",
            ct$estimate, ct$p.value))
cat("   (a high correlation would mean the ASD pattern is not ASD-specific)\n")
sink()
cat("\ndone\n")

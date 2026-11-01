script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
setwd(dirname(normalizePath(script)))
## ddmm_ad_replication.R -- locus-level replication tables for the AD analysis
##
## Table 1  DDMM against the loci GhostKnockoff reports. We use the 34
##          representative loci of He et al. (2022), Supplementary Table 2 at
##          FDR 0.05, as the GhostKnockoff result; the method is not re-run.
## Table 2  Support for the DDMM loci in the nine studies He et al. combine.
##
## Variants on a chromosome form one locus while consecutive positions are
## within 1.2 Mb, and each locus is represented by its strongest variant. Two
## loci agree when their representatives lie within 1 Mb.
##
## Run ddmm_vs_ghostknockoff_realdata.R first, then this script, from the folder
## that contains them:  Rscript ddmm_ad_replication.R

suppressMessages(library(data.table))

DDIR <- "data/realdata_section5"
OUT  <- "output"
gap    <- 1.2e6
window <- 1.0e6

## ---- He et al. (2022), Supplementary Table 2: 34 loci at FDR 0.05, hg38 ----
he <- data.table(
  chr = c(1, 1, 1, 2, 2, 4, 5, 6, 6, 6, 8, 10, 10, 11, 11, 11, 14, 14, 15, 15,
          16, 16, 16, 17, 17, 17, 17, 17, 17, 19, 19, 19, 20, 21),
  pos = c(155419060, 161185602, 178992361, 44026309, 127135234, 11026080,
          87002714, 32637301, 40783137, 47479305, 27598736, 11678621, 29966853,
          47440232, 60212842, 86089237, 52710264, 92469490, 58889786, 63277703,
          17478817, 31121341, 70258841, 388402, 49391824, 58320645, 60239372,
          63482562, 73384739, 1046077, 44908684, 51224706, 56443204, 26161943),
  gene = c("POU5F1P4", "ADAMTS4", "FAM20B", "LRPPRC", "BIN1", "CLNK", "MIR4280",
           "HLA-DQA1", "LOC101929555", "CD2AP", "CLU", "ECHDC3", "JCAD",
           "RAPSN", "MIR6503", "PICALM", "PSMC6", "SLC24A4", "SLTM", "APH1B",
           "XYLT1", "KAT8", "AARS1", "LOC105371430", "LOC102724596",
           "TSPOAP1-AS1", "SCARNA20", "ACE", "SDK2", "ABCA7", "APOE", "CD33",
           "CASS4", "APP"))

## Schwartzentruber aggregates Kunkle with UK Biobank, and ADGC and ADSP are the
## in-house sets, so none of the three is independent of Kunkle.
studies <- c(Jansen = "Jansenetal_2019sept", Kunkle = "Kunkle_etal_Stage1",
             Schwartz = "Schwartzentruber_etal", Bis = "Bis_etal_WES",
             ADGC = "Stanford2021_GWAS_AD_ADGC",
             ADSP_WES = "Stanford2021_GWAS_AD_ADSP_WES",
             ADSP_WGS = "Stanford2021_GWAS_AD_ADSP_WGS",
             LeGuan = "LeGuan_etal_WES", IGAP = "IGAP_survival")
independent <- c("Jansen", "Kunkle", "Bis", "LeGuan", "IGAP")

## ---- inputs ----
disc_file <- file.path(OUT, "ddmm_vs_ghostknockoff_realdata_discoveries.csv")
if (!file.exists(disc_file)) stop("Run ddmm_vs_ghostknockoff_realdata.R first.")
disc <- fread(disc_file)
cand <- fread(file.path(DDIR, "AD_Zscores_candidate.csv"))
setnames(cand, c("hg38.chr", "hg38.pos"), c("chr", "pos"))

## Distance from a position to the nearest row of ref on the same chromosome.
nearest <- function(chr_i, pos_i, ref) {
  same <- ref[chr == chr_i]
  if (nrow(same) == 0) Inf else min(abs(same$pos - pos_i))
}

## ---- collapse the DDMM discoveries into loci ----
ddmm <- disc[ddmm_disc == TRUE]
setorder(ddmm, chr, pos)
ddmm[, locus := cumsum(c(TRUE, diff(pos) > gap)), by = chr]
ddmm <- ddmm[, .SD[which.max(abs(Z))], by = .(chr, locus)]
n_loci <- nrow(ddmm)

## ---- Table 1 ----
he[, covered := mapply(nearest, chr, pos, MoreArgs = list(ref = cand)) <= window]
he[, replicated := mapply(nearest, chr, pos, MoreArgs = list(ref = ddmm)) <= window]
n_cov <- sum(he$covered)
n_rep <- sum(he$covered & he$replicated)

cat("Table 1. GhostKnockoff loci of He et al. (2022) vs DDMM\n")
cat(sprintf("  reported at FDR 0.05        : %d\n", nrow(he)))
cat(sprintf("  in the candidate set        : %d\n", n_cov))
cat(sprintf("  replicated by DDMM          : %d of %d (%.0f%%)\n",
            n_rep, n_cov, 100 * n_rep / n_cov))
cat("  covered but missed          : ", paste(he[covered & !replicated, gene], collapse = ", "), "\n", sep = "")
cat("  absent from the candidates  : ", paste(he[!covered, gene], collapse = ", "), "\n\n", sep = "")
fwrite(he, file.path(OUT, "ad_he_loci_replication.csv"))

## ---- Table 2 ----
## z-scores of the nine studies at each DDMM locus representative.
z_by_study <- cand[ddmm[, .(chr, pos)], on = .(chr, pos)][, ..studies]
z_by_study <- as.matrix(z_by_study)
colnames(z_by_study) <- names(studies)
z_by_study[is.na(z_by_study) | z_by_study == 0] <- NA

supported <- function(threshold, group = names(studies)) {
  hits <- abs(z_by_study[, group, drop = FALSE]) >= threshold
  rowSums(hits, na.rm = TRUE) > 0
}
report <- function(label, hits) {
  cat(sprintf("  %-38s: %2d of %d (%3.0f%%)\n", label, sum(hits), n_loci,
              100 * sum(hits) / n_loci))
}

cat(sprintf("Table 2. Support for the %d DDMM loci\n", n_loci))
report("GhostKnockoff (He, covered loci)",
       mapply(nearest, ddmm$chr, ddmm$pos,
              MoreArgs = list(ref = he[covered == TRUE])) <= window)
for (s in c("Kunkle", "Jansen", "Schwartz")) {
  report(sprintf("%s, |z| >= 1.96", s), supported(1.96, s))
}
report("Kunkle, |z| >= 3", supported(3, "Kunkle"))
report("an independent study, |z| >= 1.96", supported(1.96, independent))
report("an independent study, |z| >= 3", supported(3, independent))

fwrite(data.table(chr = ddmm$chr, pos = ddmm$pos, gene = ddmm$gene,
                  Z = round(ddmm$Z, 2),
                  n_studies = rowSums(abs(z_by_study) >= 1.96, na.rm = TRUE),
                  n_independent = rowSums(abs(z_by_study[, independent, drop = FALSE]) >= 1.96, na.rm = TRUE)),
       file.path(OUT, "ad_ddmm_locus_support.csv"))

cat(sprintf("\nWrote %s and %s\n",
            file.path(OUT, "ad_he_loci_replication.csv"),
            file.path(OUT, "ad_ddmm_locus_support.csv")))

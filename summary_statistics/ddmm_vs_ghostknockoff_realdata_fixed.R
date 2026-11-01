script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
setwd(dirname(normalizePath(script)))
## ============================================================================
##  ddmm_vs_ghostknockoff_realdata_fixed.R  --  AD real data, fixed DDMM
##
##  Same data, ridge, seed and GhostKnockoff draws as
##  ddmm_vs_ghostknockoff_realdata.R. On every chromosome the SAME (z, Sigma)
##  go to
##
##    DDMM original   ddmm_pvalues(Sigma, z, n = 1)             (manuscript)
##    DDMM fixed      ddmm_pvalues_fixed(..., n_bic = N)        (FIX 1)
##    GhostKnockoff   unchanged, M = 5, equi
##
##  Why FIX 1 matters here. With n = 1 the original nodewise penalty is
##  sqrt(2 log p), about 3.7, while every entry of the ridged LD matrix is below
##  0.9. Every nodewise fit is therefore exactly zero and the "debiased"
##  statistic is z_j - Sigma_{j,-j} beta_hat_{-j}: a tag in LD with a causal
##  variant that the lasso shrank keeps most of its marginal signal. The fixed
##  version chooses the nodewise lambda per variant by BIC with the sample size
##  N behind the z-scores. Any N of a GWAS gives nearly the same answer; the
##  script reports N = 500, 5,000 and 50,000 as a sensitivity check and uses
##  N_MAIN (default 50,000) for the headline numbers and figures.
##
##  sigma = 1 is kept (z-score convention of the manuscript), so FIX 3 does not
##  apply; FIX 2 concerns group tests and does not apply either.
##
##  Output (in output_fixed/):
##     ddmm_vs_ghostknockoff_realdata_fixed_pvalues.csv      every variant
##     ddmm_vs_ghostknockoff_realdata_fixed_discoveries.csv  union of discoveries
##     ddmm_vs_ghostknockoff_realdata_fixed_summary.csv      counts, q = 0.05 / 0.10
##     ddmm_vs_ghostknockoff_realdata_fixed.png              three-panel Manhattan
##  USAGE:  Rscript ddmm_vs_ghostknockoff_realdata_fixed.R [chr_subset]
##          (optional: RD_CKPT=<dir> makes the run resumable chromosome by
##           chromosome; RD_BUDGET=<seconds> exits with status 99 when reached)
## ============================================================================
suppressMessages({ library(data.table); library(Matrix) })
set.seed(1)

DDIR <- "data/realdata_section5"
OUT  <- "output_fixed"
dir.create(OUT, showWarnings = FALSE)
source("ddmm_sumstat_fixed.R")
source("ghostknockoff.R")

s0 <- 0.10; M <- 5; gk_method <- "equi"
N_SENS <- c(500, 5000, 50000)
N_MAIN <- as.numeric(Sys.getenv("N_MAIN", "50000"))
args <- commandArgs(trailingOnly = TRUE)
CHRS <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else 1:22

Z  <- fread(file.path(DDIR, "AD_Zscores_candidate.csv"))
GA <- fread(file.path(DDIR, "Candidate_variants_info_other.csv"))
setnames(Z, c("hg38.chr", "hg38.pos"), c("chr", "pos"))
gene_of <- function(ch, ps) { g <- GA$HH_ref[GA$chr == ch & GA$pos == ps]; if (length(g)) g[1] else NA_character_ }

CKPT <- Sys.getenv("RD_CKPT", ""); if (nzchar(CKPT)) dir.create(CKPT, recursive = TRUE, showWarnings = FALSE)
BUDGET <- as.numeric(Sys.getenv("RD_BUDGET", "Inf")); T0 <- proc.time()[["elapsed"]]

meta <- list()
for (ch in CHRS) {
  cf <- file.path(DDIR, sprintf("Cor_%d.csv", ch)); if (!file.exists(cf)) stop("Missing ", cf, "; run download_data.R first.")
  ck <- file.path(CKPT, sprintf("chr%02d.rds", ch))
  if (nzchar(CKPT) && file.exists(ck)) {
    zz <- readRDS(ck); assign(".Random.seed", zz$seed, envir = globalenv())
    meta[[length(meta) + 1]] <- zz$dt; next
  }
  if (proc.time()[["elapsed"]] - T0 > BUDGET) { cat("budget reached; resume later\n"); quit(status = 99) }

  S <- fread(cf); you_index <- S[[1]]; Sigma <- as.matrix(S[, -1, with = FALSE])
  rows_ch <- which(Z$chr == ch)
  Zc <- Z[rows_ch][you_index]
  z  <- Zc$Z; ok <- is.finite(z); z[!ok] <- 0
  Sigma <- (1 - s0) * Sigma + s0 * diag(nrow(Sigma))
  Sigma <- as.matrix(forceSymmetric(Sigma))

  ## same order of calls as the original script: DDMM (no random numbers), then GK
  pv_orig <- ddmm_pvalues(Sigma, z, n = 1, q = 0)
  pv_fix  <- lapply(N_SENS, function(N) ddmm_pvalues_fixed(Sigma, z, n = 1, q = 0, n_bic = N))
  if (!N_MAIN %in% N_SENS) pv_fix[[length(pv_fix) + 1]] <- ddmm_pvalues_fixed(Sigma, z, n = 1, q = 0, n_bic = N_MAIN)
  names(pv_fix) <- paste0("p_fix_n", c(N_SENS, if (!N_MAIN %in% N_SENS) N_MAIN))

  prelim <- ghostknockoff_prelim(Sigma, M, method = gk_method)
  Zk <- ghostknockoff_zscores(z, prelim)
  gk <- ghostknockoff_filter(z^2, Zk^2)

  dt <- data.table(chr = ch, pos = Zc$pos, Z = z, p_orig = pv_orig,
                   T_fix = attr(pv_fix[[paste0("p_fix_n", N_MAIN)]], "T"),
                   kappa = gk$kappa, tau = gk$tau)
  for (nm in names(pv_fix)) dt[[nm]] <- as.numeric(pv_fix[[nm]])
  meta[[length(meta) + 1]] <- dt
  if (nzchar(CKPT)) saveRDS(list(dt = dt, seed = get(".Random.seed", envir = globalenv())), ck)
  cat(sprintf("  chr%-2d p=%4d | tested: orig %d, fixed %d | min p: orig %.1e, fixed %.1e\n",
              ch, length(z), sum(pv_orig < 1), sum(dt[[paste0("p_fix_n", N_MAIN)]] < 1),
              min(pv_orig), min(dt[[paste0("p_fix_n", N_MAIN)]])))
}
MM <- rbindlist(meta)
setorder(MM, chr, pos)
MM[, gene := mapply(gene_of, chr, pos)]
MM[, p_fix := get(paste0("p_fix_n", N_MAIN))]
P <- nrow(MM)

bh_sel <- function(pv, lvl) {
  o <- order(pv); bel <- pv[o] <= lvl * seq_along(pv) / length(pv)
  sel <- logical(length(pv)); if (any(bel)) sel[o[seq_len(max(which(bel)))]] <- TRUE; sel
}
ko_sel <- function(kappa, tau, lvl, M) {
  tss <- sort(unique(tau[kappa == 0 & tau > 0]))
  good <- vapply(tss, function(t) (1/M)*(1 + sum(kappa>=1 & tau>=t)) / max(1, sum(kappa==0 & tau>=t)) <= lvl, logical(1))
  if (!any(good)) return(list(sel = logical(length(tau)), thr = Inf))
  thr <- min(tss[good]); list(sel = (kappa == 0 & tau >= thr), thr = thr)
}
loci_lab <- function(k) if (k == 1) "1 locus" else paste(k, "loci")
## a locus = discoveries within 1 Mb of each other on the same chromosome
n_loci <- function(sel) {
  if (!any(sel)) return(0L)
  d <- MM[sel, .(chr, pos)][order(chr, pos)]
  sum(c(TRUE, diff(d$pos) > 1e6 | diff(d$chr) != 0))
}

summ <- list()
for (q in c(0.05, 0.10)) {
  gks <- ko_sel(MM$kappa, MM$tau, q, M)$sel
  MM[[sprintf("gk_q%02d", 100 * q)]] <- gks
  MM[[sprintf("orig_q%02d", 100 * q)]] <- bh_sel(MM$p_orig, q)
  for (N in c(N_SENS, if (!N_MAIN %in% N_SENS) N_MAIN))
    MM[[sprintf("fix_n%d_q%02d", N, 100 * q)]] <- bh_sel(MM[[paste0("p_fix_n", N)]], q)
  for (col in grep(sprintf("_q%02d$", 100 * q), names(MM), value = TRUE)) {
    sel <- MM[[col]]
    summ[[length(summ) + 1]] <- data.table(q = q, method = sub(sprintf("_q%02d$", 100 * q), "", col),
      variants = sum(sel), loci = n_loci(sel), shared_with_GK = sum(sel & gks))
  }
}
SUM <- unique(rbindlist(summ))
print(SUM)
fwrite(SUM, file.path(OUT, "ddmm_vs_ghostknockoff_realdata_fixed_summary.csv"))
fwrite(MM[, !c("T_fix")], file.path(OUT, "ddmm_vs_ghostknockoff_realdata_fixed_pvalues.csv"))
keep <- MM$orig_q05 | MM$orig_q10 | MM$gk_q05 | MM$gk_q10 |
        MM[[sprintf("fix_n%d_q05", N_MAIN)]] | MM[[sprintf("fix_n%d_q10", N_MAIN)]]
fwrite(MM[keep, .(chr, pos, gene, Z, p_orig, p_fix, tau, kappa,
                  orig_q05, fix_q05 = get(sprintf("fix_n%d_q05", N_MAIN)), gk_q05,
                  orig_q10, fix_q10 = get(sprintf("fix_n%d_q10", N_MAIN)), gk_q10)],
       file.path(OUT, "ddmm_vs_ghostknockoff_realdata_fixed_discoveries.csv"))

## ---- three-panel Manhattan at q = 0.05: original DDMM, fixed DDMM, GK --------
MM <- copy(MM); MM[, x := seq_len(.N)]
bnd <- MM[, .(mid = mean(x)), by = chr]
RED <- "#c9252d"; GREY <- ifelse(MM$chr %% 2 == 0, "#9aa7b5", "#d1d8df")
g1  <- vapply(strsplit(ifelse(is.na(MM$gene), "", MM$gene), ","), function(z) z[1], "")
lab <- function(sel, y, ytext) {
  keepg <- sel & g1 != ""; if (!any(keepg)) return(invisible())
  ax <- tapply(MM$x[keepg], g1[keepg], mean); ay <- tapply(y[keepg], g1[keepg], max)
  o <- order(ax); ax <- ax[o]; ay <- ay[o]
  lx <- ax; gap <- max(MM$x) * 0.013
  for (i in seq_along(lx)[-1]) if (lx[i] - lx[i - 1] < gap) lx[i] <- lx[i - 1] + gap
  segments(ax, ay, lx, ytext, col = "grey70", lwd = 0.4, xpd = NA)
  text(lx, ytext, names(ax), srt = 90, adj = c(0, 0.5), cex = 0.5, font = 3, col = RED, xpd = NA)
}
png(file.path(OUT, "ddmm_vs_ghostknockoff_realdata_fixed.png"), width = 1900, height = 1800, res = 160)
par(mfrow = c(3, 1), mar = c(1.2, 5.2, 2.6, 1.0), oma = c(2.5, 0, 0, 0))
panel <- function(yv, sel, ylab, title, cap) {
  yv <- pmin(yv, cap)
  plot(MM$x, yv, pch = 19, cex = 0.3, col = GREY, ylim = c(0, cap * 1.35), xaxt = "n",
       xlab = "", ylab = ylab, bty = "l", main = title, adj = 0, font.main = 2)
  points(MM$x[sel], yv[sel], pch = 19, cex = 0.7, col = RED)
  lab(sel, yv, cap * 1.02)
}
s_o <- MM$orig_q05; s_f <- MM[[sprintf("fix_n%d_q05", N_MAIN)]]; s_g <- MM$gk_q05
panel(-log10(pmax(MM$p_orig, 1e-300)), s_o, expression(-log[10] ~ p[j]),
      sprintf("A   DDMM, original nodewise (%d variants, %s)", sum(s_o), loci_lab(n_loci(s_o))), 30)
panel(-log10(pmax(MM$p_fix, 1e-300)), s_f, expression(-log[10] ~ p[j]),
      sprintf("B   DDMM, per-variant BIC nodewise, N = %s (%d variants, %s)",
              format(N_MAIN, big.mark = ","), sum(s_f), loci_lab(n_loci(s_f))), 30)
panel(MM$tau, s_g, expression(tau[j]),
      sprintf("C   GhostKnockoff, M = %d (%d variants, %s)", M, sum(s_g), loci_lab(n_loci(s_g))), 40)
axis(1, at = bnd$mid, labels = bnd$chr, tick = FALSE, line = -0.5, cex.axis = 0.78)
mtext("Chromosome", side = 1, line = 1, outer = TRUE)
dev.off()
cat("figure ->", file.path(OUT, "ddmm_vs_ghostknockoff_realdata_fixed.png"), "\nDone.\n")

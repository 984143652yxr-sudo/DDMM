script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
setwd(dirname(normalizePath(script)))
## ============================================================================
##  ddmm_vs_ghostknockoff_realdata.R   --   Section 5 real-data comparison
##
##  Alzheimer's disease meta-analysis summary statistics (GuJQ5/FWER_Knockoff_filter,
##  Section 5):  7,963 candidate variants on chr1-22, Z-scores in
##  AD_Zscores_candidate.csv, per-chromosome LD matrices Cor_1..Cor_22.csv.
##
##  Both methods consume ONLY (z, Sigma) per chromosome block, then FDR is
##  controlled GENOME-WIDE at alpha.
##
##  DDMM (summary-statistic form) is implemented in ddmm_sumstat.R and called
##  here on the Z-score scale, i.e. the n = 1, sigma = 1 case of Cov(z) = Sigma.
##  Per block: Gram-form lasso, nodewise debiasing of the active set, then
##  asymptotically standard normal T_j. p-values are pooled and BH is applied
##  genome-wide.
##
##  LD regularisation: Cor blocks are near-singular (lambda_min ~ 1e-3). We apply
##  the SAME ridge shrinkage  Sigma <- (1-s0)Sigma + s0 I  (s0=0.10) to BOTH
##  methods before any inversion -- a standard LD-shrinkage step, applied equally.
##
##  Output (NEW files):
##     ddmm_vs_ghostknockoff_realdata.png
##     ddmm_vs_ghostknockoff_realdata_discoveries.csv
##     ddmm_vs_ghostknockoff_realdata_report.txt
##  USAGE:  Rscript ddmm_vs_ghostknockoff_realdata.R [chr_subset e.g. 21,22]
## ============================================================================
suppressMessages({ library(data.table); library(Matrix) })
set.seed(1)

## Paths are relative: run this script from the folder that contains it.
DDIR <- "data/realdata_section5"
OUT  <- "output"
dir.create(OUT, showWarnings = FALSE)
source("ddmm_sumstat.R")

alpha <- 0.10         # genome-wide FDR
s0    <- 0.10         # LD ridge shrinkage, applied to both methods
M     <- 5            # knockoffs per variant; a single knockoff is far too
                      # stringent under dense LD
gk_method <- "equi"   # knockoff contrast: "equi" or "sdp" (sdp needs Rdsdp)

args <- commandArgs(trailingOnly = TRUE)
CHRS <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else 1:22

## ---- load Z-scores + gene annotation ----
Z  <- fread(file.path(DDIR, "AD_Zscores_candidate.csv"))
GA <- fread(file.path(DDIR, "Candidate_variants_info_other.csv"))   # chr,pos,HH_ref(genes)
setnames(Z, c("hg38.chr", "hg38.pos"), c("chr", "pos"))
gene_of <- function(ch, ps) { g <- GA$HH_ref[GA$chr == ch & GA$pos == ps]; if (length(g)) g[1] else NA_character_ }

## ---- GhostKnockoff, as in GhostKnockoff_hello.R ----
source("ghostknockoff.R")

## =====================================================================
##  loop over chromosomes: align z with LD block, run both methods
## =====================================================================
all_pv <- c(); all_kappa <- c(); all_tau <- c(); meta <- list()
for (ch in CHRS) {
  cf <- file.path(DDIR, sprintf("Cor_%d.csv", ch)); if (!file.exists(cf)) stop("Missing ", cf, "; run download_data.R first.")
  S <- fread(cf); you_index <- S[[1]]; Sigma <- as.matrix(S[, -1, with = FALSE])
  rows_ch <- which(Z$chr == ch)                       # avoid data.table NSE on i
  Zc <- Z[rows_ch][you_index]                          # chromosome z-scores, ordered as in Cor
  z  <- Zc$Z; ok <- is.finite(z); z[!ok] <- 0
  Sigma <- (1 - s0) * Sigma + s0 * diag(nrow(Sigma))   # ridge shrinkage (both methods)
  Sigma <- as.matrix(forceSymmetric(Sigma))
  lam_min <- max(min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values), 1e-6)

  ## The candidate z-scores are already adjusted for structure upstream by the
  ## contributing studies, so no further spectral filtering is applied (q = 0).
  pv <- ddmm_pvalues(Sigma, z, n = 1, q = 0)

  prelim <- ghostknockoff_prelim(Sigma, M, method = gk_method)
  Zk <- ghostknockoff_zscores(z, prelim)
  gk <- ghostknockoff_filter(z^2, Zk^2)           # method = 'marginal'

  meta[[length(meta)+1]] <- data.table(chr = ch, pos = Zc$pos, Z = z, p_ddmm = pv,
                                        kappa = gk$kappa, tau = gk$tau)
  cat(sprintf("  chr%-2d: p=%d, lambda_min=%.4f, |Z|>5.45: %d, DDMM min-p=%.1e\n",
              ch, length(z), lam_min, sum(abs(z) > 5.45), min(pv)))
}
MM <- rbindlist(meta)

## =====================================================================
##  genome-wide FDR:  BH (DDMM)  vs  knockoff+ filter (GhostKnockoff)
## =====================================================================
P <- nrow(MM)
## DDMM: Benjamini-Hochberg at alpha over all variants
o <- order(MM$p_ddmm); bel <- MM$p_ddmm[o] <= alpha * seq_len(P) / P
ddmm_sel <- if (any(bel)) o[seq_len(max(which(bel)))] else integer(0)
## GhostKnockoff: knockoff+ FDR filter on pooled (kappa, tau)
kappa <- MM$kappa; tau <- MM$tau
ts <- sort(unique(tau[kappa == 0 & tau > 0]))
ok <- vapply(ts, function(t) {
  num <- (1/M) * (1 + sum(kappa >= 1 & tau >= t)); den <- max(1, sum(kappa == 0 & tau >= t)); num/den <= alpha
}, logical(1))
gk_sel <- if (any(ok)) which(kappa == 0 & tau >= min(ts[ok])) else integer(0)

MM[, ddmm_disc := FALSE][ddmm_sel, ddmm_disc := TRUE]
MM[, gk_disc := FALSE][gk_sel, gk_disc := TRUE]
MM[, gene := mapply(gene_of, chr, pos)]

both <- which(MM$ddmm_disc & MM$gk_disc); only_dd <- which(MM$ddmm_disc & !MM$gk_disc); only_gk <- which(MM$gk_disc & !MM$ddmm_disc)
cat(sprintf("\nGenome-wide (alpha=%.2f): DDMM=%d, GhostKnockoff=%d, shared=%d, DDMM-only=%d, GK-only=%d\n",
            alpha, length(ddmm_sel), length(gk_sel), length(both), length(only_dd), length(only_gk)))

disc <- MM[ddmm_disc | gk_disc][order(chr, pos)]
fwrite(disc, file.path(OUT, "ddmm_vs_ghostknockoff_realdata_discoveries.csv"))

## =====================================================================
##  figure: He-style two-panel Manhattan (Nat Commun 2022, He et al., Fig 3C)
##    TOP = DDMM          : -log10 p ; BH at FDR 0.05 (red) / 0.10 (blue)
##    BOT = GhostKnockoff : W = tau  ; knockoff+ at FDR 0.05 (red) / 0.10 (blue)
##    colored gene labels + variant-density colour strip beneath each panel
## =====================================================================
setorder(MM, chr, pos)
MM[, x := seq_len(.N)]
chr_alt <- MM$chr %% 2
grey    <- ifelse(chr_alt == 0, "#3f3f3f", "#9b9b9b")
xspan   <- max(MM$x)
bnd     <- MM[, .(mid = mean(x)), by = chr]
g1      <- vapply(strsplit(ifelse(is.na(MM$gene), "", MM$gene), ","), function(z) z[1], "")
RED <- "#cc0000"; BLUE <- "#1f4fd8"

## ---- selections at two FDR levels ------------------------------------------
bh_sel <- function(pv, lvl) {                       # Benjamini-Hochberg -> logical
  o <- order(pv); bel <- pv[o] <= lvl * seq_along(pv) / length(pv)
  sel <- logical(length(pv)); if (any(bel)) sel[o[seq_len(max(which(bel)))]] <- TRUE; sel
}
ko_sel <- function(kappa, tau, lvl, M) {            # knockoff+ -> list(sel, thr)
  tss <- sort(unique(tau[kappa == 0 & tau > 0]))
  good <- vapply(tss, function(t) (1/M)*(1 + sum(kappa>=1 & tau>=t)) / max(1, sum(kappa==0 & tau>=t)) <= lvl, logical(1))
  if (!any(good)) return(list(sel = logical(length(tau)), thr = Inf))
  thr <- min(tss[good]); list(sel = (kappa == 0 & tau >= thr), thr = thr)
}
dd05 <- bh_sel(MM$p_ddmm, 0.05); dd10 <- bh_sel(MM$p_ddmm, 0.10)
g05  <- ko_sel(MM$kappa, MM$tau, 0.05, M); g10 <- ko_sel(MM$kappa, MM$tau, 0.10, M)
bh05 <- if (any(dd05)) max(MM$p_ddmm[dd05]) else NA; bh10 <- if (any(dd10)) max(MM$p_ddmm[dd10]) else NA

## ---- variant density per 1 Mb (for the colour strip) -----------------------
dens <- vapply(seq_len(P), function(i) sum(MM$chr == MM$chr[i] & abs(MM$pos - MM$pos[i]) <= 5e5), 0L)

## ---- 1-D label de-overlap + colored one-label-per-gene ---------------------
spread <- function(xs, gap, lo, hi) {
  pos <- as.numeric(xs); n <- length(pos); if (n <= 1) return(pmin(pmax(pos, lo), hi))
  for (it in 1:500) {
    moved <- FALSE
    for (i in 2:n) if (pos[i] - pos[i-1] < gap) { mid <- (pos[i]+pos[i-1])/2; pos[i-1] <- mid-gap/2; pos[i] <- mid+gap/2; moved <- TRUE }
    pos[1] <- max(pos[1], lo)
    for (i in 2:n) if (pos[i] - pos[i-1] < gap) pos[i] <- pos[i-1] + gap
    if (pos[n] > hi) { pos[n] <- hi; for (i in (n-1):1) if (pos[i+1]-pos[i] < gap) pos[i] <- pos[i+1]-gap }
    if (!moved) break
  }
  pos
}
label_genes <- function(red_mask, blue_mask, yv, ytext) {     # red overrides blue per gene
  keep <- (red_mask | blue_mask) & !is.na(g1) & g1 != ""
  if (!any(keep)) return(invisible())
  gg <- g1[keep]; xx <- MM$x[keep]; yy <- yv[keep]; isr <- red_mask[keep]
  ax <- tapply(xx, gg, mean); ay <- tapply(yy, gg, max); ar <- tapply(isr, gg, any)
  o <- order(ax); ax <- ax[o]; ay <- ay[o]; ar <- ar[o]; gn <- names(ax)
  lx <- spread(ax, gap = xspan * 0.012, lo = 0, hi = xspan)
  segments(ax, ay, lx, ytext, col = "grey65", lwd = 0.4, xpd = NA)
  text(lx, ytext, gn, srt = 90, adj = c(0, 0.5), cex = 0.55, font = 3,
       col = ifelse(ar, RED, BLUE), xpd = NA)
}
he_pal <- colorRampPalette(c("#e9e9e9","#1a7a1a","#7fbf00","#ffd400","#ff8c00","#e00000"))(100)
strip <- function() {
  br <- pmin(dens / max(quantile(dens, 0.99), 1), 1); cols <- he_pal[pmax(1, ceiling(br * 100))]
  par(mar = c(0.2, 5.0, 0.2, 1.2)); plot.new(); plot.window(xlim = c(0, xspan), ylim = c(0, 1))
  rect(MM$x - 0.5, 0, MM$x + 0.5, 1, col = cols, border = NA)
}

png(file.path(OUT, "ddmm_vs_ghostknockoff_realdata.png"), width = 2000, height = 1180, res = 150)
layout(matrix(1:4, 4, 1), heights = c(10, 0.7, 10, 0.7))

## ---- TOP: DDMM (-log10 p) ----
YCAP <- 30; logp <- pmin(-log10(pmax(MM$p_ddmm, 1e-300)), YCAP)
par(mar = c(0.5, 5.0, 2.6, 1.2))
plot(MM$x, logp, pch = 19, cex = 0.35, col = grey, ylim = c(-1, YCAP * 1.5), yaxt = "n",
     xaxt = "n", xlab = "", ylab = expression(-log[10] ~ p ~ ~ "(DDMM)"), bty = "l",
     main = "DDMM  (debiased lasso + Benjamini-Hochberg;  -log10 p truncated at 30)",
     cex.main = 1.0, font.main = 1, adj = 0)
axis(2, at = seq(0, YCAP, 10))
points(MM$x[dd10], logp[dd10], pch = 19, cex = 0.7, col = BLUE)
points(MM$x[dd05], logp[dd05], pch = 19, cex = 0.7, col = RED)
if (!is.na(bh10)) abline(h = -log10(bh10), lty = 2, col = BLUE, lwd = 1.0)
if (!is.na(bh05)) abline(h = -log10(bh05), lty = 2, col = RED,  lwd = 1.0)
label_genes(dd05, dd10, logp, ytext = YCAP * 0.86)
legend("topleft", c("FDR 0.05", "FDR 0.10"), col = c(RED, BLUE), pch = 19, bty = "n", cex = 0.85)
strip()

## ---- BOTTOM: GhostKnockoff (W = tau) ----
WCAP <- 100; W <- pmin(MM$tau, WCAP)
par(mar = c(2.4, 5.0, 2.6, 1.2))
plot(MM$x, W, pch = 19, cex = 0.35, col = grey, ylim = c(-3, WCAP * 1.55), yaxt = "n",
     xaxt = "n", xlab = "Chr", ylab = expression(W == tau[j] ~ ~ "(GhostKnockoff)"), bty = "l",
     main = sprintf("GhostKnockoff  (M=%d knockoffs, knockoff+ filter;  W truncated at 100)", M),
     cex.main = 1.0, font.main = 1, adj = 0)
axis(2, at = seq(0, WCAP, 20)); axis(1, at = bnd$mid, labels = bnd$chr, tick = FALSE, line = -0.4, cex.axis = 0.75)
points(MM$x[g10$sel], W[g10$sel], pch = 19, cex = 0.7, col = BLUE)
points(MM$x[g05$sel], W[g05$sel], pch = 19, cex = 0.7, col = RED)
if (is.finite(g10$thr)) abline(h = min(g10$thr, WCAP), lty = 2, col = BLUE, lwd = 1.0)
if (is.finite(g05$thr)) abline(h = min(g05$thr, WCAP), lty = 2, col = RED,  lwd = 1.0)
label_genes(g05$sel, g10$sel, W, ytext = WCAP * 0.84)
strip()
dev.off()
cat("figure ->", file.path(OUT, "ddmm_vs_ghostknockoff_realdata.png"), "\n")

## =====================================================================
##  publication figure and results at genome-wide FDR alpha = 0.05
## =====================================================================
q_pub <- 0.05
dd_pub <- dd05
gk_pub <- g05$sel
shared_pub <- dd_pub & gk_pub
pub_ddmm_genes <- c(
  "FCER1G", "CR1", "BIN1", "PILRA", "EPHA1", "TMEM106B",
  "TSPAN14", "PLCG2", "ABCA1", "ABCA7", "APOE"
)
pub_ddmm_labels <- dd_pub & g1 %in% pub_ddmm_genes
pub_gk_labels <- gk_pub & g1 %in% c("CR1", "BIN1", "APOE")

MM[, ddmm_disc_q05 := dd_pub]
MM[, gk_disc_q05 := gk_pub]
disc_pub <- MM[ddmm_disc_q05 | gk_disc_q05][order(chr, pos)]
fwrite(
  disc_pub,
  file.path(OUT, "ddmm_vs_ghostknockoff_realdata_q05_discoveries.csv")
)

PUB_RED <- "#c9252d"
PUB_GREY <- ifelse(MM$chr %% 2 == 0, "#9aa7b5", "#d1d8df")
pub_file <- file.path(OUT, "ddmm_vs_ghostknockoff_realdata_twopanel_q05.png")
png(pub_file, width = 1900, height = 1260, res = 160)
layout(matrix(1:2, 2, 1), heights = c(1, 1))

## ---- panel A: DDMM ---------------------------------------------------
pub_logp <- pmin(-log10(pmax(MM$p_ddmm, 1e-300)), 30)
par(mar = c(1.0, 5.2, 2.8, 1.0))
plot(
  MM$x, pub_logp,
  pch = 19, cex = 0.32, col = PUB_GREY,
  ylim = c(0, 38), xaxt = "n",
  xlab = "", ylab = expression(-log[10] ~ p[j]),
  bty = "l", main = "A   DDMM", adj = 0, font.main = 2
)
points(MM$x[dd_pub], pub_logp[dd_pub], pch = 19, cex = 0.72, col = PUB_RED)
if (is.finite(bh05)) {
  abline(h = -log10(bh05), lty = 2, lwd = 1.2, col = PUB_RED)
}
label_genes(pub_ddmm_labels, rep(FALSE, P), pub_logp, ytext = 27.5)

## ---- panel B: GhostKnockoff -----------------------------------------
pub_wcap <- 40
pub_w <- pmin(MM$tau, pub_wcap)
par(mar = c(3.8, 5.2, 2.8, 1.0))
plot(
  MM$x, pub_w,
  pch = 19, cex = 0.32, col = PUB_GREY,
  ylim = c(0, 44), xaxt = "n",
  xlab = "", ylab = expression(tau[j]),
  bty = "l", main = "B   GhostKnockoff", adj = 0, font.main = 2
)
points(MM$x[gk_pub], pub_w[gk_pub], pch = 19, cex = 0.72, col = PUB_RED)
if (is.finite(g05$thr)) {
  abline(h = min(g05$thr, pub_wcap), lty = 2, lwd = 1.2, col = PUB_RED)
}
axis(
  1, at = bnd$mid, labels = bnd$chr,
  tick = FALSE, line = -0.5, cex.axis = 0.78
)
mtext("Chromosome", side = 1, line = 2.6)
label_genes(pub_gk_labels, rep(FALSE, P), pub_w, ytext = 31.0)

dev.off()
cat("publication figure ->", pub_file, "\n")

pub_report_file <- file.path(
  DIR, "ddmm_vs_ghostknockoff_realdata_q05_report.txt"
)
sink(pub_report_file)
cat("DDMM and GhostKnockoff on Alzheimer's disease summary statistics\n\n")
cat(sprintf("Aligned variants: %d on chromosomes %s\n", P, paste(range(CHRS), collapse = "-")))
cat(sprintf("Genome-wide FDR target: alpha = %.2f\n", q_pub))
cat(sprintf("GhostKnockoff simultaneous knockoffs: M = %d\n\n", M))
cat(sprintf("DDMM discoveries: %d variants\n", sum(dd_pub)))
cat(sprintf("GhostKnockoff discoveries: %d variants\n", sum(gk_pub)))
cat(sprintf("Shared discoveries: %d variants\n", sum(shared_pub)))
cat(sprintf("DDMM-only discoveries: %d variants\n", sum(dd_pub & !gk_pub)))
cat(sprintf("GhostKnockoff-only discoveries: %d variant\n", sum(!dd_pub & gk_pub)))
sink()
cat("publication report ->", pub_report_file, "\n")

## =====================================================================
##  report
## =====================================================================
sink(file.path(OUT, "ddmm_vs_ghostknockoff_realdata_report.txt"))
cat(strrep("=", 72), "\n")
cat("  DDMM (summary-statistic) vs GhostKnockoff -- AD real data (Section 5)\n")
cat(strrep("=", 72), "\n\n")
cat(sprintf("Variants: %d on chr %s.  Per-chromosome LD blocks (Cor_*.csv),\n", P, paste(range(CHRS), collapse = "-")))
cat(sprintf("ridge shrinkage s0=%.2f applied to BOTH methods. Genome-wide FDR alpha=%.2f.\n", s0, alpha))
cat(sprintf("GhostKnockoff: M=%d knockoffs, multiple-knockoff FDR filter (Gimenez-Zou\n", M))
cat("2019; He et al. 2022) -- the FWER filter of the repo is replaced by FDR.\n\n")
cat(sprintf("Discoveries:  DDMM=%d   GhostKnockoff=%d   shared=%d\n", length(ddmm_sel), length(gk_sel), length(both)))
cat(sprintf("              DDMM-only=%d   GhostKnockoff-only=%d\n\n", length(only_dd), length(only_gk)))
cat("Discovered variants (chr:pos, Z, gene, which method):\n")
pr <- disc[, .(chr, pos, Z = round(Z, 2), p_ddmm = signif(p_ddmm, 2),
               method = ifelse(ddmm_disc & gk_disc, "both", ifelse(ddmm_disc, "DDMM", "GK")), gene)]
print(pr, nrow = nrow(pr), row.names = FALSE)
cat("\n", strrep("-", 72), "\n READING\n", strrep("-", 72), "\n")
cat("* Both methods consume only (z, Sigma) and control FDR genome-wide; the\n")
cat("  discovery lists are directly comparable. Shared hits are the high-\n")
cat("  confidence AD loci both methods agree on.\n")
cat("* DDMM (debiased lasso + BH) is deterministic; GhostKnockoff resamples a\n")
cat("  knockoff (seed-fixed here). Differences concentrate at LD-dense regions.\n")
sink()
cat("report ->", file.path(OUT, "ddmm_vs_ghostknockoff_realdata_report.txt"), "\nDone.\n")

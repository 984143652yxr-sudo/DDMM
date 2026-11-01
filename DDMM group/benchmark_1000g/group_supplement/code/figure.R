## Figure of the group benchmark: columns EUR, MIX, ALL; rows locus power at amplitude 8,
## locus power at amplitude 4, realised FDR at amplitude 8, realised FDR at amplitude 4.
## DDMM (effective rank, BH) and KnockoffGWAS (knockoff+) per resolution, from the coarsest
## (res6) to the finest (res0) adjclust partition closed under r^2 >= 0.8. SuSiE has no
## fixed resolution and is drawn as horizontal reference lines (95% and 90% credible sets).
## The grey band marks res6, where every group contains a causal variant.
##   BENCH=<folder> Rscript figure.R  -> results/benchmark_group.pdf / .png
suppressMessages(library(data.table))
B <- Sys.getenv("BENCH"); S <- fread(file.path(B, "results", "summary.csv")); DZ <- fread(file.path(B, "results", "design_summary.csv"))
COH <- c("EUR", "MIX", "ALL"); RES <- sprintf("res%d", 6:0); x <- seq_along(RES)
LAB <- c(EUR = "EUR (n = 503)", MIX = "MIX (n = 504)", ALL = "ALL (n = 2,504)")
COL <- c(DDMM = "#1f4e79", KnockoffGWAS = "#d95f02", CS95 = "#555555", CS90 = "#2a9d8f")
val <- function(coh, cs, meth, rl, col) { d <- S[cohort == coh & case == cs & method == meth & rule == rl]; d[match(RES, d$resolution)][[col]] }
ref <- function(coh, cs, s, col) S[cohort == coh & case == cs & method == paste0("SuSiE_", s) & rule == s][[col]]
xlab_w <- function(coh) { z <- DZ[cohort == coh][match(RES, resolution)]
  paste0(RES, "\n", ifelse(z$median_width_kb >= 1, sprintf("%.0f kb", z$median_width_kb), "LD cl.")) }
draw <- function(file, dev) {
  if (dev == "pdf") pdf(file, width = 10.5, height = 11) else png(file, width = 2100, height = 2200, res = 200)
  par(mfrow = c(4, 3), mar = c(4.3, 4.2, 2.2, 0.8), mgp = c(2.4, 0.7, 0), cex.axis = 0.8, oma = c(0, 0, 0, 0))
  li <- 0
  for (what in c("power", "FDR")) for (cs in c("a8", "a4")) for (coh in COH) {
    yl <- if (what == "power") c(0, 1.12) else c(0, 0.3)
    plot(NA, xlim = c(0.7, 7.3), ylim = yl, xaxt = "n", yaxt = "n", xlab = "", ylab = if (what == "power") "Locus power" else "Realised FDR",
         main = sprintf("(%s) %s, amplitude %s, %s", letters[li <- li + 1], if (what == "power") "Power" else "FDR", sub("a", "", cs), LAB[coh]), cex.main = 0.92)
    rect(0.5, -1, 1.5, 2, col = "grey93", border = NA)
    axis(1, at = x, labels = xlab_w(coh), padj = 0.5, cex.axis = 0.72); axis(2, at = if (what == "power") seq(0, 1, 0.2) else seq(0, 0.3, 0.1))
    if (what == "FDR") abline(h = 0.1, col = "grey40", lty = 3)
    for (s in c("CS95", "CS90")) { v <- ref(coh, cs, s, what); if (length(v)) abline(h = v, col = COL[s], lwd = 1.6, lty = if (s == "CS95") 2 else 4) }
    lines(x, val(coh, cs, "DDMM_rank_effective", "BH", what), col = COL["DDMM"], lwd = 2.2, type = "o", pch = 16, cex = 0.8)
    lines(x, val(coh, cs, "KnockoffGWAS", "knockoff+", what), col = COL["KnockoffGWAS"], lwd = 2.2, type = "o", pch = 17, cex = 0.8)
    if (li == 1) legend("bottomleft", c("DDMM, effective rank", "KnockoffGWAS", "SuSiE, 95% CS", "SuSiE, 90% CS"), col = COL,
                        lwd = c(2.2, 2.2, 1.6, 1.6), lty = c(1, 1, 2, 4), pch = c(16, 17, NA, NA), bty = "n", cex = 0.78, bg = "white")
  }
  dev.off()
}
draw(file.path(B, "results", "benchmark_group.pdf"), "pdf"); draw(file.path(B, "results", "benchmark_group.png"), "png")

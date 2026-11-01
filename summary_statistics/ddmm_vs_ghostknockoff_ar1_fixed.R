script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
setwd(dirname(normalizePath(script)))
## ============================================================================
##  ddmm_vs_ghostknockoff_ar1_fixed.R  --  AR(1) regime sweep, fixed DDMM
##
##  Same design, seed and random stream as ddmm_vs_ghostknockoff_ar1.R. In every
##  replicate the SAME z-scores are given to four procedures:
##
##    DDMM original      ddmm_pvalues(), as in the manuscript
##    DDMM FIX 1         per-variant BIC nodewise lambda, sigma = 1
##    DDMM FIX 1 + 3     per-variant BIC nodewise lambda, sigma estimated from
##                       (s_y, r, Sigma) with the degrees-of-freedom correction
##    GhostKnockoff      unchanged
##
##  The fixed variants draw no random numbers, so the original DDMM and
##  GhostKnockoff columns are exactly those of the original script at the same
##  number of replicates.
##
##  REPS = 100, the number stated in the manuscript (the original script
##  defaults to 60).
##
##  Output (in output_fixed/):
##     ddmm_vs_ghostknockoff_ar1_fixed_results.csv   (means and SEs per cell)
##     ddmm_vs_ghostknockoff_ar1_fixed_raw.csv       (one row per replicate)
##     ddmm_vs_ghostknockoff_ar1_fixed.png
##  USAGE:  Rscript ddmm_vs_ghostknockoff_ar1_fixed.R
##          (optional: AR1_CKPT=<dir> makes the run resumable cell by cell)
## ============================================================================
suppressMessages({ library(Matrix) })
set.seed(1)

OUTDIR <- "output_fixed"
dir.create(OUTDIR, showWarnings = FALSE)
source("ddmm_sumstat_fixed.R")        # also sources the unchanged ddmm_sumstat.R

## ---- parameters (as in the original, except REPS) ---------------------------
p       <- as.integer(Sys.getenv("AR1_P", "500"))
n       <- as.integer(Sys.getenv("AR1_N", "500"))
k       <- as.integer(Sys.getenv("AR1_K", "30"))
amp.seq <- as.numeric(strsplit(Sys.getenv("AR1_AMPLITUDES", "3,7,11,15,19"), ",")[[1]])
REPS    <- as.integer(Sys.getenv("AR1_REPS", "100"))
stopifnot(p >= 5, n > 1, k > 0, k < p, REPS >= 2, all(is.finite(amp.seq)))
M       <- 5
alpha   <- 0.10
s0      <- 0.10
gk_method <- "equi"
regimes <- list(list(name = "well-conditioned LD", rho = 0.50),
                list(name = "dense LD (real-data regime)", rho = 0.95))
cat(sprintf("AR(1) regime sweep (fixed): p=%d n=%d k=%d M=%d alpha=%.2f s0=%.2f reps=%d\n",
            p, n, k, M, alpha, s0, REPS))

source("ghostknockoff.R")

setup <- function(rho) {
  J       <- abs(outer(seq_len(p), seq_len(p), "-"))
  Sigma   <- (1 - s0) * rho^J + s0 * diag(p)
  ev      <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  list(Sigma = Sigma,
       Sig_chol = chol(Sigma),
       lam_min = min(ev),
       gk_s = min(1, (1 + 1/M) * min(ev)),
       prelim = ghostknockoff_prelim(Sigma, M, method = gk_method),
       cache = nodewise_cache(Sigma))   # nodewise fits for this Sigma (FIX 1)
}

## Sigma is fixed within a regime, so the per-variant nodewise fits are
## computed once for every variant and reused in all replicates.
nodewise_cache <- function(Sigma) {
  e <- new.env()
  nw <- nodewise_bic_gram(Sigma, seq_len(ncol(Sigma)), n_bic = n)
  for (j in seq_len(ncol(Sigma))) assign(as.character(j), nw[[j]], envir = e)
  e
}
NCORES <- as.integer(Sys.getenv("NCORES", "4"))

gk_select <- function(z, G) {
  Zk <- ghostknockoff_zscores(z, G$prelim)
  which(ghostknockoff_filter(z^2, Zk^2)$q <= alpha)
}

## z = X_s' y_s / sqrt(n) with standardised y, so r = z / sqrt(n) and s_y = 1
ddmm_select_orig <- function(z, G) bh_select(ddmm_pvalues(G$Sigma, z / sqrt(n), n = n), alpha)
ddmm_select_fix1 <- function(z, G)
  bh_select(ddmm_pvalues_fixed(G$Sigma, z / sqrt(n), n = n, sigma = "unit", cache = G$cache), alpha)
ddmm_select_fix13 <- function(z, G)
  bh_select(ddmm_pvalues_fixed(G$Sigma, z / sqrt(n), n = n, sigma = "estimate", s_y = 1,
                               cache = G$cache), alpha)

metrics <- function(sel, causal) {
  if (length(sel) == 0) return(c(fdp = 0, pow = 0))
  tp <- sum(sel %in% causal); c(fdp = (length(sel) - tp) / length(sel), pow = tp / length(causal))
}

## ---- optional cell-level checkpointing --------------------------------------
CKPT <- Sys.getenv("AR1_CKPT", "")
if (nzchar(CKPT)) dir.create(CKPT, recursive = TRUE, showWarnings = FALSE)
cell_file <- function(ir, amp) file.path(CKPT, sprintf("cell_r%d_a%02d.rds", ir, amp))
BUDGET <- as.numeric(Sys.getenv("AR1_BUDGET", "Inf")); T0 <- proc.time()[["elapsed"]]

## ---- RUN --------------------------------------------------------------------
raw <- list()
for (ir in seq_along(regimes)) {
  rg <- regimes[[ir]]
  gfile <- file.path(CKPT, sprintf("setup_r%d.rds", ir))
  if (nzchar(CKPT) && file.exists(gfile)) G <- readRDS(gfile) else {
    G <- setup(rg$rho)                      # deterministic, no random numbers
    if (nzchar(CKPT)) saveRDS(G, gfile)
  }
  cat(sprintf("\n[%s] rho=%.2f lambda_min=%.4f knockoff s=%.4f\n", rg$name, rg$rho, G$lam_min, G$gk_s))
  for (amp in amp.seq) {
    if (nzchar(CKPT) && file.exists(cell_file(ir, amp))) {
      z <- readRDS(cell_file(ir, amp))
      assign(".Random.seed", z$seed, envir = globalenv())
      raw[[length(raw) + 1]] <- z$raw
      next
    }
    if (proc.time()[["elapsed"]] - T0 > BUDGET) { cat("budget reached; resume later\n"); quit(status = 99) }
    ## The random stream is consumed only by the data and by GhostKnockoff,
    ## in the same order as the original script. DDMM draws nothing, so it is
    ## evaluated afterwards, in parallel over replicates, on the stored z.
    zs <- vector("list", REPS); cs <- vector("list", REPS); gk <- matrix(0, REPS, 2)
    for (rep in seq_len(REPS)) {
      causal <- sort(sample(p, k))
      beta <- numeric(p); beta[causal] <- sample(c(-1, 1), k, TRUE) * amp / sqrt(n)
      X <- matrix(rnorm(n * p), n, p) %*% G$Sig_chol
      Y <- as.numeric(X %*% beta) + rnorm(n)
      z <- as.numeric(crossprod(scale(X), scale(Y))) / sqrt(n)
      zs[[rep]] <- z; cs[[rep]] <- causal
      gk[rep, ] <- metrics(gk_select(z, G), causal)
    }
    dd <- parallel::mclapply(seq_len(REPS), function(rep)
      c(metrics(ddmm_select_orig(zs[[rep]], G), cs[[rep]]),
        metrics(ddmm_select_fix1(zs[[rep]], G), cs[[rep]]),
        metrics(ddmm_select_fix13(zs[[rep]], G), cs[[rep]])),
      mc.cores = NCORES, mc.set.seed = FALSE)
    out <- cbind(do.call(rbind, dd), gk)
    colnames(out) <- c("DDMM_orig_fdp", "DDMM_orig_pow", "DDMM_fix1_fdp", "DDMM_fix1_pow",
                       "DDMM_fix13_fdp", "DDMM_fix13_pow", "GK_fdp", "GK_pow")
    cell <- data.frame(regime = rg$name, rho = rg$rho, lam_min = G$lam_min, s = G$gk_s,
                       amplitude = amp, rep = seq_len(REPS), out)
    raw[[length(raw) + 1]] <- cell
    if (nzchar(CKPT)) saveRDS(list(raw = cell, seed = get(".Random.seed", envir = globalenv())),
                              cell_file(ir, amp))
    cat(sprintf("  amp=%2d | orig fdp=%.3f pow=%.3f | fix1 fdp=%.3f pow=%.3f | fix1+3 fdp=%.3f pow=%.3f | GK fdp=%.3f pow=%.3f\n",
                amp, mean(out[, 1]), mean(out[, 2]), mean(out[, 3]), mean(out[, 4]),
                mean(out[, 5]), mean(out[, 6]), mean(out[, 7]), mean(out[, 8])))
  }
}
RAW <- do.call(rbind, raw)
write.csv(RAW, file.path(OUTDIR, "ddmm_vs_ghostknockoff_ar1_fixed_raw.csv"), row.names = FALSE)

cols <- grep("_(fdp|pow)$", names(RAW), value = TRUE)
agg  <- aggregate(RAW[, cols], RAW[, c("regime", "rho", "lam_min", "s", "amplitude")], mean)
se   <- aggregate(RAW[, cols], RAW[, c("regime", "rho", "amplitude")], function(x) sd(x) / sqrt(length(x)))
names(se)[-(1:3)] <- paste0(names(se)[-(1:3)], "_se")
R <- merge(agg, se, by = c("regime", "rho", "amplitude"))
R <- R[order(R$rho, R$amplitude), ]
write.csv(R, file.path(OUTDIR, "ddmm_vs_ghostknockoff_ar1_fixed_results.csv"), row.names = FALSE)
print(format(R[, c("rho", "amplitude", cols)], digits = 3), row.names = FALSE)

## ---- FIGURE: rows = {power, FDP}, cols = regimes -----------------------------
arms <- list(list(k = "DDMM_orig",  lab = "DDMM, original",       col = "#8c8c8c", pch = 1,  lty = 2),
             list(k = "DDMM_fix1",  lab = "DDMM, FIX 1 (sigma=1)", col = "#6baed6", pch = 15, lty = 1),
             list(k = "DDMM_fix13", lab = "DDMM, FIX 1+3",        col = "#1f4fd8", pch = 19, lty = 1),
             list(k = "GK",         lab = "GhostKnockoff (M=5)",  col = "#cc0000", pch = 17, lty = 1))
png(file.path(OUTDIR, "ddmm_vs_ghostknockoff_ar1_fixed.png"), width = 1700, height = 1150, res = 150)
par(mfrow = c(2, 2), mar = c(4.4, 4.8, 3.4, 1.2), oma = c(0, 0, 2.2, 0))
for (what in c("pow", "fdp")) for (rg in regimes) {
  Rr <- R[R$rho == rg$rho, ]
  ymax <- if (what == "pow") 1 else max(0.4, max(R[, grep("_fdp$", names(R))]) * 1.1)
  plot(NA, xlim = range(amp.seq), ylim = c(0, ymax), xlab = "signal amplitude",
       ylab = if (what == "pow") "power" else "FDP (mean)",
       main = sprintf("%s, rho=%.2f", rg$name, rg$rho), cex.main = 0.95)
  if (what == "fdp") abline(h = alpha, lty = 3, col = "grey40")
  for (a in arms) {
    m <- Rr[[paste0(a$k, "_", what)]]; s <- Rr[[paste0(a$k, "_", what, "_se")]]
    ok <- is.finite(s) & s > 0
    if (any(ok)) arrows(Rr$amplitude[ok], m[ok] - s[ok], Rr$amplitude[ok], m[ok] + s[ok], angle = 90, code = 3, length = 0.02, col = a$col)
    lines(Rr$amplitude, m, type = "b", col = a$col, pch = a$pch, lty = a$lty)
  }
  if (what == "pow" && rg$rho == 0.5)
    legend("topleft", sapply(arms, `[[`, "lab"), col = sapply(arms, `[[`, "col"),
           pch = sapply(arms, `[[`, "pch"), lty = sapply(arms, `[[`, "lty"), bty = "n", cex = 0.8)
}
mtext(sprintf("Summary-statistic DDMM vs GhostKnockoff, AR(1), %d replicates, target FDR %.2f (error bars: 1 SE)",
              REPS, alpha), outer = TRUE, cex = 0.95, font = 2)
dev.off()
cat("figure ->", file.path(OUTDIR, "ddmm_vs_ghostknockoff_ar1_fixed.png"), "\nDone.\n")

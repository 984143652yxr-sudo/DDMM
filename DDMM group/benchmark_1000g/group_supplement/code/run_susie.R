## ============================================================
##  run_susie.R — SuSiE reference for the group benchmark of one cohort.
##  Sex and the intercept are partialled out of y and X, as in run_ddmm.R. One SuSiE fit
##  per phenotype (susieR::susie, L = 40, the number of causal variants, other arguments
##  at their defaults). From the same fit:
##    CS95  credible sets at coverage 0.95 (the default of Wang et al. 2020 and susieR),
##    CS90  credible sets at coverage 0.90 (susie_get_cs, no refit),
##    PIP90 single variants with posterior inclusion probability > 0.9 (the rule of the
##          semi-synthetic single-variant experiment).
##  Credible sets keep the default purity filter (smallest |r| within the set >= 0.5).
##  Writes susie_selections.rds (per phenotype: variant indices of each set, PIPs, fit
##  diagnostics).
##    DDMM_PKG not needed.  Rscript run_susie.R <cohort dir> [L = 40]
## ============================================================
suppressMessages({ library(data.table); library(susieR) })
args <- commandArgs(TRUE); dir <- args[1]; L <- if (length(args) > 1) as.integer(args[2]) else 40L
design <- readRDS(file.path(dir, "design.rds")); pheno <- fread(file.path(dir, "phenotypes.tab"))
out_file <- file.path(dir, "susie_selections.rds")
done <- if (file.exists(out_file)) readRDS(out_file) else list()
X <- design$X; n <- nrow(X); sex <- design$samples$sex; C <- cbind(1, sex)
resid_on <- function(M) M - C %*% solve(crossprod(C), crossprod(C, M))
Xs <- scale(resid_on(X * 1.0)); Xs[is.na(Xs)] <- 0
t0 <- proc.time()[[3]]
for (yc in setdiff(grep("^y_", names(pheno), value = TRUE), names(done))) {
  y <- as.numeric(resid_on(matrix(pheno[[yc]], ncol = 1)))
  fit <- susie(Xs, y, L = L)
  cs95 <- susie_get_cs(fit, X = Xs, coverage = 0.95)
  cs90 <- susie_get_cs(fit, X = Xs, coverage = 0.90)
  pip <- as.numeric(fit$pip)
  done[[yc]] <- list(CS95 = unname(cs95$cs), CS90 = unname(cs90$cs), PIP90 = which(pip > 0.9),
                     purity95 = cs95$purity, purity90 = cs90$purity, pip = round(pip, 5),
                     niter = fit$niter, converged = fit$converged, sigma2 = fit$sigma2)
  saveRDS(done, out_file)
  cat(sprintf("%s %s: %d CS95, %d CS90, %d PIP>0.9, %d iterations%s  [%.0f s]\n", format(Sys.time(), "%H:%M:%S"), yc,
              length(cs95$cs), length(cs90$cs), sum(pip > 0.9), fit$niter, if (fit$converged) "" else " (not converged)", proc.time()[[3]] - t0))
}

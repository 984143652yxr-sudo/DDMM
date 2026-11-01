# DDMM group tests (corrected version, code_official/ddmm_group/R) on one
# cohort of the 1000 Genomes benchmark, on the KnockoffGWAS partitions.
#
#   Rscript run_ddmm.R data/EUR
#
# Sex is partialled out of y and X (the authors' phenotypes contain a sex
# effect; KnockoffGWAS's lasso.R uses sex as a covariate). For every phenotype
# column of phenotypes.tab and every resolution: Simes and the chi-squared
# group test at the effective rank, K = 3 and full rank; selection by BH at
# q = 0.1 and by the truncated rule of the theory (u >= u_m). Writes
# ddmm_selections.rds (selected groups per phenotype x resolution x track x rule).

suppressMessages(library(data.table))
script_dir <- function() { a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(gsub('~+~', ' ', sub('^--file=', '', a[1]), fixed = TRUE))) else 'code' }
B <- Sys.getenv('BENCH', normalizePath(file.path(script_dir(), '..')))      # this folder
PKG <- Sys.getenv('DDMM_PKG', normalizePath(file.path(B, '..', '..')))     # ddmm_group package root (benchmark_1000g/group_supplement -> ddmm_group)
for (f in c('utils.R', 'ddmm.R')) source(file.path(PKG, 'R', f))
args <- commandArgs(TRUE); dir <- args[1]

# Pilot Lasso path. glmnet ends its path at 1e-4 x lambda_max when n > p and
# at 0.01 x lambda_max when n < p. With n = 2,504 > p = 2,000 (ALL cohort) and
# variants in near-perfect LD, the long path makes cv.glmnet take more than ten
# minutes per phenotype. We fix the ratio at 0.01, the value glmnet uses for the
# toy data (n = 1,000) and the EUR cohort (n = 503), so all runs use one path.
# Apart from this argument the function is R/ddmm.R::variant_tests.
variant_tests <- function(y, xcen, nodewise, lambda_min_ratio = 0.01) {
  lam <- glmnet::cv.glmnet(xcen, y, intercept = FALSE, standardize = FALSE, nfolds = 5,
                           lambda.min.ratio = lambda_min_ratio)$lambda.min
  beta <- as.numeric(glmnet::glmnet(xcen, y, lambda = lam, intercept = FALSE,
                                    standardize = FALSE)$beta[, 1])
  S <- which(beta != 0)
  r <- as.numeric(y - xcen %*% beta)
  if (length(S) < 3) S <- marginal_screen(y, xcen, S, min_selected = 3, top = 20)
  # sigma^2 = RSS / ||(I - P_S)||_F^2 = RSS / (n - rank X_S)  (FIX 3)
  sigma <- if (length(S)) {
    d <- svd(xcen[, S, drop = FALSE], nu = 0, nv = 0)$d
    sqrt(sum(r^2) / (nrow(xcen) - sum(d > 1e-10)))
  } else sqrt(sum(r^2) / nrow(xcen))
  statistic <- numeric(ncol(xcen)); p_value <- rep(1, ncol(xcen))
  if (length(S)) {
    ensure_nodewise(nodewise, S)
    v <- as.numeric(crossprod(xcen, r))
    zy <- zz <- zx <- numeric(length(S))
    for (k in seq_along(S)) {
      j <- S[k]; f <- nodewise$fits[[j]]
      zz[k] <- f$zz; zx[k] <- f$zx; zy[k] <- v[j] - sum(f$g * v[-j])
    }
    debiased <- beta[S] + zy / zx
    statistic[S] <- debiased * zx / (sigma * sqrt(zz))
    p_value[S] <- 2 * stats::pnorm(-abs(statistic[S]))
  }
  list(p = p_value, statistic = statistic, S = S, sigma = sigma, beta = beta)
}
design <- readRDS(file.path(dir, 'design.rds'))
pheno <- fread(file.path(dir, 'phenotypes.tab'))
out_file <- file.path(dir, 'ddmm_selections.rds')
done <- if (file.exists(out_file)) readRDS(out_file) else list()

X <- design$X; n <- nrow(X); p <- ncol(X)
sex <- design$samples$sex
C <- cbind(1, sex)
resid_on <- function(M) M - C %*% solve(crossprod(C), crossprod(C, M))
Xs <- scale(resid_on(X * 1.0))
geno <- list(Xc = Xs, xcen = Xs, gram = crossprod(Xs) / n, chr = design$map$CHR, pos = design$map$BP, n = n, p = p)

# KnockoffGWAS partitions: res_1 (1 cM) ... res_7 (0 cM) -> res6 ... res0
part <- rbindlist(lapply(c(21, 22), function(ch) fread(file.path(dir, 'partitions', sprintf('part_chr%d.txt', ch)))))
part <- part[match(design$map$SNP, part$SNP)]
stopifnot(!anyNA(part$SNP))
resolutions <- sprintf('res%d', 6:0)
assignment <- lapply(1:7, function(k) {
  col <- sprintf('res_%d', k); g <- integer(p); off <- 0L
  for (ch in c(21, 22)) { i <- which(geno$chr == ch); loc <- as.integer(factor(part[[col]][i])); g[i] <- loc + off; off <- off + max(loc) }
  g
})
names(assignment) <- resolutions

t0 <- proc.time()[[3]]
operators <- lapply(assignment, function(a) if (max(tabulate(a)) == 1) NULL else group_operators(geno, a))
cat(sprintf('%s: operators for %d resolutions in %.0f s\n', basename(dir), sum(!vapply(operators, is.null, NA)), proc.time()[[3]] - t0))
nodewise <- new_nodewise_cache(geno$xcen)

RANKS <- c('effective', '3', 'full')
um <- function(m) { t <- sqrt(2 * log(m) - 2 * log(log(m))); 2 * (1 - pnorm(t)) }
truncated <- function(pv, q) {          # BH restricted to levels u >= u_m, else nothing
  sel <- bh(pv, q); if (!length(sel)) return(sel)
  if (q * length(sel) / length(pv) >= um(length(pv))) sel else integer(0)
}
ycols <- grep('^y_', names(pheno), value = TRUE)
for (yc in setdiff(ycols, names(done))) {
  y <- as.numeric(resid_on(matrix(pheno[[yc]], ncol = 1)))
  set.seed(1)
  fit <- variant_tests(y, geno$xcen, nodewise)
  residual <- y - as.numeric(geno$Xc %*% fit$beta)
  res_out <- list()
  for (r in resolutions) {
    a <- assignment[[r]]; sizes <- tabulate(a)
    if (is.null(operators[[r]])) {
      P <- matrix(fit$p, p, 4, dimnames = list(NULL, c('simes', paste0('rank_', RANKS)))); K <- matrix(1L, p, 3)
    } else {
      gt <- group_tests(operators[[r]], fit$beta, residual, fit$sigma^2, RANKS)
      P <- cbind(simes = simes(fit$p, a, sizes), gt$p); colnames(P) <- c('simes', paste0('rank_', RANKS)); K <- gt$rank
    }
    res_out[[r]] <- list(p = P, K = K,
      BH = lapply(setNames(colnames(P), colnames(P)), function(t) bh(P[, t], 0.1)),
      truncated = lapply(setNames(colnames(P), colnames(P)), function(t) truncated(P[, t], 0.1)))
  }
  done[[yc]] <- list(sigma = fit$sigma, tested = fit$S, resolutions = res_out)
  saveRDS(done, out_file)
  cat(sprintf('%s %s: %d tested, sigma %.3f  [%.0f s]\n', format(Sys.time(), '%H:%M:%S'), yc, length(fit$S), fit$sigma, proc.time()[[3]] - t0))
}
attr(done, 'assignment') <- assignment
saveRDS(done, out_file)

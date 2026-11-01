# From individual data to summary-statistic DDMM: S = X'X/n, r = X'y/n, sy = y'y/n.
script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
setwd(dirname(normalizePath(script)))
source('ddmm_sumstat_fixed.R')
set.seed(1)
n <- 200; X <- scale(matrix(rnorm(n * 40), n)); y <- drop(scale(2 * X[, 1] + rnorm(n), scale = FALSE))
S <- crossprod(X) / n; r <- drop(crossprod(X, y)) / n; sy <- mean(y^2)
p <- ddmm_pvalues_fixed(S, r, n = n, n_bic = n, sigma = 'estimate', s_y = sy)
print(which(p.adjust(p, 'BH') <= 0.1))

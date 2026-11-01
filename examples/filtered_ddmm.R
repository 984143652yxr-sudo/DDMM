# Run from the repository root: Rscript examples/filtered_ddmm.R
source('R/filters.R')
source('DDMM group/R/ddmm.R')
set.seed(1)
X <- matrix(rnorm(120 * 80), 120, 80)
y <- 2 * X[, 1] + rnorm(120)
f <- spectral_filter(X, y, 'TRIM') # or 'PCA', k = 3; or 'LMM', K = ..., eta = ...
fit <- variant_tests(f$y, f$X, new_nodewise_cache(f$X), noise_cov = f$noise_cov)
selected <- which(p.adjust(fit$p, 'BH') <= 0.1)
print(selected)

# Rebuild or check the per-panel projection caches used by DDMM / SuSiE /
# Knockoff. They depend on genotypes only (never on phenotypes) and already
# exist in data/<panel>/; run this only after the genotype data change.
#
#   Rscript tools/build_projection_cache.R check          verify all panels
#   Rscript tools/build_projection_cache.R build AFR300   build missing caches
#
# cache_NONE.rds : Q removes the intercept (Q = Q_none, n-1 orthonormal rows)
# cache_PCA.rds  : Q removes the intercept and the top 3 PCs (MIX600)
# Each holds X = Q X, Q, XtX, and the nodewise-Lasso residuals Z (BIC choice
# of the penalty for each column) with zx = colSums(Z * X), zz = colSums(Z^2).

script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
setwd(dirname(dirname(script)))

source('R/data.R')
source('R/cache.R')

check_data <- function() {
  for (name in c('AFR300', 'EAS300', 'EUR300', 'MIX300', 'MIX600')) {
    d <- read_panel(name)
    stopifnot(nrow(d$X) == d$n, ncol(d$X) == 1000, !anyDuplicated(d$ids), length(d$betas$SNP) == 1000)
    stopifnot(max(abs(tcrossprod(d$Q_none) - diag(d$n - 1))) < 1e-10)
    cache <- readRDS(file.path('data', name, 'cache_NONE.rds'))
    stopifnot(max(abs(cache$X - d$Q_none %*% d$X)) < 1e-9)
    if (name == 'MIX600') {
      cache <- readRDS(file.path('data', name, 'cache_PCA.rds'))
      stopifnot(max(abs(cache$X - d$Q_pca %*% d$X)) < 1e-9)
      stopifnot(max(abs(d$K - tcrossprod(sweep(d$U, 2, sqrt(d$values), '*')))) < 1e-8)
    }
    cat(name, 'n =', d$n, 'p =', d$p, 'ancestry:',
        paste(names(table(d$ancestry)), table(d$ancestry), collapse = '; '), 'PASS\n')
  }
}

args <- commandArgs(TRUE)
if (!length(args) || args[1] == 'check') check_data() else if (args[1] == 'build' && length(args) == 2)
  build_cache(args[2]) else stop('Usage: Rscript tools/build_projection_cache.R check | build PANEL')

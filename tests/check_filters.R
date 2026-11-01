# Run from repository root. Numerical contracts for the shared filtering interface.
source('R/filters.R'); source('DDMM group/R/ddmm.R')
set.seed(15)
X <- matrix(rnorm(90 * 25), 90); y <- 2 * X[, 1] + rnorm(90)
Z <- scale(X); H <- diag(90) - matrix(1/90, 90, 90)
f <- spectral_filter(X, y, 'PCA', k = 3)
stopifnot(max(abs(f$Q %*% f$Q - f$Q)) < 1e-10,
          abs(sum(diag(f$Q)) - 86) < 1e-8)
f <- spectral_filter(X, y, 'TRIM', rho = 0.3)
d <- svd(Z, nu = 0, nv = 0)$d; cap <- d[floor(0.3 * length(d))]
stopifnot(max(abs(svd(f$X, nu = 0, nv = 0)$d - pmin(d, cap))) < 1e-8)
K <- tcrossprod(Z) / ncol(Z)
f <- spectral_filter(X, y, 'LMM', K = K, eta = 0.4)
stopifnot(max(abs(f$noise_cov - H)) < 1e-8)
# Identity covariance must preserve the existing unfiltered test.
set.seed(18); a <- variant_tests(y, Z, new_nodewise_cache(Z))
set.seed(18); b <- variant_tests(y, Z, new_nodewise_cache(Z), diag(nrow(Z)))
stopifnot(max(abs(a$p - b$p)) < 1e-10, abs(a$sigma - b$sigma) < 1e-10)
# Check the reported statistic against z' C z after a non-idempotent filter.
f <- spectral_filter(X, y, 'TRIM', rho = 0.3)
cache <- new_nodewise_cache(f$X)
fit <- variant_tests(f$y, f$X, cache, f$noise_cov)
j <- fit$S[1]; nw <- cache$fits[[j]]
z <- f$X[, j] - drop(f$X[, -j] %*% nw$g)
r <- f$y - drop(f$X %*% fit$beta)
estimate <- fit$beta[j] + sum(z * r) / sum(z * f$X[, j])
expected <- estimate * sum(z * f$X[, j]) / (fit$sigma * sqrt(sum(z * (f$noise_cov %*% z))))
stopifnot(abs(expected - fit$statistic[j]) < 1e-8)
cat('PASS: PCA projection, TRIM cap, LMM whitening, default compatibility, filtered variance.\n')

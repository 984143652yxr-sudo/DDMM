# Spectral filters for y = X beta + noise. All filters remove the intercept.
# X is centred and column-standardised; coefficients refer to that scale.
# PCA/TRIM follow the spectral construction in the original simulation scripts.
# LMM uses covariance eigenvalues and inverse square roots (whitening).
spectral_filter <- function(X, y, method = c('PCA', 'TRIM', 'LMM'),
                            k = 3L, rho = 0.3, K = NULL, eta = 0.5,
                            noise_cov = NULL) {
  method <- match.arg(method)
  X <- as.matrix(X); n <- nrow(X)
  stopifnot(n >= 6, ncol(X) >= 2, length(y) == n, all(is.finite(X)), all(is.finite(y)))
  Z <- scale(X)
  if (any(!is.finite(Z))) stop('Remove constant genotype columns before filtering.')
  H <- diag(n) - matrix(1 / n, n, n)
  if (method == 'LMM') {
    if (is.null(K)) K <- tcrossprod(Z) / ncol(Z)
    stopifnot(all(dim(K) == c(n, n)), eta >= 0, eta < 1, all(is.finite(K)))
    if (max(abs(K - t(K))) > 1e-8 || min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) < -1e-8)
      stop('K must be symmetric positive semidefinite.')
    C <- eta * K + (1 - eta) * diag(n)
    eg <- eigen(C, symmetric = TRUE)
    F <- tcrossprod(sweep(eg$vectors, 2, 1 / sqrt(eg$values), '*'), eg$vectors)
    if (is.null(noise_cov)) noise_cov <- C
  } else {
    s <- svd(Z, nv = 0)
    rank <- sum(s$d > max(s$d) * 1e-10)
    if (method == 'PCA') {
      stopifnot(length(k) == 1, k == as.integer(k), k >= 0, k < rank)
      U <- s$u[, seq_len(k), drop = FALSE]
      F <- diag(n) - tcrossprod(U)
    } else {
      stopifnot(rho > 0, rho <= 1)
      j <- max(1L, floor(rho * rank))
      U <- s$u[, seq_len(rank), drop = FALSE]
      w <- pmin(1, s$d[j] / s$d[seq_len(rank)])
      F <- diag(n) - tcrossprod(sweep(U, 2, 1 - w, '*'), U)
    }
  }
  if (is.null(noise_cov)) noise_cov <- diag(n)
  stopifnot(all(dim(noise_cov) == c(n, n)))
  Q <- F %*% H
  list(X = Q %*% Z, y = drop(Q %*% y), Q = Q,
       noise_cov = Q %*% noise_cov %*% t(Q),
       center = attr(Z, 'scaled:center'), scale = attr(Z, 'scaled:scale'))
}

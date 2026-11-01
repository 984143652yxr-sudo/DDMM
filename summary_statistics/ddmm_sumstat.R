## ddmm_sumstat.R -- summary-statistic DDMM (Section 5)
##
## Everything below depends on the data only through
##
##     Sigma = t(X) %*% X / n     LD matrix
##     r     = t(X) %*% y / n     marginal score
##
## so no individual-level genotype is needed. Marginal z-scores are the case
## n = 1, since z ~ N(Sigma beta, Sigma) when Sigma is an LD correlation matrix.
## Sigma is a correlation matrix in all of our analyses, so the penalty weights
## D of Eq (5.3) equal the identity and the statistics are on the unit-noise
## scale (sigma = 1), as in the remark closing Section 5.

soft <- function(x, lambda) sign(x) * pmax(abs(x) - lambda, 0)

## FISTA for  min 1/2 b' A b - c' b + lambda |b|_1.  Both the initial fit
## Eq (5.3) and the nodewise regressions Eq (5.4) have this form.
lasso_gram <- function(A, c, lambda, step, iters = 300) {
  b <- numeric(length(c))
  u <- b
  t <- 1
  for (i in seq_len(iters)) {
    b_new <- soft(u - (A %*% u - c) / step, lambda / step)
    t_new <- (1 + sqrt(1 + 4 * t^2)) / 2
    u <- b_new + ((t - 1) / t_new) * (b_new - b)
    b <- as.numeric(b_new)
    t <- t_new
  }
  b
}

## PCA filter of Eq (5.2): remove the leading q eigendirections of Sigma.
deflate <- function(Sigma, r, q) {
  if (q == 0) return(list(Sigma = Sigma, r = r))
  eig <- eigen(Sigma, symmetric = TRUE)
  V <- eig$vectors[, seq_len(q), drop = FALSE]
  mu <- eig$values[seq_len(q)]
  list(Sigma = Sigma - V %*% (mu * t(V)),
       r     = r - as.vector(V %*% crossprod(V, r)))
}

## DDMM p-values from (Sigma, r): initial fit, nodewise debiasing of the active
## set, then the studentized statistics of Eq (5.7). Coordinates outside the
## active set are not tested and keep a p-value of one.
ddmm_pvalues <- function(Sigma, r, n = 1, q = 0) {
  filtered <- deflate(Sigma, r, q)
  Sigma <- filtered$Sigma
  r <- filtered$r

  p <- ncol(Sigma)
  step <- max(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  beta <- lasso_gram(Sigma, r, 1.1 * sqrt(2 * log(p) / n), step)

  lambda_node <- sqrt(2 * log(p) / n)
  pv <- rep(1, p)
  for (j in which(abs(beta) > 1e-6)) {
    rest <- seq_len(p)[-j]
    A <- Sigma[rest, rest]
    gamma <- lasso_gram(A, Sigma[rest, j], lambda_node, step)

    d <- Sigma[j, j] - sum(Sigma[j, rest] * gamma)                    # Eq (5.5)
    v <- Sigma[j, j] - 2 * sum(Sigma[rest, j] * gamma) + sum(gamma * (A %*% gamma))
    g <- r[j] - sum(gamma * r[rest]) -
      sum((Sigma[j, ] - as.vector(crossprod(gamma, Sigma[rest, ]))) * beta)

    ## beta_debiased - beta_j has standard deviation sigma sqrt(v/n)/d, so the
    ## studentized statistic divides by that. The nodewise residual satisfies
    ## d = v only when gamma is the unpenalized projection; under dense LD the
    ## L1 penalty makes d/v appreciably larger than one, and replacing d/sqrt(v)
    ## by sqrt(v) as in Section 3.2 then costs most of the power.
    beta_debiased <- beta[j] + g / d
    T <- sqrt(n) * beta_debiased * d / sqrt(max(v, 1e-8))             # Eq (5.7)
    pv[j] <- 2 * pnorm(abs(T), lower.tail = FALSE)
  }
  pv
}

## Benjamini-Hochberg, applied after pooling p-values across LD blocks.
bh_select <- function(pv, alpha) {
  ord <- order(pv)
  below <- pv[ord] <= alpha * seq_along(pv) / length(pv)
  if (any(below)) ord[seq_len(max(which(below)))] else integer(0)
}

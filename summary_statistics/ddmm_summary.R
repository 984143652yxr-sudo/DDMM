# Corrected summary-statistic DDMM with per-SNP BIC nodewise tuning.

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

nodewise_bic_gram <- function(Sigma, S, n_bic, nlambda = 30L, eps = 0.01) {
  p  <- ncol(Sigma)
  eg <- eigen(Sigma, symmetric = TRUE)
  R  <- sqrt(pmax(eg$values, 0)) * t(eg$vectors)          # R'R = Sigma
  out <- vector("list", length(S))
  for (k in seq_along(S)) {
    j   <- S[k]
    Smj <- Sigma[-j, j]
    lam_max <- max(abs(Smj)) / p
    if (!is.finite(lam_max) || lam_max <= 0) lam_max <- 1e-3
    lam <- lam_max * exp(seq(0, log(eps), length.out = nlambda))
    fit <- glmnet::glmnet(R[, -j, drop = FALSE], R[, j], lambda = lam, intercept = FALSE)
    B   <- as.matrix(fit$beta)
    rss <- Sigma[j, j] - 2 * as.numeric(crossprod(B, Smj)) +
      colSums(B * (Sigma[-j, -j] %*% B))
    b   <- which.min(n_bic * log(pmax(rss, .Machine$double.eps)) + log(n_bic) * fit$df)
    g   <- B[, b]
    out[[k]] <- list(gamma = g,
                     d = Sigma[j, j] - sum(Smj * g),
                     v = max(rss[b], 1e-8))
  }
  out
}

## Corrected DDMM p-values from (Sigma, r, n). Additional arguments:
##   n_bic  sample size used in the nodewise BIC (FIX 1)
##   sigma  "unit" (sigma = 1, z-score convention) or "estimate" (FIX 3)
##   s_y    y'y / n, needed only for sigma = "estimate"
##   cache  optional environment for nodewise fits when Sigma (and q) are
##          the same across calls; never share one across different Sigma
## Returns the p-values; attributes carry the statistic, sigma and the support.
ddmm_pvalues <- function(Sigma, r, n = 1, q = 0, n_bic = n,
                               sigma = c("unit", "estimate"), s_y = 1,
                               cache = NULL) {
  sigma <- match.arg(sigma)
  if (q > 0) {                               # s_y is deflated with the same q
    eig <- eigen(Sigma, symmetric = TRUE)
    V <- eig$vectors[, seq_len(q), drop = FALSE]; mu <- eig$values[seq_len(q)]
    s_y <- s_y - sum(as.numeric(crossprod(V, r))^2 / mu)
  }
  filtered <- deflate(Sigma, r, q)
  Sigma <- filtered$Sigma
  r <- filtered$r

  p <- ncol(Sigma)
  step <- max(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  beta <- lasso_gram(Sigma, r, 1.1 * sqrt(2 * log(p) / n), step)   # unchanged
  S <- which(abs(beta) > 1e-6)

  sig <- 1
  if (sigma == "estimate") {                 # FIX 3: df-corrected denominator
    rss_n <- s_y - 2 * sum(beta * r) + sum(beta * (Sigma %*% beta))
    sig <- sqrt(max(n * rss_n, 1e-12) / max(n - q - length(S), 1))
  }

  pv <- rep(1, p); Tst <- rep(0, p)
  if (length(S) > 0) {
    ## FIX 1. The nodewise fits depend on Sigma only, so a caller that reuses
    ## one Sigma across replicates can pass an environment as `cache`.
    if (is.null(cache)) {
      nw <- nodewise_bic_gram(Sigma, S, n_bic)
    } else {
      key <- as.character(S)
      todo <- S[!vapply(key, exists, logical(1), envir = cache, inherits = FALSE)]
      if (length(todo) > 0) {
        new <- nodewise_bic_gram(Sigma, todo, n_bic)
        for (k in seq_along(todo)) assign(as.character(todo[k]), new[[k]], envir = cache)
      }
      nw <- mget(key, envir = cache)
    }
    for (k in seq_along(S)) {
      j <- S[k]; rest <- seq_len(p)[-j]; gamma <- nw[[k]]$gamma
      d <- nw[[k]]$d; v <- nw[[k]]$v
      g <- r[j] - sum(gamma * r[rest]) -
        sum((Sigma[j, ] - as.vector(crossprod(gamma, Sigma[rest, ]))) * beta)
      beta_debiased <- beta[j] + g / d
      Tst[j] <- sqrt(n) * beta_debiased * d / (sig * sqrt(v))
      pv[j] <- 2 * pnorm(abs(Tst[j]), lower.tail = FALSE)
    }
  }
  attr(pv, "T") <- Tst; attr(pv, "sigma") <- sig; attr(pv, "support") <- S
  pv
}

## Benjamini-Hochberg, applied after pooling p-values across LD blocks.
bh_select <- function(pv, alpha) {
  ord <- order(pv)
  below <- pv[ord] <= alpha * seq_along(pv) / length(pv)
  if (any(below)) ord[seq_len(max(which(below)))] else integer(0)
}

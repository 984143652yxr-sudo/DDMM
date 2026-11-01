## ghostknockoff.R -- GhostKnockoff for comparison, following the reference
## implementation at
##   https://github.com/GuJQ5/FWER_Knockoff_filter/blob/main/Functions/GhostKnockoff_hello.R
##
## The three functions below mirror GhostKnockoff.prelim, GhostKnockoff.fit with
## method = 'marginal', and GhostKnockoff.filter. Notation follows the reference:
## diag_s, P.each, Sigma_k, V.each, V.left, kappa, tau.
##
## Two constructions of the knockoff contrast s are provided, matching the
## reference. "equi" gives every variant the same s and is bounded by the
## smallest eigenvalue of Sigma, so it loses contrast as the LD densifies.
## "sdp" solves for a per-variant s and is the reference default; it is
## typically more powerful under dense LD but requires the Rdsdp package.

is_posdef <- function(A, tol = 1e-9) {
  min(eigen(A, symmetric = TRUE, only.values = TRUE)$values) > tol * 10
}

## Shrink s until (M+1)/M Sigma - diag(s) stays positive definite.
shrink_to_psd <- function(s, G, M) {
  s_eps <- 1e-8
  while (s_eps <= 0.1 && !is_posdef((M + 1) / M * G - diag(s * (1 - s_eps), length(s)))) {
    s_eps <- s_eps * 10
  }
  s * (1 - s_eps)
}

## create.solve_equi_M: s = min((M+1)/M lambda_min, 1) for every variant.
solve_equi_M <- function(Sigma, M) {
  G <- cov2cor(Sigma)
  lambda_min <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
  s <- rep(min((M + 1) / M * lambda_min, 1), nrow(G))
  shrink_to_psd(s, G, M) * diag(Sigma)
}

## create.solve_sdp_M: maximise sum(s) subject to 0 <= s <= 1 and
## (M+1)/M G - diag(s) positive semidefinite.
solve_sdp_M <- function(Sigma, M, gaptol = 1e-6, maxit = 1000) {
  if (!requireNamespace("Rdsdp", quietly = TRUE)) {
    stop("method = 'sdp' needs the Rdsdp package; install it or use method = 'equi'.")
  }
  G <- cov2cor(Sigma)
  p <- ncol(G)
  if (!is_posdef(G)) {
    warning("LD matrix is not positive definite; knockoffs may lose power.", immediate. = TRUE)
  }

  ## Linear constraints s >= 0 and s <= 1, then the positive-definite cone.
  Al1 <- -Matrix::Diagonal(p)
  Al2 <- Matrix::Diagonal(p)
  d_As <- c(diag(p))
  As <- Matrix::Diagonal(length(d_As), x = d_As)
  As <- As[which(Matrix::rowSums(As) > 0), ]

  A <- cbind(Al1, Al2, As)
  C <- matrix(c(rep(0, p), rep(1, p), c((M + 1) / M * G)), 1)
  K <- list(s = p, l = 2 * p)

  sol <- Rdsdp::dsdp(A, rep(1, p), C, K,
                     list(gaptol = gaptol, maxit = maxit, logsummary = 0,
                          outputstats = 0, print = 0))
  if (!identical(sol$STATS$stype, "PDFeasible")) {
    warning("SDP solver returned a non-feasible solution; knockoffs may lose power.")
  }

  s <- pmin(pmax(sol$y, 0), 1)
  s <- shrink_to_psd(s, G, M)
  s[s < 0] <- 0
  if (all(s == 0)) {
    warning("SDP construction failed; knockoffs will have no power.", immediate. = TRUE)
  }
  s * diag(Sigma)
}

## GhostKnockoff.prelim
ghostknockoff_prelim <- function(Sigma, M, method = c("equi", "sdp")) {
  method <- match.arg(method)
  p <- nrow(Sigma)
  SigmaInv <- solve(Sigma)
  s <- if (method == "equi") solve_equi_M(Sigma, M) else solve_sdp_M(Sigma, M)
  diag_s <- diag(s, p)

  Sigma_k <- 2 * diag_s - diag_s %*% SigmaInv %*% diag_s
  V.each <- as.matrix(forceSymmetric(Sigma_k - diag_s))
  V <- matrix(1, M, M) %x% V.each + diag(1, M) %x% diag_s
  V.left <- t(chol(as.matrix(forceSymmetric(V)) + diag(1e-8, M * p)))

  list(P.each = diag(1, p) - diag_s %*% SigmaInv, V.left = V.left, M = M)
}

## GhostKnockoff.fit, method = 'marginal': one draw of M knockoff z-score
## vectors, returned as a p x M matrix.
ghostknockoff_zscores <- function(z, prelim) {
  p <- length(z)
  noise <- matrix(prelim$V.left %*% rnorm(ncol(prelim$V.left)), nrow = p)
  as.vector(prelim$P.each %*% z) + noise
}

## GhostKnockoff.filter: multiple-knockoff statistics and q-values.
ghostknockoff_filter <- function(T_0, T_k) {
  T.temp <- cbind(T_0, T_k)
  T.temp[is.na(T.temp)] <- 0
  M <- ncol(T_k)

  ## Ties go to the knockoff, as in the reference which.max.alt.
  which.max.alt <- function(x) {
    i <- which(x == max(x))
    if (length(i) != 1) i[2] else i[1]
  }
  kappa <- apply(T.temp, 1, which.max.alt) - 1
  tau <- apply(T.temp, 1, max) -
    apply(T.temp, 1, function(x) median(x[-which.max(x)]))

  b <- order(tau, decreasing = TRUE)
  c_0 <- kappa[b] == 0
  n_orig <- cumsum(c_0)                                   # temp_0
  ratio <- (1 + (seq_along(b) - n_orig)) / (M * pmax(1, n_orig))

  q <- rep(1, length(tau))
  if (any(tau[b] > 0)) {
    bound <- max(which(tau[b] > 0))
    for (i in seq_along(b)) {
      idx <- i:min(length(b), bound)
      if (length(idx) > 0) q[b[i]] <- min(ratio[idx]) * c_0[i] + 1 - c_0[i]
    }
    q[q > 1] <- 1
  }
  list(kappa = kappa, tau = tau, q = q)
}

# DDMM single-variant and group tests (the corrected version of 25 Sep 2026;
# see docs/METHOD.md).
#
#   pilot     beta_hat from a 5-fold cross-validated Lasso (lambda.min), the
#             same fit for every test (FIX 2)
#   sigma2    RSS / ||(I - P_S)||_F^2, i.e. RSS / (n - rank(X_S)), corrected
#             for the degrees of freedom the Lasso spends (FIX 3)
#   variant   b_j = beta_j + z_j'r / z_j'x_j,  T_j = b_j z_j'x_j / (sigma ||z_j||),
#             z_j the nodewise residual of x_j on all other variants with its
#             own lambda chosen by BIC (FIX 1)
#   group     b_g = beta_g + (Z_g'X_g)^+ Z_g'r, Z_g = X_g - X_N Gamma_g from a
#             nodewise group Lasso on the +/-40 neighbouring variants,
#             Omega_g = (Z_g'X_g)^+ Z_g'Z_g (Z_g'X_g)^+',
#             chi-squared Wald statistic on the top K within-group LD
#             directions (rank K) or on all |g| directions (full rank)

NODEWISE_PATH <- list(nlambda = 30L, eps = 0.01)      # per-variant nodewise path
GROUP_NODEWISE <- list(neighbourhood = 40L, iterations = 90L)

# ---------------------------------------------------------------------------
# Rank of the reduced-rank group statistic
#   'effective'  participation ratio of the group's LD eigenvalues,
#                K_g = round((sum l)^2 / sum l^2), bounded to [1, |g|]:
#                1 for perfectly correlated variants, |g| for independent ones
#   integer K    the top K directions, K_g = min(K, |g|)
#   'full'       no reduction: all |g| directions, chi-squared on |g| df
# ---------------------------------------------------------------------------

effective_rank <- function(ev) {
  ev <- pmax(ev, 0)
  if (sum(ev) <= 0) return(1L)
  max(1L, min(length(ev), as.integer(round(sum(ev)^2 / sum(ev^2)))))
}

parse_ranks <- function(spec) {
  ranks <- unique(tolower(strsplit(gsub('\\s', '', spec), ',', fixed = TRUE)[[1]]))
  ok <- ranks %in% c('effective', 'full') | grepl('^[1-9][0-9]*$', ranks)
  if (!all(ok)) stop('Unknown rank "', ranks[!ok][1], '". Use effective, full or a positive integer such as 3.')
  ranks
}

rank_label <- function(rank) switch(rank,
  effective = 'chi2, effective rank', full = 'chi2, full rank', paste0('chi2, K = ', rank))

track_id <- function(rank) paste0('rank_', rank)

group_rank <- function(rank, eigenvalues) {
  m <- length(eigenvalues)
  if (rank == 'full') m else if (rank == 'effective') min(effective_rank(eigenvalues), m)
  else min(as.integer(rank), m)
}

# ---------------------------------------------------------------------------
# Per-variant nodewise regressions (FIX 1), cached
# The regression of x_j on all other variants depends on the genotypes only.
# The cache keeps, for every variant needed so far, the BIC-chosen coefficient
# vector g_j and z_j'z_j, z_j'x_j, computed from the Gram matrix exactly as in
# fast_nodewise.R (Gram eigendecomposition, 30-point path).
# ---------------------------------------------------------------------------

new_nodewise_cache <- function(xcen) {
  G <- crossprod(xcen)
  eg <- eigen(G, symmetric = TRUE)
  cache <- new.env()
  cache$G <- G
  cache$R <- sqrt(pmax(eg$values, 0)) * t(eg$vectors)
  cache$n <- nrow(xcen)
  cache$fits <- vector('list', ncol(G))
  cache
}

nodewise_fit <- function(cache, j) {
  G <- cache$G; R <- cache$R; p <- ncol(G); n <- cache$n
  Gmj <- G[-j, j]
  lam_max <- max(abs(Gmj)) / p
  if (!is.finite(lam_max) || lam_max <= 0) lam_max <- 1e-3
  lam <- lam_max * exp(seq(0, log(NODEWISE_PATH$eps), length.out = NODEWISE_PATH$nlambda))
  fit <- glmnet::glmnet(R[, -j, drop = FALSE], R[, j], lambda = lam, intercept = FALSE)
  B <- as.matrix(fit$beta)
  rss <- G[j, j] - 2 * as.numeric(crossprod(B, Gmj)) + colSums(B * (G[-j, -j] %*% B))
  b <- which.min(n * log(pmax(rss, .Machine$double.eps) / n) + log(n) * fit$df)
  g <- B[, b]
  list(g = g, zz = max(rss[b], 0), zx = G[j, j] - sum(g * Gmj))
}

ensure_nodewise <- function(cache, S) {
  missing <- S[vapply(cache$fits[S], is.null, logical(1))]
  for (j in missing) cache$fits[[j]] <- nodewise_fit(cache, j)
  length(missing)
}

# ---------------------------------------------------------------------------
# Per-variant debiased Lasso: pilot, sigma and p-values
# ---------------------------------------------------------------------------

variant_tests <- function(y, xcen, nodewise) {
  lam <- glmnet::cv.glmnet(xcen, y, intercept = FALSE, standardize = FALSE, nfolds = 5)$lambda.min
  beta <- as.numeric(glmnet::glmnet(xcen, y, lambda = lam, intercept = FALSE,
                                    standardize = FALSE)$beta[, 1])
  S <- which(beta != 0)
  r <- as.numeric(y - xcen %*% beta)
  # Guard: with fewer than 3 selected variants (common under the null), test
  # the top marginal variants as well, as the shared DDMM core does.
  if (length(S) < 3) S <- marginal_screen(y, xcen, S, min_selected = 3, top = 20)

  # FIX 3: sigma^2 = RSS / ||(I - P_S)||_F^2
  sigma <- if (length(S)) {
    sv <- svd(xcen[, S, drop = FALSE], nv = 0)
    P <- diag(1, nrow(xcen)) - tcrossprod(sv$u[, seq_len(sum(sv$d > 1e-10)), drop = FALSE])
    sqrt(sum(r^2) / sum(P^2))
  } else sqrt(sum(r^2) / nrow(xcen))

  statistic <- numeric(ncol(xcen)); p_value <- rep(1, ncol(xcen))
  if (length(S)) {
    computed <- ensure_nodewise(nodewise, S)
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

# BASIL-style marginal screen (marginal_screen_signal of the shared DDMM core).
marginal_screen <- function(y, x, selected, min_selected = 3, top = 20) {
  p <- ncol(x)
  selected <- sort(unique(as.integer(selected)))
  min_selected <- min(max(1, min_selected), p)
  top <- min(p, max(min_selected, top))
  if (length(selected) >= min_selected) return(selected)
  score <- abs(as.numeric(crossprod(x, as.numeric(y)))) /
    (sqrt(colSums(x^2)) * sqrt(sum(y^2)) + 1e-15)
  score[!is.finite(score)] <- 0
  sort(unique(c(selected, order(score, decreasing = TRUE)[seq_len(top)])))
}

# ---------------------------------------------------------------------------
# Group operators (genotype-only, one set per partition)
# ---------------------------------------------------------------------------

pseudo_inverse <- function(a, tol = 1e-10) {
  s <- svd(a)
  keep <- s$d > tol * max(s$d)
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# Nodewise group Lasso Gamma_g of X_g on its neighbours X_N, by FISTA:
#   minimise tr(G' gram_NN G)/2 - tr(G' gram_Ng) + lambda sum_b sqrt(|b|) ||G_b||_F
# The penalty acts on whole neighbouring groups.
nodewise_gamma <- function(gram_nn, gram_ng, neighbour_group, lambda,
                           iter = GROUP_NODEWISE$iterations) {
  step <- max(eigen(gram_nn, symmetric = TRUE, only.values = TRUE)$values) + 1e-6
  gamma <- matrix(0, nrow(gram_ng), ncol(gram_ng)); z <- gamma; t <- 1
  blocks <- unique(neighbour_group)
  for (i in seq_len(iter)) {
    u <- z - (gram_nn %*% z - gram_ng) / step
    next_gamma <- u
    for (b in blocks) {
      rows <- neighbour_group == b
      shrink <- lambda * sqrt(sum(rows)) / step
      block_norm <- sqrt(sum(u[rows, , drop = FALSE]^2))
      next_gamma[rows, ] <- max(0, 1 - shrink / (block_norm + 1e-12)) * u[rows, , drop = FALSE]
    }
    next_t <- (1 + sqrt(1 + 4 * t * t)) / 2
    z <- next_gamma + ((t - 1) / next_t) * (next_gamma - gamma)
    gamma <- next_gamma; t <- next_t
  }
  gamma
}

group_operators <- function(geno, assignment) {
  p <- geno$p; n <- geno$n
  members_by_group <- split(seq_len(p), assignment)
  lambda <- sqrt(log(p) / n)
  nbhd <- GROUP_NODEWISE$neighbourhood
  lapply(members_by_group, function(members) {
    ch <- geno$chr[members[1]]
    neighbours <- setdiff(which(geno$chr == ch & seq_len(p) >= min(members) - nbhd &
                                  seq_len(p) <= max(members) + nbhd), members)
    Xg <- geno$Xc[, members, drop = FALSE]
    Zg <- if (!length(neighbours)) Xg else {
      gamma <- nodewise_gamma(geno$gram[neighbours, neighbours, drop = FALSE],
                              geno$gram[neighbours, members, drop = FALSE],
                              assignment[neighbours], lambda)
      Xg - geno$Xc[, neighbours, drop = FALSE] %*% gamma
    }
    inverse <- pseudo_inverse(crossprod(Zg, Xg))
    eig <- eigen(geno$gram[members, members, drop = FALSE], symmetric = TRUE)
    list(members = members, Z = Zg, inverse = inverse,
         omega = inverse %*% crossprod(Zg) %*% t(inverse),
         eigenvalues = eig$values, eigenvectors = eig$vectors)
  })
}

# ---------------------------------------------------------------------------
# Group tests for one case at one resolution
# ---------------------------------------------------------------------------

# Simes combination of the tested variants' p-values within each group;
# untested variants count as p = 1 through the group size.
simes <- function(p_variant, assignment, group_size) {
  out <- rep(1, length(group_size))
  tested <- which(p_variant < 1)
  by_group <- split(tested, assignment[tested])
  for (label in names(by_group)) {
    g <- as.integer(label)
    q <- sort(p_variant[by_group[[label]]])
    out[g] <- min(1, group_size[g] * q / seq_along(q))
  }
  out
}

group_tests <- function(operators, beta_pilot, residual, sigma2, ranks) {
  out <- matrix(1, length(operators), length(ranks), dimnames = list(NULL, ranks))
  used <- matrix(NA_integer_, length(operators), length(ranks), dimnames = list(NULL, ranks))
  for (g in seq_along(operators)) {
    op <- operators[[g]]
    b <- beta_pilot[op$members] + as.numeric(op$inverse %*% crossprod(op$Z, residual))
    for (rank in ranks) {
      k <- group_rank(rank, op$eigenvalues)
      if (rank == 'full') {
        stat <- as.numeric(t(b) %*% pseudo_inverse(op$omega) %*% b) / sigma2
      } else {
        C <- op$eigenvectors[, seq_len(k), drop = FALSE]
        pb <- as.numeric(crossprod(C, b))
        stat <- as.numeric(t(pb) %*% pseudo_inverse(t(C) %*% op$omega %*% C) %*% pb) / sigma2
      }
      out[g, rank] <- stats::pchisq(stat, df = k, lower.tail = FALSE)
      used[g, rank] <- k
    }
  }
  list(p = out, rank = used)
}

# Benjamini-Hochberg selection at level q (indices of selected groups).
bh <- function(p_value, q) {
  o <- order(p_value)
  below <- p_value[o] <= q * seq_along(p_value) / length(p_value)
  if (any(below)) sort(o[seq_len(max(which(below)))]) else integer(0)
}

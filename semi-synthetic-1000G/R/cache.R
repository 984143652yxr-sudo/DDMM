# Genotype-only caches, rebuilt locally from panel.rds.
make_nodewise_cache <- function(X, Q) {
  n <- nrow(X)
  p <- ncol(X)
  residuals <- lapply(seq_len(p), function(j) {
    target <- X[, j]
    nuisance <- X[, -j, drop = FALSE]
    upper <- max(abs(crossprod(nuisance, target))) / p
    if (!is.finite(upper) || upper <= 0) upper <- 1e-3
    lambda <- upper * exp(seq(0, log(.01), length.out = 30)) * sqrt(p / n)
    fit <- glmnet::glmnet(nuisance, target, lambda = lambda, intercept = FALSE)
    errors <- target - as.matrix(nuisance %*% fit$beta)
    rss <- colSums(errors^2)
    bic <- n * log(pmax(rss, .Machine$double.eps) / n) + log(n) * fit$df
    errors[, which.min(bic)]
  })
  Z <- do.call(cbind, residuals)
  result <- list(X = X, Q = Q, XtX = crossprod(X), Z = Z, zx = colSums(Z * X), zz = colSums(Z^2))
  stopifnot(all(is.finite(Z)), min(abs(result$zx)) > 1e-8)
  result
}

build_cache <- function(panel, projection = NULL) {
  path <- file.path('data', panel)
  d <- read_panel(panel)
  destination <- file.path(path, 'cache_NONE.rds')
  if ((is.null(projection) || projection == 'NONE') && !file.exists(destination))
    saveRDS(make_nodewise_cache(d$Q_none %*% d$X, d$Q_none), destination)
  if (file.exists(destination)) cat('cache_NONE ready:', panel, '\n')
  if (!is.null(d$Q_pca) && (is.null(projection) || projection == 'PCA')) {
    destination <- file.path(path, 'cache_PCA.rds')
    if (!file.exists(destination)) saveRDS(make_nodewise_cache(d$Q_pca %*% d$X, d$Q_pca), destination)
    cat('cache_PCA ready:', panel, '\n')
  }
}

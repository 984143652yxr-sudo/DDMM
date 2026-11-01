# The four selection procedures compared in both experiments, and scoring.
#
# DDMM, SuSiE and Knockoff+ are fitted to the transformed data (Q y, Q X)
# held in the projection cache: Q removes the intercept (projection NONE,
# ancestry experiment) or the intercept and the top 3 genotype PCs
# (projection PCA, kinship600). BOLT-LMM is run on the raw genotypes and raw
# phenotype, with PC1-3 as covariates in the kinship experiment.
# Every fit returns list(selected = <SNP indices>, valid = TRUE/FALSE, error).

# Shared CV-Lasso pilot used by DDMM (debiasing) and SuSiE (initialisation).
fit_pilot <- function(X, y, seed) {
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  lambda <- sqrt(mean(y^2)) * sqrt(log(p) / n) * (100:1) / 100 * 3
  fit <- glmnet::cv.glmnet(X, y, lambda = lambda, intercept = FALSE,
                           family = 'gaussian', nfolds = min(10, n))
  beta <- as.numeric(stats::coef(fit, s = 'lambda.min'))[-1]
  residual <- drop(y - X %*% beta)
  active <- which(beta != 0)
  rank <- if (length(active)) qr(X[, active, drop = FALSE])$rank else 0
  stopifnot(n > rank)
  list(beta = beta, residual = residual, sigma2 = sum(residual^2) / (n - rank), df = n - rank)
}

# DDMM: coordinatewise debiasing with nodewise-Lasso residuals z_j
#   b_j = beta_j + z_j'(y - X beta) / z_j'x_j,
#   se_j = sqrt(sigma2 * z_j'z_j) / |z_j'x_j|,
# two-sided normal p-values, Benjamini-Hochberg at level q.
fit_ddmm <- function(cache, pilot, q = .1) {
  correction <- drop(crossprod(cache$Z, pilot$residual)) / cache$zx
  estimate <- pilot$beta + correction
  se <- sqrt(pilot$sigma2 * cache$zz) / abs(cache$zx)
  pvalue <- 2 * pnorm(-abs(estimate / se))
  stopifnot(all(is.finite(pvalue)))
  list(selected = which(p.adjust(pvalue, 'BH') <= q),
       pvalue = pvalue, estimate = estimate, valid = TRUE, error = '')
}

# SuSiE (sufficient statistics), L effects, select PIP > cutoff.
# Four starts (default / pilot-initialised, each with estimated or pilot
# residual variance); the start with the highest final ELBO is kept.
fit_susie <- function(cache, y, pilot, L = 40, cutoff = .9) {
  X <- cache$X
  p <- ncol(X)
  initial_indices <- head(order(abs(pilot$beta), decreasing = TRUE), L)
  initial_indices <- initial_indices[pilot$beta[initial_indices] != 0]
  initial <- if (length(initial_indices))
    susieR::susie_init_coef(initial_indices, pilot$beta[initial_indices], p) else NULL
  # Mute only susieR's note that a pilot start with < L effects is padded to L.
  quiet <- function(expr) withCallingHandlers(expr, message = function(m)
    if (grepl('is greater the number of effects', conditionMessage(m))) invokeRestart('muffleMessage'))
  fits <- lapply(1:4, function(start) quiet(
    susieR::susie_suff_stat(
      XtX = cache$XtX, Xty = drop(crossprod(X, y)), yty = sum(y^2), n = nrow(X),
      X_colmeans = rep(0, p), y_mean = 0, L = L,
      residual_variance = if (start > 2) pilot$sigma2 else NULL,
      estimate_residual_variance = TRUE,
      s_init = if (start %% 2 == 0) initial else NULL,
      max_iter = 1000, tol = 1e-4)))
  objective <- vapply(fits, function(f) tail(f$elbo, 1), numeric(1))
  best <- fits[[which.max(objective)]]
  probability <- best$pip
  valid <- isTRUE(best$converged) && all(is.finite(probability))
  list(selected = which(probability > cutoff), probability = probability,
       valid = valid, error = if (valid) '' else 'SuSiE did not converge',
       best_start = which.max(objective), elbo = objective)
}

# Knockoff+ with SNP-level HMM knockoff copies (see R/knockoffs.R).
# Importance: Lasso coefficient difference on the genotype scale,
# W_j = |b_j| - |b_j~|, with the Knockoff+ threshold (offset = 1) at level q.
fit_knockoff <- function(augmented, y, seed, q = .1) {
  set.seed(seed)
  p <- ncol(augmented) / 2
  permutation <- sample.int(2 * p)
  folds <- sample(rep(1:5, length.out = nrow(augmented)))
  fit <- glmnet::cv.glmnet(augmented[, permutation], y, foldid = folds,
                           nlambda = 100, intercept = FALSE, standardize = TRUE)
  coefficients <- numeric(2 * p)
  coefficients[permutation] <- as.numeric(coef(fit, s = 'lambda.min'))[-1]
  coefficients <- coefficients * sqrt(colMeans(augmented^2))
  W <- abs(coefficients[seq_len(p)]) - abs(coefficients[p + seq_len(p)])
  threshold <- knockoff::knockoff.threshold(W, fdr = q, offset = 1)
  list(selected = which(W >= threshold), W = W, threshold = threshold,
       valid = TRUE, error = '')
}

# BOLT-LMM-inf association p-values (--lmmInfOnly) on the raw data, then BH
# at level q over the tested SNPs. A BOLT failure (e.g. "heritability
# estimate is close to 1" at n = 300) is recorded as an invalid fit; its logs
# stay in <case>/bolt/.
fit_bolt <- function(panel, y, use_pc_covariates, directory, q = .1) {
  executable <- bolt_path()
  if (!nzchar(executable)) stop('Set BOLT_BIN to an installed BOLT-LMM executable.')
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  phenotype_file <- file.path(directory, 'phenotype.txt')
  utils::write.table(data.frame(FID = panel$ids, IID = panel$ids, y = y), phenotype_file,
                     row.names = FALSE, quote = FALSE)
  output <- file.path(directory, 'association.stats')
  args <- c(paste0('--bfile=', file.path('data', panel$panel, 'bolt')),
            paste0('--phenoFile=', phenotype_file), '--phenoCol=y', '--lmmInfOnly',
            '--numThreads=1', '--geneticMapFile=data/genetic_map.txt',
            paste0('--statsFile=', output))
  if (isTRUE(use_pc_covariates)) args <- c(args,
    paste0('--covarFile=', file.path('data', panel$panel, 'bolt_pcs.txt')),
    '--qCovarCol=PC{1:3}')
  writeLines(c(executable, args), file.path(directory, 'arguments.txt'))
  status <- system2(executable, vapply(args, shQuote, character(1)),
                    stdout = file.path(directory, 'stdout.log'),
                    stderr = file.path(directory, 'stderr.log'))
  failure <- function(message) list(selected = integer(), valid = FALSE, error = message)
  if (status != 0 || !file.exists(output))
    return(failure(paste('BOLT exit status', status, '; see saved logs')))
  tab <- utils::read.delim(output)
  if (!'P_BOLT_LMM_INF' %in% names(tab)) return(failure('Missing P_BOLT_LMM_INF column'))
  pvalue <- tab$P_BOLT_LMM_INF[match(panel$map$ID, tab$SNP)]
  if (length(pvalue) != panel$p || anyNA(pvalue) || any(!is.finite(pvalue)))
    return(failure('Some tested variants lack finite BOLT p-values'))
  list(selected = which(p.adjust(pvalue, 'BH') <= q), pvalue = pvalue, valid = TRUE, error = '')
}

# Exact-SNP scoring. power = TP / #causal (NA in the null case);
# FDP = FP / max(1, R), so an empty selection has FDP 0 and, in the null
# case, FDP = 1{R > 0}.
score_fit <- function(fit, truth) {
  if (!fit$valid) return(data.frame(valid = FALSE, R = NA, TP = NA, FP = NA,
                                    power = NA, FDP = NA, error = fit$error))
  selected <- fit$selected
  tp <- sum(selected %in% truth)
  fp <- length(selected) - tp
  data.frame(valid = TRUE, R = length(selected), TP = tp, FP = fp,
             power = if (length(truth)) tp / length(truth) else NA_real_,
             FDP = fp / max(1, length(selected)), error = '')
}

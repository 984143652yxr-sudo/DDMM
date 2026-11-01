# Portable adaptation of extended_model/ and family_structure_simulation.R.
# Simulation: y = X beta + u + epsilon; u ~ N(0, eta K), epsilon ~ N(0, (1-eta) I).
source('../R/filters.R')
source('../DDMM group/R/ddmm.R')

cli_settings <- function() {
  args <- commandArgs(TRUE)
  if ('--help' %in% args) {
    cat('Options: --smoke; --n 1000 --p 1200 --reps 20 --generations 20 --seed 1\n')
    quit(status = 0)
  }
  smoke <- '--smoke' %in% args; args <- setdiff(args, '--smoke')
  opt <- list(n = if (smoke) 80L else 1000L, p = if (smoke) 40L else 1200L,
              reps = if (smoke) 1L else 20L, generations = if (smoke) 3L else 20L, seed = 1L)
  if (length(args) %% 2) stop('Options need values; use --help.')
  if (length(args)) for (j in seq(1, length(args), by = 2)) {
    key <- sub('^--', '', args[j]); if (!key %in% names(opt)) stop('Unknown option: ', args[j])
    opt[[key]] <- as.integer(args[j + 1])
  }
  stopifnot(all(is.finite(unlist(opt))), opt$n >= 40, opt$p >= 20, opt$reps >= 1, opt$generations >= 2)
  opt
}

simulate_design <- function(n, p, family = FALSE, generations = 20L) {
  a <- bnpsd::admix_prop_1d_linear(n_ind = n, k_subpops = 5,
          bias_coeff = 0.5, coanc_subpops = 1:5, fst = 0.1)
  obj <- bnpsd::draw_all_admix(a$admix_proportions, a$coanc_subpops, m_loci = 3L * p)
  G <- obj$X
  if (family) {
    pedigree <- simfam::sim_pedigree(n, generations)
    fam <- simfam::prune_fam(pedigree$fam, pedigree$ids[[generations]])
    colnames(G) <- pedigree$ids[[1]]
    G <- simfam::geno_last_gen(G, fam, pedigree$ids)
  }
  raw <- t(G); maf <- pmin(colMeans(raw) / 2, 1 - colMeans(raw) / 2)
  keep <- which(maf >= 0.05 & apply(raw, 2, sd) > 0)
  if (length(keep) < p) stop('Too few polymorphic SNPs; increase the simulated candidate pool.')
  X <- raw[, keep[seq_len(p)], drop = FALSE]
  K <- 2 * popkin::popkin(X, loci_on_cols = TRUE)
  eg <- eigen((K + t(K)) / 2, symmetric = TRUE)
  K <- tcrossprod(sweep(eg$vectors, 2, sqrt(pmax(eg$values, 0)), '*'))
  K <- K / mean(diag(K))
  beta <- numeric(p); causal <- sample.int(p, max(2L, floor(0.1 * p)))
  beta[causal] <- rnorm(length(causal))
  e <- eigen(K, symmetric = TRUE)
  y <- drop(X %*% beta) + sqrt(0.5) * drop(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(n))) + sqrt(0.5) * rnorm(n)
  list(X = X, y = y, K = K, beta = beta, causal = causal)
}

run_comparison <- function(kind) {
  opt <- cli_settings()
  pkgs <- c('bnpsd', 'popkin', 'glmnet', 'ggmix', if (kind == 'family') 'simfam' else 'plmmr')
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Install: ', paste(missing, collapse = ', '))
  out <- file.path('results', kind); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  records <- list()
  for (replicate in seq_len(opt$reps)) {
    set.seed(opt$seed + replicate - 1L)
    d <- simulate_design(opt$n, opt$p, kind == 'family', opt$generations)
    train <- sample.int(opt$n, floor(0.8 * opt$n)); test <- setdiff(seq_len(opt$n), train)
    X <- d$X[train, , drop = FALSE]; y <- d$y[train]; K <- d$K[train, train]
    # Resampling can make a rare SNP constant in the training subset.
    keep <- which(apply(X, 2, sd) > 0); X <- X[, keep, drop = FALSE]
    if (ncol(X) < 5) stop('Too few varying training SNPs.')
    truth <- which(d$beta[keep] != 0)
    add <- function(method, beta, selected) {
      stopifnot(length(beta) == length(keep), all(is.finite(beta)))
      intercept <- mean(y) - sum(colMeans(X) * beta)
      pred <- intercept + drop(d$X[test, keep, drop = FALSE] %*% beta)
      fp <- sum(!selected %in% truth); tp <- sum(selected %in% truth)
      records[[length(records) + 1L]] <<- data.frame(replicate, method, discoveries = length(selected),
        FDP = fp / max(1, length(selected)), power = if (length(truth)) tp / length(truth) else NA_real_,
        MSPE = mean((d$y[test] - pred)^2), L1 = sum(abs(beta - d$beta[keep])))
    }
    gg <- ggmix::gic(ggmix::ggmix(x = X, y = y, kinship = K))
    b <- as.numeric(coef(gg))[-1]
    add(if (kind == 'family') 'LMM-Lasso (ggmix)' else 'ggmix', b, which(b != 0))
    eta <- as.numeric(gg$ggmix_fit$eta[, gg$lambda.min.name]); eta <- min(0.99, max(0, eta))
    C <- eta * K + (1 - eta) * diag(nrow(X))
    if (kind == 'admixture') {
      design <- plmmr::create_design(X = X, y = y)
      cv <- plmmr::cv_plmm(design, nfolds = 5, seed = opt$seed + replicate)
      b <- as.numeric(coef(cv))[-1]
      add('plmmr', b, which(b != 0))
      f <- spectral_filter(X, y, 'PCA', k = 3)
      cv <- glmnet::cv.glmnet(f$X, f$y, intercept = FALSE, standardize = FALSE, nfolds = 5)
      b <- as.numeric(coef(cv, s = 'lambda.min'))[-1] / f$scale
      add('PCA-Lasso', b, which(b != 0))
    }
    for (method in if (kind == 'family') c('LMM', 'TRIM') else c('PCA', 'LMM', 'TRIM')) {
      f <- spectral_filter(X, y, method, k = 3, K = K, eta = eta, noise_cov = C)
      fit <- variant_tests(f$y, f$X, new_nodewise_cache(f$X), noise_cov = f$noise_cov)
      sel <- which(p.adjust(fit$p, 'BH') <= 0.1)
      # Prediction uses the selected debiased coefficients, on the raw SNP scale.
      b <- numeric(length(keep))
      residual <- f$y - drop(f$X %*% fit$beta)
      if (length(sel)) {
        cache <- new_nodewise_cache(f$X); ensure_nodewise(cache, sel)
        for (j in sel) {
          nw <- cache$fits[[j]]; z <- f$X[, j] - drop(f$X[, -j, drop = FALSE] %*% nw$g)
          b[j] <- (fit$beta[j] + sum(z * residual) / nw$zx) / f$scale[j]
        }
      }
      add(paste0('DDMM-', method), b, sel)
    }
    write.csv(do.call(rbind, records), file.path(out, 'replicates.csv'), row.names = FALSE)
    cat(kind, 'replicate', replicate, '/', opt$reps, '\n')
  }
  saveRDS(opt, file.path(out, 'settings.rds'))
  cat('Run: Rscript plot.R', kind, '\n')
}

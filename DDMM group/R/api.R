# Programmatic interface: the kept (corrected) DDMM group tests on any
# genotype matrix and any partition, outside the two toy experiments.
#
#   source('R/ddmm.R'); source('R/api.R')
#   out <- ddmm_group_test(y, X, groups)                              # effective rank (default)
#   out <- ddmm_group_test(y, X, groups, rank = 3)                    # fixed K = 3
#   out <- ddmm_group_test(y, X, groups, rank = c('effective', 3, 'full'))
#   out$groups        one row per group: size, rank used, p-values, BH selections
#   out$variants      per-variant debiased Lasso statistic and p-value
#
# Arguments
#   y        phenotype, length n (centred inside)
#   X        n x p genotype (dosage) matrix, columns in genomic order;
#            missing values are imputed by the column mean, then standardised
#   groups   length-p group labels (any type); a group should be a contiguous
#            run of columns, as in every LD-based partition
#   rank     'effective' (participation ratio of the group's LD eigenvalues),
#            a positive integer K, 'full', or a vector of these; one chi-squared
#            track per entry
#   chr      optional length-p chromosome labels: group-nodewise neighbours are
#            taken within the chromosome only
#   q        Benjamini-Hochberg level within the partition
#   seed     seed set before the cross-validated pilot Lasso (folds are random)

ddmm_group_test <- function(y, X, groups, rank = 'effective', chr = NULL, q = 0.1, seed = 1) {
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(y) == n, length(groups) == p, !anyNA(y), !anyNA(groups))
  ranks <- parse_ranks(paste(rank, collapse = ','))

  missing <- which(is.na(X), arr.ind = TRUE)
  if (nrow(missing)) X[missing] <- colMeans(X, na.rm = TRUE)[missing[, 'col']]
  Xs <- scale(X); Xs[is.na(Xs)] <- 0                       # monomorphic columns -> 0
  Xc <- scale(Xs); Xc[is.na(Xc)] <- 0                      # as load_genotypes(), so that
  geno <- list(Xc = Xc, xcen = scale(Xs, scale = FALSE),    # run.R is reproduced exactly
               gram = crossprod(Xc) / n,
               chr = if (is.null(chr)) rep(1L, p) else chr, n = n, p = p)

  labels <- unique(groups)
  assignment <- match(groups, labels)                      # 1..G in order of appearance
  sizes <- tabulate(assignment, length(labels))

  set.seed(seed)
  fit <- variant_tests(y - mean(y), geno$xcen, new_nodewise_cache(geno$xcen))
  residual <- (y - mean(y)) - as.numeric(geno$Xc %*% fit$beta)

  p_group <- cbind(simes = simes(fit$p, assignment, sizes))
  used <- NULL
  if (any(sizes > 1)) {
    gt <- group_tests(group_operators(geno, assignment), fit$beta, residual, fit$sigma^2, ranks)
    p_group <- cbind(p_group, gt$p); used <- gt$rank
  } else {                                                 # single-variant partition
    p_group <- cbind(p_group, matrix(fit$p, p, length(ranks)))
    used <- matrix(1L, p, length(ranks))
  }
  tracks <- c('simes', track_id(ranks))
  colnames(p_group) <- tracks

  first <- vapply(split(seq_len(p), assignment), min, 1L)
  last <- vapply(split(seq_len(p), assignment), max, 1L)
  out <- data.frame(group = labels, size = sizes, first = first, last = last)
  for (i in seq_along(ranks)) out[[paste0('K_', ranks[i])]] <- used[, i]
  for (t in tracks) out[[paste0('p_', t)]] <- p_group[, t]
  for (t in tracks) {
    sel <- logical(length(labels)); sel[bh(p_group[, t], q)] <- TRUE
    out[[paste0('selected_', t)]] <- sel
  }
  list(groups = out,
       variants = data.frame(column = seq_len(p), group = groups, tested = seq_len(p) %in% fit$S,
                             statistic = fit$statistic, p = fit$p),
       sigma = fit$sigma, q = q, ranks = ranks)
}

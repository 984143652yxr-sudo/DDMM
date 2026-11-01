# Genotype-only description of each cohort of the group benchmark, per resolution of the
# adjclust partitions closed under r^2 >= 0.8: groups, mean size, median width, effective
# rank (participation ratio of the group's LD spectrum), LD variance kept by the effective
# rank and by K = 3, groups that contain a causal variant, the realised FDP of the trivial
# rule that declares every group a discovery, loci whose two causal variants fall in
# different groups, the share of variant pairs with r^2 >= 0.8 split between groups (0 by
# construction) and the largest r^2 between a causal variant and a variant of a group with
# no causal variant. Also LD around the causal variants.
#   Rscript describe_design.R EUR MIX ALL  -> results/design_summary.csv, results/causal_ld.csv
suppressMessages(library(data.table))
script_dir <- function() { a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(gsub('~+~', ' ', sub('^--file=', '', a[1]), fixed = TRUE))) else 'code' }
B <- Sys.getenv('BENCH', normalizePath(file.path(script_dir(), '..')))      # this folder
PKG <- Sys.getenv('DDMM_PKG', normalizePath(file.path(B, '..', '..')))     # ddmm_group package root
SCHEMES <- commandArgs(TRUE)
source(file.path(PKG, 'R', 'ddmm.R'))
out <- list(); ld <- list()
for (coh in SCHEMES) {
  dir <- file.path(B, 'data', coh)
  design0 <- readRDS(file.path(dir, 'design.rds')); R2 <- cor(design0$X * 1.0)^2; diag(R2) <- 0; R2[outer(design0$map$CHR, design0$map$CHR, '!=')] <- 0
  hi <- which(R2 >= 0.8 & upper.tri(R2), arr.ind = TRUE)
  design <- readRDS(file.path(dir, 'design.rds'))
  X <- scale(design$X * 1.0); X[is.na(X)] <- 0; n <- nrow(X)
  part <- rbindlist(lapply(c(21, 22), function(ch) fread(file.path(dir, 'partitions', sprintf('part_chr%d.txt', ch)))))
  part <- part[match(design$map$SNP, part$SNP)]
  causal_idx <- match(design$causal$SNP, design$map$SNP)
  for (k in 1:7) {
    col <- sprintf('res_%d', k); g <- integer(ncol(X)); off <- 0L
    for (ch in c(21, 22)) { i <- which(design$map$CHR == ch); loc <- as.integer(factor(part[[col]][i])); g[i] <- loc + off; off <- off + max(loc) }
    sizes <- tabulate(g); m <- length(sizes)
    stats <- t(vapply(seq_len(m), function(j) {
      mem <- which(g == j)
      ev <- if (length(mem) == 1) 1 else pmax(eigen(crossprod(X[, mem, drop = FALSE]) / n, symmetric = TRUE, only.values = TRUE)$values, 0)
      K <- effective_rank(ev)
      c(K = K, kept_eff = sum(ev[seq_len(K)]) / sum(ev), kept_3 = sum(ev[seq_len(min(3, length(ev)))]) / sum(ev))
    }, numeric(3)))
    has_causal <- seq_len(m) %in% unique(g[causal_idx])
    out[[length(out) + 1]] <- data.table(cohort = coh, resolution = sprintf('res%d', 7 - k),
      percent = rev(c(100, 50, 20, 10, 5, 2, 1))[k], groups = m, mean_size = mean(sizes),
      median_width_kb = median(tapply(design$map$BP, g, function(b) (max(b) - min(b)) / 1e3)),
      mean_eff_rank = mean(stats[, 'K']), kept_eff = mean(stats[, 'kept_eff']), kept_3 = mean(stats[, 'kept_3']),
      groups_with_causal = sum(has_causal), trivial_FDP = mean(!has_causal),
      loci_split = sum(tapply(g[causal_idx], design$causal$Locus, function(v) length(unique(v)) > 1)),
      r2_pairs_split = mean(g[hi[, 1]] != g[hi[, 2]]),
      max_r2_causal_to_null_group = if (all(g %in% unique(g[causal_idx]))) NA_real_ else max(R2[causal_idx, !(g %in% unique(g[causal_idx])), drop = FALSE]))
  }
  # LD of each causal variant with its neighbours (|r| with variants within 100 kb)
  R <- NULL
  for (j in causal_idx) {
    nb <- setdiff(which(design$map$CHR == design$map$CHR[j] & abs(design$map$BP - design$map$BP[j]) <= 1e5), j)
    r <- abs(cor(X[, j], X[, nb]))
    ld[[length(ld) + 1]] <- data.table(cohort = coh, SNP = design$map$SNP[j], maf = min(mean(design$X[, j]) / 2, 1 - mean(design$X[, j]) / 2),
      neighbours_100kb = length(nb), max_r = max(r), n_r2_above_0.5 = sum(r^2 > 0.5))
  }
}
dir.create(file.path(B, 'results'), showWarnings = FALSE)
d <- rbindlist(out); fwrite(d, file.path(B, 'results', 'design_summary.csv'))
l <- rbindlist(ld); fwrite(l, file.path(B, 'results', 'causal_ld.csv'))
print(d[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
print(l[, .(median_max_r = median(max_r), mean_n_r2_0.5 = mean(n_r2_above_0.5), median_maf = median(maf)), by = cohort])

## ============================================================
##  make_adjclust_partitions.R — nested LD-block partitions of one cohort, built with the
##  KnockoffZoom grouping rule and closed under strong LD.
##
##  Step 1 (KnockoffZoom, utils/cluster.R). Within each chromosome, the dissimilarity
##    1 - |r_jk| between genotype dosages is clustered by adjacency-constrained
##    hierarchical clustering with Ward-type linkage (adjclust::adjClust), so every
##    cluster is a run of adjacent variants.
##  Step 2 (KnockoffZoom, utils/partition.R). The dendrogram is cut into
##    round(0.01 * res * p_chr) groups for res = 100, 50, 20, 10, 5, 2, 1 (per cent of
##    the variants), the seven resolutions of the KnockoffZoom UK Biobank analysis.
##    Cutting by the number of groups follows the merge order, so the partitions are
##    nested even when adjClust reports non-monotone merge heights.
##  Step 3 (LD closure, as make_ld_partitions.R of the ld_aware benchmark). Variants with
##    r^2 >= THR (0.8, the tagging threshold of de Bakker et al. 2005) are never split.
##    At every resolution the groups are joined with the connected components of the
##    graph r^2 >= THR and closed to intervals of adjacent variants, which keeps the
##    partitions contiguous (required by snpknock2) and nested.
##  Output: partitions/part_chr{21,22}.txt in the KnockoffGWAS format (SNP res_7 ... res_1,
##  res_7 = finest, labels increasing along the chromosome) and partition_summary.csv.
##    Rscript make_adjclust_partitions.R <cohort dir> [r2 threshold = 0.8]
## ============================================================
suppressMessages({ library(data.table); library(adjclust) })
args <- commandArgs(TRUE); dir <- args[1]; THR <- if (length(args) > 1) as.numeric(args[2]) else 0.8
RES_PCT <- c(100, 50, 20, 10, 5, 2, 1)                # res0 ... res6, columns res_7 ... res_1
dir.create(file.path(dir, "partitions"), showWarnings = FALSE)
design <- readRDS(file.path(dir, "design.rds")); X <- design$X * 1.0; map <- design$map

uf_join <- function(a, b) {                       # finest partition coarser than a and b
  n <- length(a); parent <- seq_len(n)
  find <- function(i) { while (parent[i] != i) { parent[i] <<- parent[parent[i]]; i <- parent[i] }; i }
  unite <- function(idx) if (length(idx) > 1) { r0 <- find(idx[1]); for (i in idx[-1]) { ri <- find(i); if (ri != r0) parent[ri] <<- r0 } }
  for (g in split(seq_len(n), a)) unite(g); for (g in split(seq_len(n), b)) unite(g)
  vapply(seq_len(n), find, 1L)
}
interval_closure <- function(lab) {               # merge groups whose index ranges overlap
  repeat { lo <- tapply(seq_along(lab), lab, min); hi <- tapply(seq_along(lab), lab, max)
    o <- order(lo); lo <- lo[o]; hi <- hi[o]; grp <- integer(length(lo)); g <- 1L; cur <- hi[1]; grp[1] <- 1L
    for (i in seq_along(lo)[-1]) { if (lo[i] <= cur) grp[i] <- g else { g <- g + 1L; grp[i] <- g }; cur <- max(cur, hi[i]) }
    new <- grp[match(lab, as.integer(names(lo)))]
    if (max(new) == length(unique(lab))) return(new); lab <- new }
}
ld_clusters <- function(R2, thr) {                # connected components of r^2 >= thr
  A <- R2 >= thr; diag(A) <- FALSE; p <- ncol(R2); m <- integer(p); k <- 0L
  for (s in seq_len(p)) if (!m[s]) { k <- k + 1L; q <- s; m[s] <- k
    while (length(q)) { nb <- which(A[q[1], ] & !m); m[nb] <- k; q <- c(q[-1], nb) } }
  m
}
relabel <- function(g) as.integer(factor(g, levels = unique(g)))
split_pairs <- function(g, R2, thr) { A <- R2 >= thr; A[lower.tri(A, diag = TRUE)] <- FALSE; w <- which(A, arr.ind = TRUE); sum(g[w[, 1]] != g[w[, 2]]) }

summ <- list()
for (ch in c(21, 22)) {
  idx <- which(map$CHR == ch); stopifnot(all(diff(map$BP[idx]) > 0))
  R <- suppressWarnings(cor(X[, idx])); R[is.na(R)] <- 0; R2 <- R^2
  fit <- suppressMessages(adjClust(as.dist(1 - abs(R)), type = "dissimilarity"))   # Step 1
  cl <- ld_clusters(R2, THR)
  out <- data.table(SNP = map$SNP[idx]); prev <- NULL
  for (i in seq_along(RES_PCT)) {
    k <- round(0.01 * RES_PCT[i] * length(idx))
    g0 <- if (k >= length(idx)) seq_along(idx) else relabel(stats::cutree(fit, k = k))      # Step 2
    stopifnot(all(diff(g0) %in% c(0L, 1L)))                                               # contiguous
    g <- relabel(interval_closure(uf_join(g0, cl)))                                       # Step 3
    col <- sprintf("res_%d", 8 - i); out[[col]] <- g
    w <- tapply(map$BP[idx], g, function(b) (max(b) - min(b)) / 1e3)
    summ[[length(summ) + 1]] <- data.table(chr = ch, resolution = sprintf("res%d", i - 1), percent = RES_PCT[i],
      groups_adjclust = max(g0), pairs_split_adjclust = split_pairs(g0, R2, THR), groups = max(g),
      pairs_split = split_pairs(g, R2, THR), singletons = sum(tabulate(g) == 1), max_size = max(tabulate(g)),
      median_width_kb = median(w), mean_size = length(idx) / max(g))
  }
  for (k in 7:2) stopifnot(all(tapply(out[[sprintf("res_%d", k - 1)]], out[[sprintf("res_%d", k)]], function(z) length(unique(z))) == 1))
  fwrite(out, file.path(dir, "partitions", sprintf("part_chr%d.txt", ch)), sep = " ")
}
S <- rbindlist(summ)
tot <- S[, .(percent = percent[1], groups_adjclust = sum(groups_adjclust), pairs_split_adjclust = sum(pairs_split_adjclust),
             groups = sum(groups), pairs_split = sum(pairs_split), singletons = sum(singletons), max_size = max(max_size),
             mean_size = 2000 / sum(groups)), by = resolution]
## median width over the groups of both chromosomes
widths <- rbindlist(lapply(c(21, 22), function(ch) { P <- fread(file.path(dir, "partitions", sprintf("part_chr%d.txt", ch)))
  rbindlist(lapply(7:1, function(k) { bp <- map$BP[match(P$SNP, map$SNP)]
    data.table(resolution = sprintf("res%d", 7 - k), width_kb = as.numeric(tapply(bp, P[[sprintf("res_%d", k)]], function(b) (max(b) - min(b)) / 1e3))) })) }))
tot <- merge(tot, widths[, .(median_width_kb = median(width_kb), mean_width_kb = mean(width_kb), max_width_kb = max(width_kb)), by = resolution], by = "resolution")
fwrite(S, file.path(dir, "partitions", "partition_summary_by_chr.csv")); fwrite(tot, file.path(dir, "partitions", "partition_summary.csv"))
cat(sprintf("%s: r^2 >= %.2f clusters closed, %d variants\n", basename(dir), THR, ncol(X))); print(tot)

# Partitions of the 2000 variants into groups, one per resolution, ordered from
# coarse to fine. Both schemes depend on the genotypes only.
#
#   published  (KnockoffGWAS, the setting the package uses) the nested
#              partitions shipped with KnockoffGWAS: DDMM tests the groups the
#              published KnockoffGWAS results refer to
#   adjclust   (KnockoffZoom, archived in legacy/knockoffzoom) contiguity-
#              constrained clustering of the |r| LD matrix per chromosome, cut
#              to KnockoffZoom's group counts; this matches the counts but only
#              part of KnockoffZoom's group boundaries, which its shipped files
#              do not fully record
#
# A resolution whose groups are all single variants ('SNP', 'res0') is tested by
# the per-variant p-values for every track.

relabel <- function(g) as.integer(factor(g, levels = unique(g)))

# adjclust tree of one chromosome; heights are forced increasing because
# adjClust can return non-monotone merge heights.
ld_tree <- function(x) {
  if (ncol(x) < 2) return(NULL)
  s <- abs(suppressWarnings(stats::cor(x)))
  s[is.na(s)] <- 0
  fit <- suppressWarnings(suppressMessages(adjclust::adjClust(s)))
  structure(list(merge = fit$merge,
                 height = cummax(pmax(fit$height, 0)) + cumsum(rep(1e-9, ncol(x) - 1)),
                 order = seq_len(ncol(x)), labels = NULL, method = 'adjClust'),
            class = 'hclust')
}

build_partitions <- function(cfg, geno) {
  chrs <- sort(unique(geno$chr))
  trees <- NULL
  if (cfg$partition == 'adjclust') {
    trees <- lapply(chrs, function(ch) ld_tree(geno$X[, geno$chr == ch, drop = FALSE]))
    names(trees) <- as.character(chrs)
  } else if (cfg$partition == 'published') {
    part <- data.table::rbindlist(lapply(cfg$chromosomes, function(ch) {
      tab <- data.table::fread(file.path(cfg$data_dir, sprintf(cfg$partition_files, ch)))
      tab$CHR <- ch
      tab
    }))
    part <- part[match(geno$snp, part$SNP), ]
    stopifnot(!anyNA(part$SNP))
  } else stop('Unknown partition scheme ', cfg$partition)

  resolutions <- lapply(cfg$resolutions, function(res) {
    if (isTRUE(res$single_variant)) {
      g <- seq_len(geno$p)
    } else {
      g <- integer(geno$p); offset <- 0L
      for (ch in chrs) {
        idx <- which(geno$chr == ch)
        local <- if (cfg$partition == 'adjclust') {
          k <- min(res$groups[[as.character(ch)]], length(idx))
          tree <- trees[[as.character(ch)]]
          if (is.null(tree) || k <= 1) rep(1L, length(idx)) else relabel(stats::cutree(tree, k = k))
        } else as.integer(factor(part[[res$column]][idx]))
        g[idx] <- local + offset
        offset <- offset + max(local)
      }
      if (cfg$partition == 'adjclust') g <- relabel(g)
    }
    sizes <- as.integer(table(factor(g, levels = seq_len(max(g)))))
    list(name = res$name, assignment = g, sizes = sizes,
         single_variant = all(sizes == 1),
         width = if (all(sizes == 1)) 'single SNP' else
           sprintf('%.3f Mb', stats::median(tapply(geno$pos, g, function(z) diff(range(z)))) / 1e6))
  })
  names(resolutions) <- vapply(resolutions, `[[`, '', 'name')
  list(resolutions = resolutions, trees = trees)
}

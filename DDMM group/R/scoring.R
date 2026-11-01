# Scoring discoveries against the 20 causal loci.
#
# A discovery is the genomic interval spanned by its group. It is a true
# discovery when the interval contains a causal variant; a locus is found when
# some discovery contains one of its causal variants.
#   power = loci found / 20      FDP = false discoveries / max(1, discoveries)
# In the null case there is no causal variant, so FDP = 1{any discovery}.

group_intervals <- function(selected, assignment, geno) {
  if (!length(selected)) return(data.frame(chr = integer(), start = numeric(), end = numeric()))
  do.call(rbind, lapply(selected, function(k) {
    m <- which(assignment == k)
    data.frame(chr = geno$chr[m[1]], start = min(geno$pos[m]), end = max(geno$pos[m]))
  }))
}

score_intervals <- function(iv, causal, n_loci) {
  if (!nrow(iv)) return(data.frame(discoveries = 0L, TP = 0L, FDP = 0, loci_found = 0L,
                                   power = if (is.null(causal)) NA_real_ else 0, width_kb = NA_real_))
  if (is.null(causal)) return(data.frame(discoveries = nrow(iv), TP = 0L, FDP = 1, loci_found = 0L,
                                         power = NA_real_, width_kb = mean(iv$end - iv$start) / 1e3))
  hit <- matrix(vapply(seq_len(n_loci), function(l) {
    s <- causal[causal$Locus == l, ]
    vapply(seq_len(nrow(iv)), function(i)
      any(s$CHR == iv$chr[i] & s$BP >= iv$start[i] & s$BP <= iv$end[i]), logical(1))
  }, logical(nrow(iv))), nrow = nrow(iv))
  tp <- sum(rowSums(hit) > 0)
  loci <- sum(colSums(hit) > 0)
  data.frame(discoveries = nrow(iv), TP = tp, FDP = (nrow(iv) - tp) / nrow(iv),
             loci_found = loci, power = loci / n_loci, width_kb = mean(iv$end - iv$start) / 1e3)
}

# Published knockoff discoveries for the shipped phenotype, as intervals.
comparator_intervals <- function(cfg, layout) {
  out <- lapply(layout$resolutions, function(name) {
    file <- cfg$comparator$overrides[[name]] %||% sprintf(cfg$comparator$files, name)
    d <- data.table::fread(file.path(cfg$data_dir, file))
    if (!nrow(d)) data.frame(chr = integer(), start = numeric(), end = numeric())
    else data.frame(chr = d$CHR, start = d$BP.min, end = d$BP.max)
  })
  names(out) <- layout$resolutions
  out
}

# The same intervals scored like the DDMM discoveries.
comparator_scores <- function(cfg, layout) {
  iv <- comparator_intervals(cfg, layout)
  do.call(rbind, lapply(layout$resolutions, function(name)
    cbind(data.frame(case = 'shipped', replicate = 0L, resolution = name, track = 'knockoff',
                     groups = layout$groups[[name]], mean_rank = NA_real_),
          score_intervals(iv[[name]], layout$causal, layout$n_loci))))
}

# Genotype-only description of an experiment, saved once per run directory so
# reports never need the genotypes.
make_layout <- function(ctx) {
  res <- ctx$parts$resolutions
  list(experiment = ctx$cfg$experiment, chr = ctx$geno$chr, pos = ctx$geno$pos,
       causal = as.data.frame(ctx$truth$causal), n_loci = ctx$truth$n_loci,
       resolutions = names(res),
       groups = lapply(res, function(r) length(r$sizes)),
       widths = vapply(res, `[[`, '', 'width'),
       assignments = lapply(res, `[[`, 'assignment'),
       trees = ctx$parts$trees,
       row_labels = rev(ifelse(vapply(res, `[[`, NA, 'single_variant'), 'single SNP',
                               sprintf('%s (%s)', names(res), vapply(res, `[[`, '', 'width')))))
}

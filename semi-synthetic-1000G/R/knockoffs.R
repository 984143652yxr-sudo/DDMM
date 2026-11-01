# HMM knockoff copies (snpknock2, Sesia, Sabatti & Candes).
#
# One copy per panel x global replicate, generated from the genotypes only,
# before any phenotype exists; the signal and null cases of that replicate
# share it. It is stored in the replicate directory (knockoff.rds, raw
# snpknock2 output in hmm/) and re-read when the replicate is resumed.
# snpknock2 is deterministic given the seed: the same replicate always gets
# the same copy (checked to be bit-identical on macOS and Linux builds).

SNPKNOCK2_OPTIONS <- c('--K', '20', '--cluster_size_min', '100', '--cluster_size_max', '500',
                       '--hmm-rho', '1', '--hmm-lambda', '0.001', '--windows', '0')

get_knockoff <- function(panel, replicate, replicate_dir, cfg) {
  seed <- knockoff_seed(cfg, panel$panel, replicate)
  path <- file.path(replicate_dir, 'knockoff.rds')
  if (file.exists(path)) {
    saved <- readRDS(path)
    stopifnot(saved$seed == seed)
  } else if (file.exists(file.path(replicate_dir, 'knockoff_path.txt'))) {
    # Legacy kinship600 layout: pointer to a shared copy under data/MIX600/.
    saved <- readRDS(readLines(file.path(replicate_dir, 'knockoff_path.txt'))[1])
  } else {
    saved <- generate_knockoff(panel, replicate, replicate_dir, seed)
  }
  stopifnot(saved$replicate == replicate, identical(saved$panel, panel$panel),
            identical(saved$ids, panel$ids), identical(saved$snps, panel$map$ID),
            identical(dim(saved$X), dim(panel$D)))
  saved$X
}

generate_knockoff <- function(panel, replicate, replicate_dir, seed) {
  executable <- snpknock2_path()
  if (!nzchar(executable)) stop('Set SNPKNOCK_BIN to the snpknock2 executable.')
  work <- file.path(replicate_dir, 'hmm')
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(normalizePath(work), 'hmm')
  input <- file.path('data', panel$panel)
  arguments <- c('--haps', file.path(input, 'hap'), '--map', 'data/genetic_map_chr1.txt',
                 '--part', file.path(input, 'partitions.txt'),
                 '--pc', file.path(input, 'hmm_pcs.txt'),
                 SNPKNOCK2_OPTIONS, '--n_threads', '2', '--seed', as.character(seed),
                 '--compute-references', '--estimate-hmm', '--generate-knockoffs',
                 '--out', prefix)
  writeLines(c(executable, arguments), file.path(work, 'arguments.txt'))
  status <- system2(executable, vapply(arguments, shQuote, ''),
                    stdout = file.path(work, 'snpknock2.log'),
                    stderr = file.path(work, 'snpknock2_error.log'))
  if (status != 0) stop('snpknock2 failed; see ', file.path(work, 'snpknock2.log'))

  # Read original + knockoff dosages and align alleles with panel$D.
  stem <- paste0(prefix, '_res0')
  plink <- snpStats::read.plink(paste0(stem, '.bed'), paste0(stem, '.bim'), paste0(stem, '.fam'))
  dosage <- methods::as(plink$genotypes, 'numeric')
  stopifnot(!anyDuplicated(plink$fam$member))
  dosage <- dosage[match(panel$ids, plink$fam$member), , drop = FALSE]
  original <- dosage[, match(panel$map$ID, colnames(dosage)), drop = FALSE]
  copy <- dosage[, match(paste0(panel$map$ID, '.k'), colnames(dosage)), drop = FALSE]
  stopifnot(!anyNA(original), !anyNA(copy), identical(dim(copy), dim(panel$D)))
  for (j in seq_len(panel$p)) {
    if (all(original[, j] == panel$D[, j])) next
    if (!all(2 - original[, j] == panel$D[, j])) stop('Allele alignment failed: SNP ', j)
    copy[, j] <- 2 - copy[, j]
  }
  saved <- list(X = copy, seed = seed, panel = panel$panel, replicate = replicate,
                ids = panel$ids, snps = panel$map$ID)
  save_rds_atomic(saved, file.path(replicate_dir, 'knockoff.rds'))
  saved
}

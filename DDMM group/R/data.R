# Genotypes, causal annotation and phenotypes of the KnockoffZoom / KnockoffGWAS
# toy data (Sesia et al. 2020, 2021). Both repositories ship byte-identical
# genotypes (n = 1000, p = 2000 on chromosomes 21 and 22) and phenotype, with
# 40 causal variants in 20 loci.

load_genotypes <- function(cfg) {
  plink <- lapply(cfg$chromosomes, function(ch) snpStats::read.plink(
    file.path(cfg$data_dir, 'data/genotypes', sprintf('example_chr%d.bed', ch))))
  dosage <- do.call(cbind, lapply(plink, function(x) methods::as(x$genotypes, 'numeric')))
  map <- do.call(rbind, lapply(plink, `[[`, 'map'))
  missing <- which(is.na(dosage), arr.ind = TRUE)
  if (nrow(missing)) dosage[missing] <- colMeans(dosage, na.rm = TRUE)[missing[, 'col']]

  X <- scale(dosage)                       # standardised genotypes
  Xc <- scale(X); Xc[is.na(Xc)] <- 0       # design of the group tests
  list(X = X, Xc = Xc,
       xcen = scale(X, scale = FALSE),     # design of the per-variant test
       gram = crossprod(Xc) / nrow(X),
       chr = map$chromosome, pos = map$position, snp = rownames(map),
       n = nrow(X), p = ncol(X))
}

load_truth <- function(cfg, geno) {
  causal <- data.table::fread(file.path(cfg$data_dir, 'data/phenotypes/causal.txt'))
  index <- match(causal$SNP, geno$snp)
  stopifnot(!anyNA(index))
  beta <- numeric(geno$p)
  beta[index] <- cfg$phenotype$effect_scale * causal$Sign * causal$Scale
  list(causal = causal, beta = beta, n_loci = max(causal$Locus))
}

# Phenotype of one case, centred.
#   shipped : the phenotype published with the data
#   signal  : y = X beta + e,  e ~ N(0, noise_sd^2), with the effects recovered
#             from the shipped phenotype (h2 about 0.71); a new e per replicate
#   null    : y = e, the same e as the signal case of that replicate
# The random state is left where the cross-validated pilot Lasso continues
# from, so a replicate is reproduced exactly by its ID.
case_phenotype <- function(case, replicate, cfg, geno, truth) {
  if (case == 'shipped') {
    y <- data.table::fread(file.path(cfg$data_dir, 'data/phenotypes/phenotypes.tab'))$y
    set.seed(cfg$phenotype$shipped_seed)
    return(as.numeric(scale(y, scale = FALSE)))
  }
  set.seed(cfg$phenotype$seed_base + replicate)
  signal <- if (case == 'signal') geno$X %*% truth$beta else 0
  y <- as.numeric(signal + stats::rnorm(geno$n, 0, cfg$phenotype$noise_sd))
  y - mean(y)
}

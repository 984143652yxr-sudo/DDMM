# Phenotype model. All four methods are fitted to the SAME simulated phenotype
# of a given panel x replicate x amplitude case.
#
#   y = a * X beta0 + u + epsilon,   epsilon ~ N(0, sigma_epsilon^2 I)
#   u ~ N(0, sigma_g^2 K) only when kinship_effect is true (kinship600),
#       drawn as U diag(sqrt(values)) z with K = U diag(values) U^T.
#
# beta0 = panel$betas$SNP holds the 20 causal SNP effects (fixed per panel);
# genotypes X are fixed; epsilon and u are fresh for every replicate and
# amplitude. The seed rule is documented in R/config.R.

simulate_phenotype <- function(panel, amplitude, replicate, cfg) {
  seed <- phenotype_seed(cfg, panel$panel, replicate, amplitude)
  set.seed(as.integer(seed), kind = "L'Ecuyer-CMRG")
  epsilon <- cfg$sigma_epsilon * rnorm(panel$n)
  u <- if (isTRUE(cfg$kinship_effect))
    cfg$sigma_g * drop(panel$U %*% (sqrt(panel$values) * rnorm(ncol(panel$U))))
  else rep(0, panel$n)
  beta <- amplitude * panel$betas$SNP
  y <- drop(panel$X %*% beta) + u + epsilon
  list(y = y, beta = beta, u = u, epsilon = epsilon, seed = seed)
}

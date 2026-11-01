# Data sources

- **Toy example:** KnockoffGWAS example genotypes, phenotypes, partitions, and
  published discovery lists (Sesia et al., 2021, *PNAS* 118, e2105841118).
- **1000G benchmark:** 2,000 SNPs from chromosomes 21 and 22 of the 1000 Genomes
  Phase 3 panel. EUR has 503 individuals; MIX has 504 (168 each EUR/AFR/EAS);
  ALL has 2,504. Outcomes are simulated: 40 causal variants in 20 loci,
  amplitudes 8 and 4, and a global null, with 20 replicates per setting.
- **Saved selections:** DDMM, KnockoffGWAS, and SuSiE outputs retained as inputs
  for rescoring. SuSiE uses L = 40, 95%/90% credible sets, and PIP > 0.9
  (Wang et al., 2020, *JRSS B* 82, 1273–1300).

The compact `design.rds` files contain genotypes, variant/sample metadata,
and simulated causal truth. Raw VCFs, phased haplotypes, and knockoff matrices
are excluded. Upstream data and software retain their own attribution and terms.

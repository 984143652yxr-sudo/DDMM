# Sources and release provenance

The release preserves the supplied `ddmm_group` analysis code and cached
results. It includes two distinct data sources:

- **KnockoffGWAS toy example:** small example genotypes, published partitions,
  phenotype/causal files, and comparator discovery/statistic files under
  `data/knockoffgwas/`. DDMM fits are under `results/knockoffgwas/`.
- **1000 Genomes Phase 3 subsets:** EUR, MIX, and ALL panels of 2,000 SNPs
  from chromosomes 21 and 22. `design.rds` retains genotype dosages and public
  sample metadata; phenotypes and causal effects were simulated by the supplied
  benchmark scripts. `raw_preparation/` documents this construction.

`input_checksums.csv` records each bundled input's relative path, byte count,
and MD5 checksum. The reproduction entry point checks those hashes before
analysis. Input copies retain their original content; filesystem symbolic
links have been materialized for portability.

Packaging changes: a portable reproduction entry point, completeness and
reference-result checks, directory names containing spaces supported by R
script-path handling, concise documentation, and separation of reference
outputs from regenerated outputs. Statistical fitting and scoring formulas
are unchanged. Machine-specific jobs, caches, logs, and archived experiments
are omitted.

## References used by the supplied implementation

- Sesia, M., Bates, S., Candès, E., Marchini, J. and Sabatti, C. (2021).
  False discovery rate control in genome-wide association studies with
  population structure. *PNAS* 118(40), e2105841118. (KnockoffGWAS.)
- Sesia, M., Katsevich, E., Bates, S., Candès, E. and Sabatti, C. (2020).
  Multi-resolution localization of causal variants across the genome.
  *Nature Communications* 11, 1093. (Multi-resolution grouping.)
- Wang, G., Sarkar, A., Carbonetto, P. and Stephens, M. (2020).
  A simple new approach to variable selection in regression, with application
  to genetic fine mapping. *JRSS B* 82, 1273–1300. (SuSiE.)
- Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate.
  *JRSS B* 57(1), 289–300.
- van de Geer, S., Bühlmann, P., Ritov, Y. and Dezeure, R. (2014).
  On asymptotically optimal confidence regions and tests for
  high-dimensional models. *Annals of Statistics* 42(3), 1166–1202.

Upstream data and software retain their own attribution and terms; no new
license is assigned to third-party material by this packaging step.

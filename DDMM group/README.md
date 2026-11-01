# DDMM group

Run terminal commands from `DDMM group/`.

## 1. 1000G: simulate phenotypes and compare group tests

**Real genotypes → simulated phenotypes → fitted selections → power and FDR.**

The included `benchmark_1000g/group_supplement/data/<panel>/design.rds`
contains the genotype matrix `X` (individuals × SNPs), variant metadata,
and 40 causal SNPs in 20 loci. Each panel has 2,000 SNPs from chromosomes 21–22:
EUR (503 individuals), MIX (504), and ALL (2,504).

For fixed standardized genotypes, simulate

$$y=z\theta+X_{\mathrm{causal}}\beta+\epsilon,\qquad
\beta_j=\frac{a\,s_j b_j}{\sqrt n},\qquad
\epsilon\sim N(0,\sigma^2 I).$$

Here `z` is standardized sex, $\theta=10/\sqrt n$, and $s_j,b_j$ are the stored
causal signs/scales. Each replicate also generates a global null with $\beta=0$.

Install once in R:

```r
install.packages(c("data.table", "glmnet", "susieR"))
```

Run one replicate (two signal strengths and a global null):

```sh
Rscript run_1000g.R --cohort EUR --replicates 1 --out results/1000g_EUR
```

For 20 replicates:

```sh
Rscript run_1000g.R --cohort EUR --replicates 1:20 --out results/1000g_EUR_20
```

DDMM compares **Simes, effective-rank, rank-3, and full-rank chi-square tests**
across seven SNP-group resolutions. SuSiE provides a credible-set comparison.
A selection is a discovered SNP group or credible set. Scores compare these
with the simulated causal loci: power is the fraction of loci detected;
FDR is the average false discovery proportion across replicates.

**Output under `--out`:** `data/EUR/` contains new phenotypes and fitted selections;
`results/` contains `records.csv`, `summary.csv`, `contrasts.csv`, and `tables.md`.
The main score uses genomic interval overlap; `FDR_member` also checks exact causal-SNP membership.

| Option | Meaning | Default |
|---|---|---|
| `--cohort` | Genotype panel: `EUR`, `MIX`, `ALL` | `EUR` |
| `--replicates` | Noise replicate IDs, e.g. `1:20` | `1` |
| `--amplitudes` | Positive integer effect multipliers | `8,4` |
| `--noise-sd` | Noise standard deviation $\sigma$ | `1` |
| `--q` | DDMM FDR level | `0.1` |
| `--seed` | Noise seed = this value + replicate ID | `1000` |

Reuse the same command to resume; use a new `--out` folder after changing settings.
Genotypes and causal locations stay fixed across replicates.
[Original 1000G downloads](DATA_SOURCES.md#download-original-1000g-data) are optional: the required genotype subset is already included.

### Shortcut: score previously fitted selections

```sh
Rscript reproduce.R
```

This skips simulation and fitting, and scores the included DDMM, SuSiE, and
KnockoffGWAS discoveries. Output: `benchmark_1000g/group_supplement/results/`
(CSV scores, tables, PNG/PDF figures). Fresh KnockoffGWAS fitting requires phased
haplotypes and external software; `run_1000g.R` fits DDMM and SuSiE.

## 2. KnockoffGWAS example: compare group-test statistics

Use the separate bundled example genotypes and their seven published group partitions.

Install additional dependencies in R:

```r
install.packages(c("jsonlite", "BiocManager"))
BiocManager::install("snpStats")
```

```sh
Rscript run.R knockoffgwas --shipped --replicates 1:20 \
  --rank effective,3,full --out results/group_tests
```

Fits Simes and the three chi-square variants on the supplied phenotype plus
20 simulated signal/null replicates. Fits: `results/group_tests/`;
automatic tables and figures: `reports/group_tests/`.

Tune `--replicates` and `--rank` in the command; tune `q`,
`phenotype.effect_scale`, and `phenotype.noise_sd` in
[config/knockoffgwas.json](config/knockoffgwas.json).
Use a new `--out` folder after changing settings.

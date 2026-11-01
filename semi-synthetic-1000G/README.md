# Semi-synthetic 1000G

Fixed genotypes and 20 causal effects (`panel.rds`: `X`, `betas$SNP`),
fresh simulated phenotypes per replicate.

| Experiment | Panels | Model |
|---|---|---|
| `ancestry` | AFR300, EAS300, EUR300, MIX300 | $y=aX\beta_0+\epsilon$ |
| `kinship` | MIX600 | $y=aX\beta_0+u+\epsilon$, $u\sim N(0,K)$ |

## Run

```r
install.packages(c("glmnet", "susieR", "jsonlite"))
```

From this folder, run DDMM and SuSiE without external programs:

```sh
Rscript run.R ancestry --replicates 1 --methods DDMM,SuSiE --out results/ancestry_r
Rscript run.R kinship --replicates 1 --methods DDMM,SuSiE --out results/kinship_r
```

Genotype-only nodewise caches are built automatically on the first run.
Use `--replicates 1:50` for repeated simulations; repeat a command to resume.
Settings are in `config/ancestry.json` and `config/kinship600.json`.
Use a new `--out` directory after changing methods or settings.

**Output:** fitted records in `results/`; CSV/LaTeX tables and figures in `reports/`.
To regenerate a report without refitting:

```sh
Rscript report.R ancestry --input results/ancestry_r
```

## Add HMM Knockoff+ and BOLT-LMM

```r
install.packages(c("knockoff", "BiocManager"))
BiocManager::install("snpStats")
```

Install `snpknock2` and BOLT-LMM separately, then set their executable paths:

```sh
export SNPKNOCK_BIN=/path/to/snpknock2
export BOLT_BIN=/path/to/bolt
Rscript run.R ancestry --replicates 1:50
Rscript run.R kinship --replicates 1:50
```

All five panels derive from public 1000 Genomes Phase 3 data; phenotypes are simulated.
BOLT can fail at these small sample sizes; failed fits are marked invalid and
counted separately in `supplement_fit_counts.csv`.

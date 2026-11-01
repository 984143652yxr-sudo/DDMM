# DDMM group

Both workflows use the included data. Run terminal commands from `DDMM group/`.

## 1. 1000G benchmark: regenerate scores and tables

Scores the saved DDMM, SuSiE, and KnockoffGWAS selections for EUR, MIX, and ALL.

```sh
Rscript -e 'install.packages("data.table", repos="https://cloud.r-project.org")'
Rscript reproduce.R
```

**Output:** `benchmark_1000g/group_supplement/results/` — CSV scores, `tables.md`, LaTeX tables, and `benchmark_group.png` / `.pdf`.

## 2. Compare group-level test statistics

Fits DDMM on the bundled KnockoffGWAS example across seven group resolutions. Compares Simes with effective-rank, rank-3, and full-rank chi-square tests; reports power, false discovery proportion, and global-null error rates.

Install once in R:

```r
install.packages(c("data.table", "jsonlite", "glmnet", "BiocManager"))
BiocManager::install("snpStats")
```

Run the supplied phenotype and 20 simulated signal/null replicates:

```sh
Rscript run.R knockoffgwas --shipped --replicates 1:20 \
  --rank effective,3,full --out results/group_tests
```

**Output:** fitted results in `results/group_tests/`; CSV summaries, LaTeX tables, and figures in `reports/group_tests/`. Reports are generated automatically.

To regenerate reports from existing fits:

```sh
Rscript report.R knockoffgwas --input results/group_tests --output reports/group_tests
```

### Parameters to tune

| Parameter | Where | Meaning |
|---|---|---|
| `--replicates 1:20` | Command | 20 simulated signal datasets and 20 global-null datasets; use `1` for a quick run. |
| `--rank effective,3,full` | Command | Chi-square ranks: LD-based effective rank, fixed rank 3, or full group dimension. Simes is always included. |
| `q` | [Config](config/knockoffgwas.json) | BH FDR level within each resolution; default `0.1`. |
| `phenotype.effect_scale` | Config | Simulated effect multiplier; default `-0.256`. |
| `phenotype.noise_sd` | Config | Simulated noise standard deviation; default `1.06`. |
| `--out results/group_tests` | Command | Fit directory. Choose a new directory after changing model settings. |

[Data sources](DATA_SOURCES.md)

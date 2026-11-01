# DDMM group

## 1000G benchmark: generate scores and tables

```sh
Rscript -e 'install.packages("data.table", repos="https://cloud.r-project.org")'
Rscript reproduce.R
```

Uses the included EUR, MIX, and ALL designs and saved DDMM, SuSiE, and
KnockoffGWAS selections. No 1000G download or model refitting is required.
These small selection files are scoring inputs; tables and figures are created when you run the script.

**Output:** `benchmark_1000g/group_supplement/results/` — `tables.md`, CSV scores,
LaTeX tables, and `benchmark_group.png` / `.pdf`.

## Toy example: fit and report

Install once:

```r
install.packages(c("data.table", "jsonlite", "glmnet", "BiocManager"))
BiocManager::install("snpStats")
```

```sh
Rscript run.R knockoffgwas --shipped --no-report
Rscript report.R knockoffgwas
```

**Output:** fits in `results/knockoffgwas/`, tables and figures in `reports/knockoffgwas/`.
For repeated simulations, add `--replicates 1:20` to `run.R`.

## Settings and files

[config/knockoffgwas.json](config/knockoffgwas.json) controls the **toy example**:

| Setting | Meaning |
|---|---|
| `q` | BH false discovery rate level; default 0.1 |
| `phenotype.effect_scale`, `noise_sd` | Simulated signal strength and noise SD |
| `phenotype.seed_base` | Replicate seed is this value + replicate number |
| `output` | Folder for fits; use a new folder after changing settings |

Keep data paths and partition columns unchanged for the supplied example.
Use `--rank effective`, `--rank 3`, or `--rank full` to choose the group-test rank.
The 1000G workflow uses saved DDMM/KnockoffGWAS selections at q = 0.1;
this toy config does not change them.

- `R/`: DDMM implementation and reporting functions.
- `data/knockoffgwas/`: toy genotypes, phenotypes, partitions, and comparator inputs.
- `benchmark_1000g/group_supplement/`: panel data, comparator inputs, and analysis scripts.
- `input_checksums.csv`: input integrity checks used by `reproduce.R`.

[Data sources](DATA_SOURCES.md)

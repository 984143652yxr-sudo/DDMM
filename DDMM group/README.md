# DDMM group

Group-level tests based on the debiased Lasso, compared with KnockoffGWAS and
SuSiE on a compact 1000 Genomes benchmark. The original KnockoffGWAS toy-example
workflow is included separately.

## 1. Reproduce the 1000G scores, tables, and figure

Requires R and **data.table**. No 1000G download, external executable, or model refit is required.

```sh
# Run from DDMM group; quote this directory name when changing directories.
Rscript -e 'install.packages("data.table", repos="https://cloud.r-project.org")'
Rscript reproduce.R --check-only
Rscript reproduce.R
```

The entry point also works when invoked by its full path from another directory.
It checks input hashes and completeness, then runs:

1. `score.R`: scores the saved DDMM, SuSiE, and KnockoffGWAS discoveries.
2. `describe_design.R`: summarizes panel sizes, LD, and group partitions.
3. `diag_false.R`: computes LD diagnostics for false discoveries.
4. `make_tables.R`: generates Markdown and LaTeX tables.
5. `figure.R`: generates the comparison figure in PNG and PDF.

Generated files go into `benchmark_1000g/group_supplement/results/`.
CSV values are checked against `reference_results/` with tolerance `1e-8`;
Markdown/LaTeX tables are compared line by line. Figures are regenerated but
are not required to be byte-identical across graphics platforms. The wrapper
stops on missing inputs or mismatches instead of silently omitting a comparator.

**Read the results:** [tables](benchmark_1000g/group_supplement/reference_results/tables.md) ·
[figure](benchmark_1000g/group_supplement/reference_results/benchmark_group.png) ·
[benchmark design and scoring](benchmark_1000g/group_supplement/README.md).

## 2. Included files

| Path | Purpose |
|---|---|
| `R/` | Current DDMM implementation, scoring, reporting, and callable API |
| `run.R`, `report.R`, `config/` | Original toy-example fit/report entry points |
| `reproduce.R`, `input_checksums.csv` | Complete 1000G rescore workflow and input integrity checks |
| `benchmark_1000g/group_supplement/data/{EUR,MIX,ALL}/` | Compact panel designs, phenotypes, partitions, and cached selections |
| `benchmark_1000g/group_supplement/knockoffs/` | Small KnockoffGWAS discovery lists needed for comparison |
| `benchmark_1000g/group_supplement/code/` | Rescoring, tables, figures, and optional fit/partition code |
| `benchmark_1000g/group_supplement/reference_results/` | Original result snapshots for inspection and comparison |
| `benchmark_1000g/group_supplement/raw_preparation/` | Source scripts documenting the upstream design construction |
| `data/knockoffgwas/`, `results/knockoffgwas/` | Small toy dataset and saved fits used by `report.R` |
| `docs/METHOD.md` | Statistical method and implementation details |

Each panel includes `partitions/`, `design.rds`, `phenotypes.tab`,
`variants_qc.txt`, `samples.*`, `ddmm_selections.rds`, and `susie_selections.rds`.
`design.rds` contains the compact genotype matrix as well as variant metadata,
sample metadata, and simulated causal truth. The phenotypes are simulated.
All copied symbolic links have been replaced by actual files.

Excluded: raw VCFs and phased haplotype files, knockoff genotype matrices,
executables, nodewise caches, cluster launch scripts, local logs/history,
and superseded experimental workflows. `raw_preparation/` is documentation
for rebuilding upstream inputs; it requires external raw data and is not run
by `reproduce.R`.

## 3. Original toy example: run.R and report.R

This is a different experiment from the EUR/MIX/ALL benchmark. Its small
genotypes, published partitions, discovery files, and cached fits are bundled.

```sh
Rscript -e 'install.packages(c("jsonlite", "data.table"), repos="https://cloud.r-project.org")'
# Rebuild the toy report from the saved fits; no fitting required:
Rscript report.R knockoffgwas
```

The report is written to `reports/knockoffgwas/`. To fit the toy data afresh:

```r
install.packages(c("glmnet", "jsonlite", "data.table"))
install.packages("BiocManager")
BiocManager::install("snpStats")
```

```sh
# A new output folder prevents resuming the bundled cached fits.
Rscript run.R knockoffgwas --shipped --replicates 1:20 --out results/new_toy_run
```

Use `Rscript run.R --help` and `Rscript report.R --help` for options.

## 4. Optional refitting

The main reproducibility command recomputes scores from existing selections.
Refitting the 1000G benchmark is a separate, substantially more expensive task:
`code/run_ddmm.R` needs `glmnet`; `code/run_susie.R` needs `susieR`;
rebuilding partitions needs `adjclust`. These scripts can use the compact
`design.rds` files without downloading 1000G. Fit scripts resume existing
selection files, so work on a copy and remove the corresponding cached
selection file there before a fresh fit. Keep partitions and selections from
the same run together. Rebuilding KnockoffGWAS discoveries requires the external
KnockoffGWAS software and phased inputs, which are outside this small release.

References and provenance are recorded in [DATA_SOURCES.md](DATA_SOURCES.md).
The executed checks and software versions are recorded in [docs/VALIDATION.md](docs/VALIDATION.md).

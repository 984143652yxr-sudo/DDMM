# DDMM

Code and reproducibility files for DDMM group-level inference.

Start with **[DDMM group](DDMM%20group/README.md)**. The package includes R code,
compact benchmark inputs, saved method selections, and reference tables.
Scores and tables can be recomputed without downloading 1000 Genomes or rerunning model fits.

```sh
git clone https://github.com/984143652yxr-sudo/DDMM.git
cd "DDMM/DDMM group"
Rscript -e 'install.packages("data.table", repos="https://cloud.r-project.org")'
Rscript reproduce.R
```

See the [reference benchmark tables](DDMM%20group/benchmark_1000g/group_supplement/reference_results/tables.md).

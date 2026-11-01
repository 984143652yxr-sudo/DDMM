# Summary-statistic DDMM

Two R files:

- `ddmm_summary.R`: self-contained corrected DDMM, including Gram-Lasso,
  PCA filtering, per-SNP BIC nodewise tuning, noise estimation, and BH selection.
- `run_ad.R`: downloads public AD data, runs corrected DDMM and GhostKnockoff,
  and produces tables, a comparison figure, and locus-support summaries.

```r
install.packages(c("glmnet", "data.table", "Matrix"))
```

```sh
Rscript run_ad.R          # all chromosomes
Rscript run_ad.R 21,22    # small execution check
```

Results are generated in `output/`; downloaded data stay in `data/`.
Only the corrected DDMM results and the GhostKnockoff comparator are reported.
The chromosome subset is an execution check, not the complete AD analysis.

For other datasets, $S=X^TX/n$, $r=X^Ty/n$, $s_y=y^Ty/n$:

```r
source("ddmm_summary.R")
p <- ddmm_pvalues(Sigma=S, r=r, n=n, n_bic=n, sigma="estimate", s_y=s_y)
```
Use `sigma="estimate"` to estimate noise. For marginal z-scores use `r=z`,
`n=1`, `sigma="unit"`, and set `n_bic` separately.

The AD runner uses `N_BIC=50000` by default, a nodewise tuning setting rather
than a shared GWAS sample-size claim. Example: `N_BIC=5000 Rscript run_ad.R`.
It uses LD shrinkage 0.1 and five equicorrelated GhostKnockoff copies.

[Public AD source](https://github.com/GuJQ5/FWER_Knockoff_filter/tree/main/Section%205)

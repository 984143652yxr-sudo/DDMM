# Summary-statistic DDMM

Inputs: $S=X^TX/n$, $r=X^Ty/n$, and $s_y=y^Ty/n$ when estimating the noise scale.
The R argument for $S$ is named `Sigma`.

```r
install.packages(c("glmnet", "data.table", "Matrix"))
```

Run from this folder:

```sh
Rscript example.R
# Small AR(1) execution check:
AR1_REPS=2 AR1_P=40 AR1_N=80 AR1_K=3 AR1_AMPLITUDES=3,7 NCORES=1 Rscript ddmm_vs_ghostknockoff_ar1_fixed.R
# Full simulation (100 replicates by default):
Rscript ddmm_vs_ghostknockoff_ar1_fixed.R
```

| File | Purpose |
|---|---|
| `ddmm_sumstat.R` | Original implementation |
| `ddmm_sumstat_fixed.R` | Adds per-SNP BIC nodewise tuning and optional noise estimation |
| `ghostknockoff.R` | GhostKnockoff comparator |
| `*_ar1_fixed.R` | AR(1) comparison, original and corrected DDMM on the same draws |
| `*_realdata_fixed.R` | AD comparison |
| `ddmm_ad_replication_fixed.R` | Locus overlap/support tables |

The versions without `_fixed` retain the original analyses. Outputs are written
to `output_fixed/` or `output/` and are excluded from Git.

## Alzheimer's disease data

Public source: [GhostKnockoff Section 5](https://github.com/GuJQ5/FWER_Knockoff_filter/tree/main/Section%205).
The downloader pins an upstream commit and saves candidate z-scores, annotations,
and LD matrices under `data/realdata_section5/`.

```sh
Rscript download_data.R
Rscript ddmm_vs_ghostknockoff_realdata_fixed.R
Rscript ddmm_ad_replication_fixed.R fix
```

For a short execution check, give `21,22` to the first two commands.
That subset is not the full AD analysis; use all chromosomes for pooled testing.
The discovery data are preselected candidates, not an unfiltered genome-wide scan.

On the z-score scale, use `n = 1` and unit noise; `n_bic` controls nodewise
penalty tuning separately. The AD script explores `n_bic = 500, 5000, 50000`
and defaults to 50000 (`N_MAIN`); these are sensitivity settings, not a claim
about a shared GWAS sample size. GhostKnockoff uses `equi` by default;
`Rdsdp` is needed for the optional `sdp` construction.

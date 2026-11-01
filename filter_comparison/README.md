# Filter comparison

Portable adaptations of `extended_model/` and `family_structure_simulation.R`.
Both generate $y=X\beta+u+\epsilon$ with ancestry structure (`bnpsd`);
`family.R` also generates pedigrees (`simfam`).

```r
install.packages(c("bnpsd", "simfam", "popkin", "glmnet", "ggmix", "plmmr"))
```

Run from this folder:

```sh
Rscript admixture.R --smoke
Rscript family.R --smoke
Rscript plot.R admixture
Rscript plot.R family
```

For the larger design, replace `--smoke` with `--n 1000 --p 1200 --reps 20`.
Family generations default to 20; `--seed` controls reproducibility.

| Simulation | Methods |
|---|---|
| Admixture | ggmix, plmmr, PCA-Lasso, DDMM-PCA, DDMM-LMM, DDMM-TRIM |
| Family | LMM-Lasso (ggmix), DDMM-LMM, DDMM-TRIM |

**Output:** `results/<simulation>/replicates.csv`, then `summary.csv` and
`comparison.pdf` after plotting. These files are generated locally.

An 80/20 split is used for fitting/prediction. The phenotype has equal kinship
and iid noise components. For DDMM, the kinship matrix is estimated from
simulated genotypes and the mixture weight is estimated by ggmix; the same
estimated noise covariance is passed through each filter. Discovery uses BH
at 0.1; Lasso baselines select nonzero coefficients. MSPE measures fixed-effect
prediction on held-out outcomes. Smoke runs check execution, not FDR calibration.

These cleaned scripts use the current shared DDMM core and covariance-aware
standard errors; they are not a byte-for-byte replay of the exploratory scripts.

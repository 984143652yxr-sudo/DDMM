# Spectral filters

`filters.R` returns `(QX, Qy)` and the transformed noise covariance.
All filters centre the data; SNP columns are standardised first.

| `method` | Filter before centring | Setting |
|---|---|---|
| `PCA` | $I-U_kU_k^T$ | `k = 3` |
| `TRIM` | singular values $d_i$ are capped at $d_{\lfloor\rho r\rfloor}$ | `rho = 0.3`, $r=\mathrm{rank}(X)$ |
| `LMM` | $[\eta K+(1-\eta)I]^{-1/2}$ | `K`, `eta` |

Run the short example from the repository root:

```sh
Rscript -e 'install.packages("glmnet")'
Rscript examples/filtered_ddmm.R
```

Pass `f$noise_cov` into `variant_tests`: filtering changes the noise covariance
and the standard errors. PCA/TRIM assume iid input noise by default; use
`noise_cov` to supply another covariance shape. LMM defaults to its stated covariance.
The shared test estimates the remaining variance scale. This interface uses one
filter for the pilot and debiasing; it does not reproduce the older two-filter pilot variants.

# DDMM

| Module | Start here |
|---|---|
| [Spectral filters](R/README.md) | PCA, TRIM, LMM and a short DDMM example |
| [DDMM group](DDMM%20group/README.md) | Group tests and compact 1000G rescoring inputs |
| [Filter comparison](filter_comparison/README.md) | Admixture and family simulations |
| [Semi-synthetic 1000G](semi-synthetic-1000G/README.md) | Five genotype panels, simulation runners, and supplied summary tables |
| [Summary statistics](summary_statistics/README.md) | Corrected DDMM and a single AD analysis runner |

```sh
git clone --depth 1 https://github.com/984143652yxr-sudo/DDMM.git
cd DDMM
```

Each module lists its dependencies and run commands. Generated fits, figures,
and caches stay local. The two requested semi-synthetic summary-table folders
are included; public AD data are downloaded separately.

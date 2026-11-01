# 1000G group benchmark

This directory bundles the compact inputs and saved selections needed to
recompute the three-panel comparison. Run `Rscript reproduce.R` from the
`DDMM group` directory. The entry point sets `BENCH` and `DDMM_PKG` for the
individual scripts automatically.

## Panels and simulation

| Panel | Individuals | Composition | Variants |
|---|---:|---|---:|
| EUR | 503 | European-ancestry panel | 2,000 |
| MIX | 504 | 168 EUR + 168 AFR + 168 EAS | 2,000 |
| ALL | 2,504 | All Phase 3 individuals | 2,000 |

Each panel retains 1,000 common biallelic SNPs on chromosome 21 and 1,000 on
chromosome 22. Selection begins at 30 Mb (GRCh37), with MAF ≥ 0.05 within the
panel and spacing ≥ 2 kb. Retained variants and partitions can differ by panel.
The MIX sample draw uses seed 2019.

The saved design has 40 causal variants in 20 loci. There are 20 replicates
at amplitude 8, 20 at amplitude 4, and 20 global-null phenotypes per panel.
The construction is

$$y = 10\,\operatorname{scale}(Z)/\sqrt n + \operatorname{scale}(G)\beta + \varepsilon,
\qquad \varepsilon\sim N(0,I),$$

with $Z=2-\mathrm{sex}$ and causal coefficients
$\beta_j=a\,\mathrm{Sign}_j\,\mathrm{Scale}_j/\sqrt n$, for $a\in\{8,4\}$.
Noncausal coefficients are zero; the null omits the genetic term.
Replicate $r$ uses noise seed $1000+r$. These are simulated outcomes on
real genotype subsets. Sex and an intercept are residualized in the DDMM
and SuSiE fitting scripts.

## Methods and scores

- **Groups:** seven nested adjacency-constrained LD partitions, `res6` (coarse)
  through `res0` (fine), closed under $r^2\ge0.8$. Panel-specific partitions
  are shared by DDMM and KnockoffGWAS.
- **DDMM:** Simes and group chi-squared tests at effective rank, rank 3, and
  full rank; BH at $q=0.1$ and a separate truncated selection rule.
- **KnockoffGWAS:** saved knockoff+ discoveries, with knockoff seed $2019+r$.
- **SuSiE:** one fit per phenotype with $L=40$; 95% and 90% credible sets,
  plus variants with PIP > 0.9. Credible sets use minimum absolute correlation
  purity 0.5. This fit uses the known simulation signal count as its $L$ choice.

For the main comparison, each discovery is scored as its genomic interval:

$$\mathrm{FDP}=\frac{\#\{\text{discovery intervals containing no causal variant}\}}
{\max(1,\#\text{discoveries})},\qquad
\mathrm{power}=\frac{\#\text{causal loci covered}}{20}.$$

Realized FDR is the replicate mean of FDP. Under the global null, FDP is the
indicator of any discovery. For SuSiE, `FDP_member` additionally checks whether
a causal variant belongs to the actual set, since a credible set can be
noncontiguous. `FDP_untagged` is an LD-tagging diagnostic, separate from the
main interval definition. At very coarse resolution, declaring all groups
can already cover all causal loci: read power alongside FDR, width, and the
`All FDP` baseline in the tables.

## Data dictionary

| Input | Contents / use |
|---|---|
| `design.rds` | `X` (individual × SNP dosage matrix), `map`, `samples`, `causal`, `n_loci`, `amplitudes`, `replicates`, `h2` |
| `partitions/part_chr{21,22}.txt` | Variant-to-group assignments at seven resolutions |
| `phenotypes.tab` | Sample identifiers/covariates and 60 simulated phenotype columns |
| `samples.tsv`, `samples.txt` | Public-panel sample metadata/identifiers |
| `variants_qc.txt` | Retained variant identifiers |
| `ddmm_selections.rds` | Selected groups by phenotype, resolution, selection rule, and test |
| `susie_selections.rds` | Credible-set variant indices, PIP selections, and fit diagnostics |
| `knockoffs/<panel>/s<seed>/results/*_discoveries.txt` | KnockoffGWAS discovery intervals; 1,260 small files total |

Saved SNP indices refer to the matching panel's `design.rds`; group indices
refer to that panel's partitions. These files should be kept together.

## Outputs and scope

`reference_results/` contains the supplied result snapshots. Running the
wrapper writes a separate `results/` directory, compares all six CSVs and
four Markdown/LaTeX tables with the snapshots, and records `sessionInfo.txt`.
The PNG/PDF figure is rebuilt from the summaries.

This verifies the scoring and reporting pipeline from saved discoveries.
It does not rerun the model fits or regenerate knockoffs. Optional DDMM/SuSiE
fits and partition generation use the bundled compact design; rebuilding
the raw phased inputs and KnockoffGWAS fits requires additional data/software.

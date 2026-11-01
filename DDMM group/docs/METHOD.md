# DDMM group method

## Group tests

For a phenotype y (centred) and standardised genotypes X:

1. **Pilot.** 5-fold cross-validated Lasso, λ = lambda.min (glmnet); support Ŝ,
   residual r = y − Xβ̂. If |Ŝ| < 3, the top 20 variants by marginal
   correlation are added to the tested set (guard).
2. **Noise level.** σ̂² = ‖r‖² / ‖(I − P_Ŝ)‖²_F = ‖r‖² / (n − rank X_Ŝ).
3. **Single variant** (j ∈ Ŝ). z_j is the nodewise Lasso residual of x_j on all
   other variants, λ chosen per variant by BIC on a 30-point path;
   b_j = β̂_j + z_jᵀr / z_jᵀx_j, T_j = b_j z_jᵀx_j / (σ̂‖z_j‖), two-sided normal
   p-value. Untested variants have p = 1.
4. **Group g.** Z_g = X_g − X_N Γ_g, with Γ_g a nodewise group Lasso of X_g on
   the ±40 neighbouring variants (λ = √(log p / n), 90 FISTA steps, penalty on
   whole neighbouring groups). b_g = β̂_g + (Z_gᵀX_g)⁺ Z_gᵀ r, using the same
   pilot, and Ω_g = (Z_gᵀX_g)⁺ Z_gᵀZ_g (Z_gᵀX_g)⁺ᵀ.
5. **Reduced-rank Wald test.** With C_K the top K_g eigenvectors of the
   group's LD matrix, W_g = (C_Kᵀb_g)ᵀ (C_KᵀΩ_gC_K)⁺ (C_Kᵀb_g) / σ̂², compared
   with χ²_{K_g}; `--rank` sets K_g.
6. **Simes.** min_k |g| p_(k) / k over the tested variants of the group.
7. **Selection.** BH at q = 0.10 within each resolution, separately per track.
   At res0, where every group is one variant, all tracks use the
   single-variant p-values.

**Scoring.** A discovery is the interval from the first to the last variant
of its group. It is true when the interval contains a causal variant. A locus
is found when some discovery contains one of its causal variants; locus power
= loci found / 20, and FDP = false discoveries / max(1, discoveries). Under the
global null FDP = 1{any discovery}, so its mean is the probability of any
discovery.

**Phenotypes.** `shipped` is `phenotypes.tab`. Replicate *r* draws
y = Xβ + ε, ε ~ N(0, 1.06²), with seed 1000 + r and β = −0.2560 × Sign × Scale
on the 40 causal variants (the effects recovered from the shipped phenotype,
h² ≈ 0.71). Its global-null partner is y = ε with the same ε.

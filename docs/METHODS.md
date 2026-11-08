# XCAL and XCAL-t: theory, estimators and limits (v0.2.0)

## 1. Setup and null model
A summary-statistic TWAS test (S-PrediXcan form) combines:
- expression weights `w` on the allele-dosage scale, trained in ancestry *A*;
- GWAS z-scores `z` from ancestry *G*;
- an LD reference panel *r* with SNP standard deviations `s_r` and correlation matrix `R_r`.

The statistic is

    Z = v'z / sqrt(v' R_r v),    v = w ∘ s_r.

**Polygenic null.** Assume per-standardized-SNP effects β with E[β] = 0, Var(β_j) = h²/M, a GWAS of size N, and LD Σ_G in the GWAS population. Marginal z-scores then satisfy

    E[z] = 0,    Cov(z) = a Σ_G + c Σ_G²,    c = N h²/M.

Here *a* is the LD-score-regression intercept (a = 1 without confounding). Under H0 the gene's genetically regulated expression has no effect beyond this background, so

    Var(v'z) = a v'Σ_G v + c v'Σ_G² v.

Var(Z) separates into:
- an LD-mismatch factor, v'Σ_G v / v'R_r v;
- a polygenic factor, 1 + (c/a) q_G, where q_G = v'Σ_G² v / v'Σ_G v.

q_G is a weight-specific LD score **of the GWAS population**. The variance-control slope of Liang et al. (AJHG 2026) is Φ = q_G / M. A Φ computed from European LD therefore mis-scales the correction whenever the GWAS ancestry has different LD.

The covariance result needs only second moments of β. It does **not** make Z Gaussian. When few SNPs carry effects (sparse architectures), v'z is a scale mixture with heavy tails, and no variance-based correction (XCAL, variance control, genomic control, or even the oracle variance) can control extreme-tail error (Section 6).

## 2. The XCAL statistic

    Z_XCAL = v_t'z / sqrt(D̂),    D̂ = a L̂ + c Q̂,    v_t = w ∘ ŝ_t

- *t* is an LD panel **of the GWAS ancestry**, with n individuals and p polymorphic SNPs in the gene window.
- R̂_t is the sample correlation matrix: genotypes centered and scaled by panel means and SDs.
- L̂ = v_t'R̂_t v_t.
- Q̂ estimates v_t'Σ²v_t (Section 3).
- c and a come from LD-score regression on the GWAS, with ancestry-matched LD scores (Section 5).

## 3. Estimating the quadratic form v'Σ²v

**Uncorrected plug-in.** For Gaussian rows with known mean, S = X'X/n satisfies E[S²] = (1 + 1/n) Σ² + (tr Σ / n) Σ. With tr Σ = p for a correlation matrix, the plug-in ‖R̂v‖² overstates v'Σ²v by about (p/n) v'Σv. For 1000 Genomes-sized panels (n ≈ 500, p ≈ 400–4,000 SNPs per window) this bias is of the same order as the quantity itself.

**Moment-corrected estimator.**

    Q̃ = [ ‖R̂v‖² − (p/n) v'R̂v ] / (1 + 1/n)

The identity is exact for Gaussian data with known means. The sample-correlation construction estimates means and scales, and genotypes are discrete, so Q̃ is only *approximately* unbiased. We measured its bias directly (Supplementary Table S1; `scripts/00_estimator_bias.py`; and `table_estimator_genotype_panels.csv`, which compares 500-person panels against 20,000-person LD on simulated genotypes):
- Gaussian data: relative bias within about ±2%.
- Simulated genotype panels: see the genotype table.

**Positivity safeguard (v0.2.0).** Q̃ can be negative in finite samples. We use

    Q̂ = max( Q̃ , L̂² / ‖v‖² ).

The floor is the Cauchy–Schwarz bound v'Σ²v ≥ (v'Σv)² / v'v, which holds for every positive semi-definite Σ. It is therefore consistent and almost never binds. Floor activation rates are reported in every benchmark output.

**Correction from v0.1.0.** v0.1.0 used the floor max(Q̃, L̂), i.e. q̂ ≥ 1. True q can be below 1 when weights have opposite signs in strong LD. Take Σ = [[1, 0.8], [0.8, 1]] and v = (1, −1): q = 0.2, but the old floor forced q = 1 for any panel size. That floor biased Q upward by up to about 20% for opposite-sign weights. A unit test now checks this example.

## 4. XCAL-t: accounting for LD-panel noise

**Derivation.** Write D = a v_t'Σv_t + c v_t'Σ²v_t, the exact null variance of the numerator given the panel. The GWAS z is independent of the panel. So, conditional on the panel, v_t'z / sqrt(D) ~ N(0, 1), whatever the panel. This standardized numerator is therefore independent of the ratio D̂/D, which depends only on the panel.

Approximate D̂/D by a scaled chi-square χ²_ν/ν, matching its first two moments (Satterthwaite):
- E[D̂] ≈ D (Section 3);
- Var(D̂) = τ² D², which gives ν = 2/τ².

Then

    Z_XCAL = [v_t'z / sqrt(D)] / sqrt(D̂/D) ≈ N(0,1) / sqrt(χ²_ν/ν) = t_ν.

**Estimating ν.** τ² is estimated with a delete-group jackknife over panel individuals (K = 20 random groups):

    τ̂² = [(K−1)/K] Σ_k (D̂_(k) − D̄)² / D̂².

Each delete-group replicate D̂_(k) **recomputes the panel means, SDs, correlation matrix and v = w ∘ ŝ** from the retained individuals. `xcal.core.jackknife_moments` does this by re-centering and re-scaling the standardized panel, which is algebraically identical to recomputing from raw genotypes (unit-tested). The simulation and real-genotype pipelines call this same function.

**Behaviour and caveats.**
- The jackknife variance estimate is upward-biased (Efron–Stein inequality). The resulting ν is therefore conservative.
- ν is clipped to [4, 10⁵]. The share of models at the lower clip is reported (`table_xcalt_df.csv`).
- The approximation ignores the dependence between v̂ and R̂, which come from the same panel. The empirical benchmarks measure the net effect.
- p-value: `p_xcal_t = 2 · P(t_ν > |Z_XCAL|)`.

## 5. Estimating c (and a)
`xcal.estimate_c` regresses χ² = z² on LD scores ℓ_j = Σ_k r²_jk, using the standard bias correction r̃² = r² − (1 − r²)/(n − 2). LD scores come from an ancestry-matched panel (`xcal.ld_scores`).
- Weights: 1 / (max(ℓ, 1) · (a + cℓ)²), iterated three times.
- Block-jackknife standard errors over contiguous SNP blocks.
- The intercept can be fixed (`fix_intercept=1.0`).

E[χ²_j] = a + c ℓ_j under the model in Section 1, so the slope is c = N h²/M on the same SNP set as the LD scores. In practice:
- Compute LD scores and estimate c on the same SNP set (for example HapMap3 SNPs), using a panel of the GWAS ancestry.
- Pass the estimated intercept as `intercept` when it clearly exceeds 1.

`scripts/08_estimated_c_benchmark.py` evaluates the complete workflow with c estimated from each simulated genome.

## 6. Benchmarks and what they establish
| Benchmark | Truth | Evidence type |
|---|---|---|
| Simulation null (`03`) | Covariance of the simulated z-scores, built from 20,000-person LD. Exact; no estimator involved. | **Counted** false positives at α down to 2.5e-6, with region-bootstrap CIs. Analytic predictions are shown beside the counts. |
| Individual-level GWAS (`04`) | Realized GWAS z-scores (N = 10,000) | Counted rejections at α = 0.05 and 0.01; empirical/predicted variance |
| Estimated c (`08`) | As in `03`, c re-estimated from each genome | Counted rejections at α ≥ 1e-4; predictions for smaller α |
| Real genotypes (`realdata/`) | Independent half-panel LD, estimated | **Model-based predictions only** (no real phenotypes) |

Architectures: π ∈ {1, 0.1, 0.01, 0.001} is the fraction of SNPs with non-zero effect, at the same total heritability.

## 7. Assumptions and limits
- **Architecture.** Calibration is designed and validated for polygenic backgrounds. Sparse architectures produce heavy-tailed nulls (Section 1) that variance calibration cannot fix (quantified in `03`).
- **Samples.** eQTL and GWAS samples are assumed independent. Sample overlap is not evaluated.
- **Traits.** Continuous traits with an additive model. Case-control traits are not evaluated.
- **Ancestry.** A homogeneous GWAS ancestry with a matching panel. For admixed GWAS, local-ancestry-aware LD is needed (e.g. CADET). The admixed-American real-genotype results are descriptive only.
- **Tails.** Gaussian-tail behaviour of z is assumed, and is checked only by simulation.

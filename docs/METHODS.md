# XCAL-t: theory, estimators and limits (v0.3.0)

For version identity and current qualifications, see [Submission review](SUBMISSION_REVIEW.md). Historical numerical results below are those of v0.3.0.

## 1. Setup and null model
A summary-statistic TWAS test (S-PrediXcan form) combines:
- expression weights **w** on the allele-dosage scale;
- GWAS z-scores **z** from a population G;
- an LD reference panel.

With per-SNP genotype SDs **s**, the numerator is v'z with v = w ∘ s.

**Polygenic null.** Assume mutually uncorrelated standardized per-SNP effects with mean 0 and common variance h²/M, a GWAS of size N, LD Σ_G in the GWAS population and LD-score-regression intercept a. Then

    Cov(z) = a Σ_G + c Σ_G²,   c = N h² / M,
    Var(v'z) = a v'Σ_G v + c v'Σ_G² v.

For the standard statistic Z = v'z / sqrt(v'R_r v), with R_r the reference-panel LD,

    Var(Z) = a · [v'Σ_G v / v'R_r v] · [1 + (c/a) q_G],   q_G = v'Σ_G² v / v'Σ_G v.

**Relation to earlier work.** Liang, Nyasimi & Im (2026, Eq. 3) already give the general slope Φ = q_G / M for an arbitrary LD matrix. The extensions here are:
1. quantifying the error made when q is computed in European LD but applied to a GWAS of another ancestry;
2. estimating q_G from a finite ancestry-matched panel;
3. a reference distribution that accounts for the noise in that estimate;
4. counted-false-positive benchmarks.

The variance result uses only second moments of the effects. When few SNPs carry effects, v'z is a scale mixture with heavier-than-Gaussian tails. Section 6 describes what this implies for the tested null.

## 2. The XCAL statistic

    Z_XCAL = v'z / sqrt(D̂),   D̂ = a L̂ + c Q̂,   L̂ = v'R̂v

Here R̂ is the sample correlation matrix of an LD panel of the GWAS ancestry, with n individuals and p polymorphic SNPs in the window.

**Recommended scaling: GWAS SDs.** Use v = w ∘ s_GWAS, with s_GWAS = sqrt(2 f (1 − f)) from the GWAS allele frequencies (`sd_from_eaf`). For a fixed retained SNP set and externally supplied frequencies, v does not depend on the LD-estimation panel, and the target D = a v'Σv + c v'Σ²v is a fixed number.

**SNPs the panel cannot describe.** A weight SNP that is missing or monomorphic in the LD panel is removed from v (numerator and denominator alike), and the lost share of squared weight is reported as coverage. Keeping it in the numerator while the panel denominator ignores it underestimates D; with panel SDs this happens implicitly because ŝ = 0, but with GWAS SDs it must be done explicitly (`xcal_gene` does so; unit test `test_snp_monomorphic_in_panel_is_dropped_from_numerator_with_gwas_sds`). In the simulation benchmark no gene model had weight on such a SNP; in 1000 Genomes half-panels it occurs.

**Fallback: panel SDs.** If GWAS frequencies are unavailable, v = w ∘ ŝ_panel. Holding this panel-derived vector fixed during deletion is an additional approximation, not a conditional-validity theorem; the package warns about it.

## 3. Estimating v'Σ²v
**Uncorrected plug-in.** For an uncentered Gaussian sample covariance matrix with known population standardization, E[S²] = (1 + 1/n) Σ² + (tr Σ / n) Σ. So the plug-in ‖R̂v‖² overstates v'Σ²v by about (p/n) v'Σv.

**Moment-corrected estimator.**

    Q̃ = [‖R̂v‖² − (p/n) v'R̂v] / (1 + 1/n),
    Q̂ = max(Q̃, L̂²/‖v‖²).

The floor is the Cauchy–Schwarz bound v'Σ²v ≥ (v'Σv)²/v'v, which holds for every positive semi-definite Σ.

**Measured accuracy** for the sample-correlation construction (`scripts/00_estimator_bias.py`; `table_estimator_genotype_panels.csv`):

| Setting | Relative bias of Q̂ | Floor activation |
|---|---|---|
| 48 Gaussian settings (n = 250–1,000; same-sign, opposite-sign, dense and low-eigenvalue weights) | −4.5% to +2.8% | up to 19.4% (hardest setting) |
| 500-person simulated genotype panels vs 20,000-person LD (3,540 models per ancestry) | −0.4% to −1.0% (uncorrected plug-in: +5.2% to +9.0%) | 0% |

This is an empirical characterization of approximate unbiasedness, not a theorem. The v0.1.0 floor at the linear form (which forced q ≥ 1) was biased by up to +113% in the Gaussian grid. The uncorrected plug-in was biased by up to +171%.

## 4. XCAL-t: the reference distribution
**Derivation (GWAS-SD scaling, a and c treated as known).** The GWAS z is independent of the panel, so v'z / sqrt(D) ~ N(0, 1) exactly under the Gaussian null. Its distribution does not involve the panel. D is fixed, and D̂ depends on the panel only through R̂. Hence

    Z_XCAL = [v'z / sqrt(D)] / sqrt(D̂ / D),

with numerator and denominator independent. Approximating D̂/D by χ²_ν / ν, matching its mean (≈ 1, Section 3) and variance τ² = Var(D̂)/D², gives ν = 2/τ² (Satterthwaite) and Z_XCAL ≈ t_ν.

**Estimating τ².** τ² is estimated with a delete-group jackknife over panel individuals (K = 20 random groups). **v is held fixed** and each replicate recomputes the correlation matrix of the retained individuals (`jackknife_moments`). The jackknife therefore targets Var(D̂) for the fixed D (as an approximation; the delete-group jackknife is not exact for this statistic).

**The one-SNP example.** With p = 1, R̂ ≡ 1 and Q̂ = v² (the floor binds), so D̂ = (a + c)v² for every panel. The jackknife variance is then 0, ν reaches its upper clip, and the t reference is effectively normal under the stated Gaussian assumptions. This is a unit test.

v0.2 recomputed v inside the jackknife. That estimated Var(D̂) for a target that itself moved with the panel, and in this example produced needless finite ν.

**Status of the approximation.**
- The scaled-χ² form is an approximation.
- The delete-group jackknife variance is not guaranteed to be conservative for this statistic.
- ν is clipped to [4, 10⁵]; the share of models at each clip is reported.

The rule is therefore an approximation whose type-I behaviour we **validate by counting** (Section 6). We do not claim it is proven conservative.

**What the counting showed (v0.3.0).** With known a and c, XCAL-t was slightly liberal: counted ratios 1.02–1.09 at α = 10⁻⁴ and 0.98–1.22 at α = 2.5e-6, in line with the analytic prediction for the same panels (1.00–1.08 and 1.01–1.13). The excess persisted with 5,000-person panels (model-based: 1.01–1.07 at 5,000; 0.94–1.10 over 500–5,000), which is consistent with contributions from the small negative bias of Q̂ and the scaled-χ² approximation, although this experiment does not isolate their contributions. With a and c estimated leave-one-chromosome-out and their SEs propagated, approximate model-based Gaussian ratios in the nonzero infinitesimal settings were 0.99–1.03 at α = 10⁻⁴ and 0.98–1.05 at α = 2.5e-6.

**Panel-SD fallback.** With v = w ∘ ŝ_panel, the same computation holds v fixed at the panel SDs. The jackknife does not in general establish the ratio-variance target or conditional validity when the vector comes from that same panel. The benchmark reports this variant separately.

**Uncertainty in a and c.** With standard errors se_a and se_c (from `estimate_c`), and the estimates obtained independently of the tested z (Section 5),

    τ² = [Var_jk(D̂) + Q̂² se_c² + L̂² se_a²] / D̂²

(diagonal first-order approximation; a–c covariance and interactions with LD-moment estimation are omitted). A full first-order expression would also contain 2 L̂ Q̂ Cov(â,ĉ); its omission has no universal conservative direction. Users should treat this as an approximation and fit a and c leaving out the chromosome being tested.

## 5. Estimating c and a
`estimate_c` regresses z² on LD scores ℓ_j = Σ_k r̃²_jk, with r̃² = r² − (1 − r²)/(n − 2), computed from a panel of the GWAS ancestry (`ld_scores`).
- Weights: 1/(max(ℓ, 1)(a + cℓ)²), iterated.
- Block-jackknife SEs and the a–c covariance are returned.
- Degenerate input (fewer than 1,000 SNPs, or no variation in LD scores) raises an error.
- With LD scores from a few hundred individuals, the slope is attenuated and the intercept absorbs part of the polygenic signal, so a and c should be used **together**.

**Intercept floor.** The benchmark (`08`) and the Biobank Japan application use max(â, 1); a fitted intercept below one is treated as one. This makes the workflow slightly conservative when there is no polygenic signal.

**Recommended workflow.** Leave-one-chromosome-out: test genes on chromosome k with (â, ĉ) fitted on all other chromosomes, so the fitted parameters are independent of the tested z-scores. `scripts/08_estimated_c_benchmark.py` evaluates this workflow, including confounded GWAS with a = 1.1.

## 6. Benchmarks: what each one establishes
| Benchmark | Truth | Evidence |
|---|---|---|
| Simulation null (`03`) | Covariance of the simulated z-scores, built from 20,000-person LD (exact) | **Counted** rejections, about 10⁸ tests per setting, down to α = 2.5e-6, with 2,000 region-bootstrap resamples. Genomic control is counted exactly in a second deterministic pass using λ from the first pass (the median χ² is read from a 2,000-point grid of counts). |
| Individual-level GWAS (`04`) | Realized z-scores (N = 10,000) | Counted at α = 0.05 and 0.01 |
| Estimated a and c (`08`) | As in `03`; a and c re-estimated leave-one-chromosome-out from each simulated genome | Counted at α ≥ 1e-4 with two-way (genome × region) bootstrap CIs. At smaller α, model-based Gaussian rates are computed for infinitesimal settings. Pilot SEs from the first five evaluated genomes excluding pseudo-chromosome 0 are reused across chromosomes and genomes; this prevents interpreting the entire calculation as exact conditional calibration. |
| Panel size (`06`), power (`07`) | Exact truth from 20,000-person LD | Model-based (analytic) |
| Real genotypes (`realdata/`) | Independent half-panel LD (estimated) | Model-based predictions, no phenotypes |
| Real GWAS (`realdata/gwas_application.py`) | Unknown | Descriptive comparison of methods |

**Sparse architectures.** The tested null redraws the causal configuration (fraction π of SNPs with non-zero effect) in every replicate. Dividing by the exact unconditional SD and using a Gaussian reference then fails in the extreme tail. This shows that variance calibration with a Gaussian reference is insufficient under these point-normal architectures. It does not rule out architecture-aware (for example mixture-tail) inference, and it is not a no-association null conditional on a fixed causal configuration.

## 7. Assumptions and limits
- **Architecture:** polygenic background (Section 6 for sparse).
- **Samples:** independent eQTL and GWAS samples.
- **Traits:** continuous traits.
- **Ancestry:** a homogeneous GWAS ancestry with a matching panel. Admixed GWAS need local-ancestry LD (e.g. CADET).
- **Reference distribution:** Gaussian z-score tails, and the approximate t reference of Section 4.


Review clarification: the fixed-vector argument assumes a fixed retained SNP set. Panel-dependent exclusions and use of panel-derived SDs require additional approximation. The diagonal parameter-uncertainty adjustment omits a–c covariance and interactions with estimated LD moments. Genomic-control bootstrap intervals hold the fitted lambda fixed.

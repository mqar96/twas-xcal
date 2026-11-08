# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.2.1] - 2026-11-08
### Removed
- Obsolete `scripts/03_null_calibration.py` (v0.1.0 analytic benchmark; replaced by `03_null_benchmark.py`, and it imported functions removed in 0.2.0).
- Build artefacts (`twas_xcal.egg-info/`) accidentally included in the 0.2.0 tag.

No change to code behaviour or results.

## [0.2.0] - 2026-11-08
Corrective release responding to an independent audit of v0.1.0. All results were regenerated from scratch.

### Fixed
- **Quadratic estimator floor.** v0.1.0 truncated Q̂ at the linear form (`max(est, lin)`), forcing q ≥ 1 and biasing Q upward by up to about 20% for opposite-sign weights. The floor is now the Cauchy–Schwarz bound L̂²/‖v‖², a true lower bound that almost never binds. The untruncated estimate and floor activation are reported.
- **Single jackknife implementation.** The real-genotype pipeline deleted rows from an already-standardized matrix and kept the full-panel weights. Both pipelines now use `xcal.core.jackknife_moments`, which re-centers and re-scales the retained individuals and recomputes v = w ∘ ŝ. It is unit-tested to equal recomputation from raw genotypes.
- **Genomic-control comparator.** v0.1.0 divided by the median predicted variance. It now uses the conventional λ = median(χ²)/0.4549, computed from the simulated test statistics.
- **Effective dimension p** in the bias correction now counts only polymorphic SNPs.
- **Dependencies.** `bed-reader` is in the new `[real]` extra; the `[sim]` extra requires `scikit-learn>=1.7`. SPDX license metadata.

### Added
- `xcal_gene` / `run_xcal` return XCAL-t inference: statistic, Satterthwaite df, `p_xcal_t`, coverage diagnostics and status. Inputs are validated (duplicates, dimensions, non-finite values), and allele harmonization is available via `harmonize`.
- `estimate_c` returns block-jackknife standard errors and the intercept, with an optional fixed intercept. `ld_scores` provides bias-corrected LD scores.
- **Counted-false-positive benchmark** (`scripts/03`): the simulated GWAS z-scores have an exactly known covariance built from 20,000-person LD. False positives are counted down to α = 2.5e-6 with region-bootstrap CIs, under infinitesimal and sparse architectures.
- Estimated-c end-to-end benchmark (`scripts/08`) and estimator bias study (`scripts/00`).
- Individual-level validation now includes XCAL-t and reports gene–phenotype evaluations correctly.
- Real-genotype pipeline: command-line paths, configuration hash checked on resume, pinned DivPop commit with blob checksums, complete 1000 Genomes MD5 list (`verify_inputs.py`), pooling over all gene models.
- One-command reproduction (`scripts/run_all.sh`) with a results manifest. CI builds and tests the wheel outside the checkout, tests minimum dependency versions, and runs a smoke reproduction.

### Changed
- Claims are narrowed to what the evidence supports (see README and `docs/METHODS.md`).

## [0.1.0] - 2026-11-07
### Added
- `xcal.core`: XCAL gene-level statistic (`xcal_gene`, `run_xcal`), finite-sample-unbiased
  quadratic estimator, LD-score-regression estimate of the polygenic scale `c = N·h²/M`.
- `xcal.calibration`: XCAL-t delete-block jackknife with Satterthwaite t-reference;
  baselines (standard TWAS, European-Φ variance control).
- Multi-ancestry TWAS null simulator (stdpopsim OutOfAfrica_3G09; YRI, CEU, CHB) and
  scripts 01–07 reproducing all simulation results.
- Real-data pipeline (`realdata/`) for 12,453 MESA DivPop models × 1000 Genomes phase 3.
- Result tables and figures for simulation (`results/`) and real data (`results/real/`).
- Unit tests and continuous integration.

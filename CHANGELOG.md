# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.3.0] - 2026-11-08
Second corrective release, responding to the author's re-audit of v0.2.1. All results were regenerated from a clean pipeline with per-step provenance.

### Changed (statistical method)
- **XCAL-t reference distribution (re-audit F1).** Weights are now scaled by **GWAS allele-frequency SDs** (`snp_sd`, `sd_from_eaf`; a `gwas.eaf` or `gwas.sd` column in `run_xcal`). The target null variance D is then fixed. The delete-group jackknife **holds v fixed** and recomputes only the panel correlation, so it estimates Var(D̂/D) as the derivation requires. In the one-SNP case the jackknife variance is now 0 and the reference is normal (unit test). The panel-SD fallback is retained with a warning and reported separately.
- **Uncertainty in a and c (F2).** `xcal_from_moments`, `xcal_gene` and `run_xcal` accept `se_c` and `se_a` and add Q̂²se_c² + L̂²se_a² to the denominator variance. `estimate_c` returns the a–c covariance and rejects degenerate input.
- The approximation is described as an approximation. The claim that it is provably conservative was removed.
- Effect on results: XCAL-t is now close to nominal and slightly liberal in the extreme tail (counted 1.02–1.09 at α = 1e-4, 0.98–1.22 at α = 2.5e-6 with known a and c; 0.98–1.05 at α = 2.5e-6 with LOCO-estimated a and c and propagated SEs). v0.2.x reported conservative ratios (0.68–0.90) because its jackknife over-estimated the denominator variance.

### Changed (benchmarks)
- Estimated-parameter benchmark (F2): **leave-one-chromosome-out** fitting, 300 genomes per setting, confounded GWAS (a = 1.1), propagated SEs and two-way bootstrap CIs. The invalid same-data tail prediction was removed. Conditional rates are reported only where the fitted parameters are independent of the tested statistics.
- **Genomic control is now counted exactly (F6)** in a second deterministic pass at the λ estimated in the first pass. 2,000 bootstrap resamples for every method.
- The comparator is labelled "adapted variance control (EUR q)" (F3).
- Figure 2 and Figure 5 whiskers are the 5th–95th percentiles, as stated in the legends (F12).

### Fixed (software)
- Real-genotype benchmark and Biobank Japan application: weight SNPs monomorphic in the LD panel are now dropped from the GWAS-SD numerator, as `xcal_gene` already did. Found during pre-release verification of v0.3.0 (a few models had 10²–10⁴-fold variance mismatches); both analyses were rerun. The simulation benchmark is unaffected (no model had weight on a panel-monomorphic SNP).
- Parameter and input validation (F7): finite c ≥ 0, intercept > 0, SEs ≥ 0, 2 ≤ K with n ≥ max(20, 2K), one-dimensional finite vectors, and coverage in [0, 1].
- `harmonize` (F8) writes the panel alleles into the effect/other columns, keeps the original labels, is idempotent, normalizes case and rejects duplicate panel SNPs.
- `matplotlib>=3.9` (F9). All text outputs are UTF-8, and `run_all.sh` sets `PYTHONUTF8=1`.
- Resumable steps check configuration hashes. The real-genotype configuration hashes genotype files, model and extra files, and source code (F10). Per-region checkpoints were added to scripts 04, 06, 07 and 08, and atomic writes to 02.
- Provenance (F5): `scripts/provenance.py` records the command, git revision and dirty state, imported module path and version, source hashes, environment and timing for every step and attempt. `MANIFEST.json` includes these records. The v0.2.x manifest's "0.1.0" entry was stale editable-install metadata.

### Added
- Real-GWAS application (`realdata/gwas_application.py`): BioBank Japan height (N = 159,095) with 1000 Genomes East Asian LD and MESA models. a and c are fitted leave-one-chromosome-out; chromosomes are processed in 20-Mb windows (±3 Mb) to bound memory; per-chromosome checkpoints.
- CI: Windows wheel job, minimum optional-dependency job, and a smoke test that follows the two-pass benchmark (F9).
- `05_analyze_null.py --combined` reproduces every null-benchmark table from the archived per-gene raw file alone.

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

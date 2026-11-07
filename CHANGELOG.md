# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

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

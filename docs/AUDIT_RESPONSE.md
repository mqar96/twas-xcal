# Response to the independent audit of v0.1.0 (7 November 2026)

Two audit reports reviewed v0.1.0 (commit `0c0bf90`) and the submission package. Every finding below was checked against the code and accepted. v0.2.0 addresses each one, and all results were regenerated from scratch.

## Software and analysis audit (12 findings)
| # | Finding | Response in v0.2.0 |
|---|---|---|
| 1 (P1) | `max(est, lin)` floor forces q ≥ 1, so Q is persistently biased upward (counterexample Σ=[[1,.8],[.8,1]], v=(1,−1)) | Floor replaced by the Cauchy–Schwarz bound L̂²/‖v‖². The untruncated estimate and floor activation are reported. A bias study covers 48 Gaussian settings plus genotype panels: bias −4.5% to +2.8% (Gaussian) and −0.2% to −0.7% (genotypes). The v0.1 floor reached up to +113%. The counterexample is now a unit test (`test_audit_counterexample_q_below_one`). p counts polymorphic SNPs only. |
| 2 (P1) | Real-genotype jackknife did not re-center or re-scale and kept v fixed | One implementation, `xcal.core.jackknife_moments`, used by both pipelines. It re-centers, re-scales and recomputes v, and is unit-tested to equal recomputation from raw genotypes. Real-data results regenerated. |
| 3 (P1) | Extreme tails were model-based; no empirical XCAL-t validation | New benchmark with exactly known simulated null covariance (20,000-person LD). False positives are counted for every method, including XCAL-t, down to α = 2.5e-6, with about 10⁸ tests per cell and region-bootstrap CIs. Infinitesimal and sparse architectures. Analytic predictions are shown next to the counts. The real-genotype analysis is labelled model-based throughout. |
| 4 (P1) | README claims exceeded the evidence ("any ancestry"; XCAL vs XCAL-t; missing derivation) | Claims narrowed. XCAL (normal reference) is reported as inflated with 500-person panels, and only XCAL-t is recommended. The XCAL-t derivation is in `docs/METHODS.md` §4. Sparse-architecture limits are stated prominently. |
| 5 (P1) | Estimated-c workflow not benchmarked; `estimate_c` lacked SEs and diagnostics | `estimate_c` returns block-jackknife SEs and the intercept, supports a fixed intercept, and no longer has the unused argument. `ld_scores` was added. The end-to-end benchmark (`scripts/08`) shows that LDSC with a free intercept attenuates c by 4–8%. Using (â, ĉ) together restores calibration, and this is the documented workflow. |
| 6 (P1) | `run_all.sh` omitted scripts 06 and 07; no source for the power tables | `run_all.sh` regenerates every table and figure (scripts 00–10) and writes `results/MANIFEST.json`. CI runs a smoke reproduction. |
| 7 (P1) | `bed_reader` undeclared; scikit-learn ≥ 1.4 incompatible with `alphas=20` | `[real]` extra with `bed-reader`; `[sim]` requires `scikit-learn>=1.7`. CI tests minimum versions of the core dependencies. |
| 8 (P2) | "80,000 tests" double-counted | Reported as gene models × phenotype replicates (160,000 correlated evaluations per cell in the new run), with region-bootstrap CIs. |
| 9 (P2) | Weak input contracts; no XCAL-t p-value API | `run_xcal` and `xcal_gene` validate duplicates, dimensions and finite values. They report coverage, exclusions and status, and return the statistic, df and `p_xcal_t`. `harmonize` aligns alleles. |
| 10 (P2) | Real-data reproducibility depended on manual edits; unpinned DivPop; missing SAS .fam checksum | Command-line paths, automatic directory creation, and a configuration hash checked before resuming. DivPop pinned to commit `85100c2` with git-blob checks. A complete 15-file MD5 list from the upstream Zenodo record, verified by `verify_inputs.py`. |
| 11 (P2) | Thin tests; CI tested an editable install with a path insert | 11 tests (estimator bias, counterexample, floor, jackknife equivalence, XCAL-t calibration, API validation, LDSC recovery, GC, harmonization). CI builds the wheel, installs it outside the checkout and asserts it is imported from site-packages. |
| 12 (P2) | Headline numbers mixed settings; macro-averaged pooling | Every headline cites its α and N·h² range. Real-data results pool all gene models, with a macro average reported separately. Numbers were cross-checked by an independent agent against the tables. |
| P3 | Deprecated license metadata | SPDX `license = "MIT"` with `setuptools>=77`. |

## Submission-package audit (17 findings, main items)
| Finding | Response |
|---|---|
| Estimator floor, jackknife, extreme-tail validation | As above |
| GC comparator scaled by median predicted variance | Now the conventional λ = median(χ²)/0.455 from the simulated statistics (`genomic_control`) |
| "XCAL-t matches oracle power" | Removed. XCAL-t's relative power loss is reported by setting (≤ 3% when power ≥ 0.6; up to 23% when power < 0.01). |
| Supplementary table S5 averaged identifiers | Workbook rebuilt directly from analysis tables (S1–S16) with a README sheet |
| Figure 5 plotted α = 1e-4 against a 2.5e-6 headline | The real-genotype figure shows both α = 1e-4 and α = 2.5e-6 |
| Abstract over 200 words | 192 words |
| Outdated submission portal | Wiley Research Exchange (authors.wiley.com/journal/GEPI) |
| No AI-use disclosure | Added to the manuscript and README |

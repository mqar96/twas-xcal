# Audit responses

> **Current qualification (9 November 2026):** the historical v0.3.0 claims below must be read with [the final review clarifications](SUBMISSION_REVIEW.md). In particular, reused pilot SEs prevent treating the estimated-parameter rates as exact conditional calibration; the panel-SD fallback has no general conditional-validity result.

## Re-audit of v0.2.1 (8 November 2026): response in v0.3.0
The author re-audited v0.2.1 and found 12 issues: two P1 and ten P2/P3. All are addressed below. Results were regenerated from scratch with a provenance record for every step (`results/provenance/`); the records show which steps ran on an uncommitted working tree and the source hashes that ran.

| ID | Finding | Response in v0.3.0 |
|---|---|---|
| F1 (P1) | The jackknife estimated Var(D̂) while v = w∘ŝ moved with the panel, so it did not target Var(D̂/D). One-SNP counterexample. | Weights are now scaled by GWAS allele-frequency SDs, so D is fixed. The jackknife holds v fixed and recomputes only the panel correlation, so it estimates Var(D̂/D). The one-SNP case now gives zero jackknife variance and a normal reference (`test_reaudit_one_snp_counterexample_gives_normal_reference`). The t rule is described as an approximation validated by counting; the claimed conservativeness argument was removed. The panel-SD fallback is reported separately in every benchmark. |
| F2 (P1) | Estimated-parameter tail predictions conditioned on parameters fitted from the same z-scores; SEs not propagated; counted experiment too small. | Leave-one-chromosome-out fitting makes the fitted (â, ĉ) independent of the tested statistics. SEs are propagated into ν. There are 300 genomes per setting (about 10⁶ tests), confounded settings (a = 1.1) and two-way genome × region bootstrap CIs. Conditional rates are reported only in the independent-fit design. The paper limits its extreme-tail claims to what these support. |
| F3 | Liang et al. Eq. 3 already gives the general formula; the comparator is an adaptation. | Eq. 3 is credited explicitly. The comparator is labelled "adapted variance control (EUR q)" everywhere, with its construction stated. Conclusions are phrased as "transferring a European-derived q", not as a test of the published fitting procedure. |
| F4 | Thresholds missing next to oracle ratios; conservative vs calibrated; sparse-null wording. | The α is stated wherever an oracle ratio appears. "Conservative control" and "nominal calibration" are distinguished. The sparse result is stated for the tested point-normal null with a redrawn causal configuration. |
| F5 | Manifest recorded twas-xcal 0.1.0. | That entry was stale editable-install metadata. Every step is now wrapped by `scripts/provenance.py` (git revision, dirty flag, imported module path and version, source hashes, environment, command, timing). The manifest embeds these records. |
| F6 | GC grid-interpolated, 300 bootstraps. | GC is counted exactly in a second deterministic pass, with 2,000 bootstraps for every method. |
| F7 | Invalid c, intercept and K accepted. | Validation and targeted tests added. |
| F8 | Harmonized values kept old allele labels. | Labels are updated, originals kept, calls are idempotent; case is normalized and duplicate panel SNPs rejected; tests added. |
| F9 | Matplotlib floor; Windows cp1252 failure. | `matplotlib>=3.9`; UTF-8 for all text output; `PYTHONUTF8=1` in `run_all.sh`. CI now also runs a Windows wheel job and a minimum optional-dependency job (bed-reader 1.0, matplotlib 3.9, pandas 2.2, numpy 1.24, scipy 1.10). |
| F10 | One-command scope; raw-file reproduction; resume validation. | `run_all.sh` is documented as simulation-only, with the real-genotype and real-GWAS steps documented separately. `05_analyze_null.py --combined` rebuilds every null table from the archived raw file. Configuration-hash resume checks cover `03`, the real-genotype calibration and the Biobank Japan application (the latter keys on summary-statistics and .bim contents, .bed size, model file names and source code); scripts 04 and 06–08 resume from checkpoints without a configuration check (documented). |
| F11 | Reviewer and assistance statements. | The AI-use statement states that the author audited the earlier releases and describes the AI assistance; no colleague is credited. |
| F12 | Legend whiskers, METHODS numbers and factorization, CADET authors, release notes, metadata, presentation. | Whiskers set to 5th–95th percentiles; METHODS corrected; all reference author lists regenerated from Crossref (APA 7); old release notes marked superseded; repository About text updated; equations, line numbers and figure fonts improved. |

### What changed in the results
- **XCAL-t is now close to nominal rather than conservative.** v0.2.1 reported 0.68–0.90 of nominal at α = 2.5e-6. With the corrected jackknife, counted ratios are 1.02–1.09 at α = 10⁻⁴ and 0.98–1.22 at α = 2.5e-6 (known a and c), matching the analytic prediction. The excess persists with 5,000-person panels, so it is reported as a property of the approximation.
- **Estimated a and c (LOCO, SEs propagated)**: conditional ratios 0.99–1.03 at α = 10⁻⁴ and 0.98–1.05 at α = 2.5e-6, including confounded GWAS.
- **Genomic control (exact counts)**: 1.1–6.0× at α = 10⁻⁴ and up to 26× at α = 2.5e-6.
- **Sparse architectures**: XCAL-t 0.97–1.01× the exact-variance oracle at α = 10⁻⁴ (v0.2.1: 0.85–0.97×).
- **Power**: XCAL-t loses at most 6.5% relative to the oracle (v0.2.1: up to 23%).
- **Real genotypes** (model-based; noisy half-panel truth puts standard TWAS at 1.28–1.52× at α = 2.5e-6 without signal): XCAL-t 1.10–1.75× at α = 2.5e-6 in the homogeneous panels and up to 3.0× in the admixed American panel (v0.2.1: 0.71–1.02× and 1.46×). The v0.2.1 values benefited from the over-estimated jackknife variance.
- **Pre-release verification found one more defect**: with GWAS-SD scaling, the real-genotype and Biobank Japan scripts kept weight SNPs that are monomorphic in the panel in the numerator. Fixed (as `xcal_gene` already did), unit-tested, and both analyses rerun. The simulations were unaffected.
- **Biobank Japan height**: λ 2.75 (standard) → 1.09 (XCAL-t).
- **A second independent verification** (AI agent, before release) reproduced every null-benchmark and real-genotype table from the raw files and checked every number in the manuscript and README. It found wording and rounding errors, undisclosed choices (intercept floored at one; GC median from a grid), a missing configuration check in the application, and an unequal intercept between methods in the Biobank Japan comparison. All were corrected, and the application was rerun with the same (â, ĉ) for every correction.

## Independent audit of v0.1.0 (7 November 2026)

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
| 12 (P2) | Headline numbers mixed settings; macro-averaged pooling | Every headline cites its α and N·h² range. Real-data results pool all gene models, with a macro average reported separately. Numbers were cross-checked against the tables by an independent AI verification agent. |
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

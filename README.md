# twas-xcal

[![tests](https://github.com/mqar96/twas-xcal/actions/workflows/tests.yml/badge.svg)](https://github.com/mqar96/twas-xcal/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24028643.svg)](https://doi.org/10.5281/zenodo.24028643)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)

**Ancestry-aware null calibration for transcriptome-wide association studies (TWAS).**

`twas-xcal` implements **XCAL-t**: a summary-statistic TWAS test whose null variance includes the polygenic background of the *GWAS* ancestry, with a reference distribution that accounts for LD-panel noise. It also ships a benchmark that **counts** false positives against an exactly known simulated null.

> **v0.2.0 is a corrective release.** v0.1.0 had an undocumented estimator floor, two different jackknife implementations, a non-standard genomic-control comparator, and tail error rates that were predicted rather than counted. All of these are fixed and every result was regenerated. See [CHANGELOG](CHANGELOG.md).

## What the evidence shows (and what it does not)
Results below come from 3,540 trained gene models per ancestry, 500 simulated 1-Mb regions, and about 10⁸ counted null tests per setting. The simulated GWAS z-scores have an exactly known covariance built from 20,000-person LD, so no estimator is involved in the truth.

1. **The polygenic inflation slope Φ is ancestry-specific.** Relative to European LD, the median Φ ratio is 1.04 for East Asian LD (5–95%: 0.62–1.77) and 0.63 for African LD (0.37–1.08).
2. **Under a highly polygenic (infinitesimal) background**, European-Φ variance control is miscalibrated outside Europeans. Genomic control fails everywhere.

   Counted type-I error ÷ nominal at α = 2.5 × 10⁻⁶ (N·h² ≈ 6k–247k):

   | GWAS ancestry | Standard TWAS | Genomic control | Variance control (EUR Φ, EUR LD) | Variance control (EUR Φ, matched LD) | XCAL (normal ref.) | **XCAL-t** |
   |---|---|---|---|---|---|---|
   | European | 14–63,600 | 1.6–26 | 1.2–1.8 | 1.2–1.8 | 1.2–1.8 | **0.71–0.90** |
   | East Asian | 15–67,000 | 1.4–25 | 3.8–61 | 1.6–35 | 1.1–1.8 | **0.68–0.88** |
   | African | 6–34,000 | 1.2–13 | 8–16 | 0.41–1.46 | 1.2–1.9 | **0.81–0.89** |

   With a 500-person LD panel, plain XCAL with a normal reference is **not** calibrated in the tail (up to 1.9×). Use XCAL-t.
3. **Sparse architectures break every variance-based correction**, including the oracle that knows the true null variance. With π = 0.01 (about 4 causal SNPs per Mb), even the oracle is 26–60× inflated at α = 10⁻⁴. These false positives are genes near causal variants (the "null" of a TWAS is not "no nearby signal"). At α = 10⁻⁴, **XCAL-t stays at 0.85–0.97× the oracle** in every ancestry and architecture. European-Φ variance control ranges from 0.38× to 9.2× the oracle.

So in the three simulated ancestries, XCAL-t is calibrated under polygenic backgrounds and stays close to the exact-variance oracle under sparse ones. When the trait has strong local signals, extreme-tail TWAS hits still need locus-level follow-up (conditional analysis, fine-mapping or colocalization).

**Real genotypes (1000 Genomes × 12,478 MESA models; model-based predictions, no real phenotypes).**
- Median Φ ratio relative to European LD:
  - African 0.61 (5–95%: 0.30–1.27)
  - South Asian 0.90
  - East Asian 0.94
  - Admixed American 1.22
- At N·h² = 250k and α = 2.5 × 10⁻⁶, European-Φ variance control with EUR LD is predicted at:
  - 11.8× (African), 5.3× (South Asian), 19.6× (East Asian) and 402× (admixed American);
  - 4.7× in European data, the half-panel noise floor.
- XCAL-t is predicted at 0.71–1.02× in the African, South Asian, East Asian and European panels.
- In the admixed American panel XCAL-t reaches 1.46× at the highest polygenicity, consistent with admixture LD.

See [`realdata/README.md`](realdata/README.md).

Not validated: case-control traits, eQTL/GWAS sample overlap, and admixed GWAS (local-ancestry LD is needed).

## Method in one line
```
Z_XCAL = v'z / sqrt( a·v'R̂v + c·Q̂ ),   v = w ∘ ŝ   (panel of the GWAS ancestry)
p = 2·P(t_ν > |Z_XCAL|),               ν = 2/τ̂²    (delete-group jackknife of the denominator)
```
- Q̂ estimates v′Σ²v with a moment correction for panel size and a Cauchy–Schwarz positivity floor. Measured bias is within ±5% (Gaussian settings) and −0.2% to −0.7% (simulated genotype panels).
- c = N·h²/M and the intercept a come from LD-score regression.

Derivation, estimator, jackknife and limits: [`docs/METHODS.md`](docs/METHODS.md).

## Installation
```bash
pip install "twas-xcal @ git+https://github.com/mqar96/twas-xcal@v0.2.0"
# extras for the benchmarks:
pip install "twas-xcal[sim,real] @ git+https://github.com/mqar96/twas-xcal@v0.2.0"
```
Requires Python ≥ 3.10.

## Quick start
```python
import pandas as pd
from xcal import harmonize, ld_scores, estimate_c, run_xcal

# 1. align weights and GWAS z-scores to the panel's counted allele (a1); drops strand-ambiguous SNPs
weights, rep_w = harmonize(weights_raw, panel_alleles)          # columns: gene, snp, effect_allele, other_allele, weight
gwas, rep_g = harmonize(gwas_raw, panel_alleles, value_col="z")

# 2. polygenic scale c (and intercept) by LD-score regression, with LD scores from a panel of the GWAS ancestry
l = ld_scores(panel_G_ldsc, panel_pos_ldsc)                     # same SNP set as the regression
fit = estimate_c(gwas_ldsc_z, l)                                # dict: c, se_c, intercept, se_intercept

# 3. XCAL-t for every gene on a chromosome
res = run_xcal(weights, gwas, panel_G, panel_snps, panel_pos,
               c=fit["c"], intercept=max(fit["intercept"], 1.0))
res[["gene", "n_snps_used", "weight_coverage", "z_xcal", "df", "p_xcal_t", "status"]]
```
`p_xcal_t` is the recommended p-value. The output also reports coverage (share of squared weight retained), the number of SNPs dropped for being absent from the panel or the GWAS, and the estimator floor flag.

## Reproducing the results
```bash
bash scripts/run_all.sh work 2        # every simulation table and figure in results/ (~3.5 h on 2 cores)
python -m pytest tests                # unit tests (incl. the audit counterexample)
```
`results/MANIFEST.json` lists SHA-256 hashes of all outputs and the package versions used. Real-genotype steps (inputs checksum-verified, model commit pinned) are in [`realdata/README.md`](realdata/README.md).

| Output | Script |
|---|---|
| Counted null calibration, sparse architectures, GC, Φ, estimator check | `03_null_benchmark.py` → `05_analyze_null.py` |
| Gaussian estimator bias | `00_estimator_bias.py` |
| Individual-level GWAS check | `04_individual_level_validation.py` |
| Panel size, power, estimated c | `06`, `07`, `08` → `09_analyze_extra.py` |

## Citation
> Rao, M. Q. A. (2026). *twas-xcal: ancestry-aware null calibration for transcriptome-wide association studies* (v0.2.0). Zenodo. https://doi.org/10.5281/zenodo.24028643

See [`CITATION.cff`](CITATION.cff). The concept DOI above always resolves to the latest version.

**Upstream data and methods:**
- 1000 Genomes Project Consortium (2015), *Nature* 526:68.
- MESA DivPop models: Mogil et al. (2018), *PLoS Genet* 14:e1007586.
- Variance control: Liang, Nyasimi & Im (2026), *AJHG* 113:276, [doi:10.1016/j.ajhg.2025.12.014](https://doi.org/10.1016/j.ajhg.2025.12.014).
- LD-score regression: Bulik-Sullivan et al. (2015), *Nat Genet* 47:291.
- Simulation software: msprime (Baumdicker et al., 2022) and stdpopsim (Adrion et al., 2020).

The MIT license covers this code only. Data and models keep their own terms.

## Author
**Muhammad Qudduse Ahmad Rao**, Independent Researcher, Lahore, Pakistan · ORCID [0009-0001-8731-7080](https://orcid.org/0009-0001-8731-7080) · qudduse1234@gmail.com

Software development and analyses were carried out with substantial assistance from an AI coding assistant (Anthropic's Claude), under the author's direction and verification. See the manuscript's AI-use statement.

Contributions and issues are welcome: [`CONTRIBUTING.md`](CONTRIBUTING.md).

## License
[MIT](LICENSE)

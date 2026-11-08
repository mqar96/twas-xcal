# twas-xcal

[![tests](https://github.com/mqar96/twas-xcal/actions/workflows/tests.yml/badge.svg)](https://github.com/mqar96/twas-xcal/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24028643.svg)](https://doi.org/10.5281/zenodo.24028643)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)

**Ancestry-aware null calibration for transcriptome-wide association studies (TWAS).**

`twas-xcal` implements **XCAL-t**: a summary-statistic TWAS test whose null variance includes the polygenic background of the *GWAS* ancestry, with a reference distribution that accounts for LD-panel noise and for uncertainty in the estimated polygenicity. It also ships a benchmark that **counts** false positives against an exactly known simulated null.

> **v0.3.0** responds to the author's re-audit of v0.2.1. Weights are now scaled by GWAS allele-frequency SDs and the jackknife holds them fixed; the standard errors of the estimated intercept and slope are propagated; polygenicity is estimated leave-one-chromosome-out; genomic control is counted exactly; every step records provenance; and the method is applied to Biobank Japan height. All results were regenerated. See [CHANGELOG](CHANGELOG.md) and [docs/AUDIT_RESPONSE.md](docs/AUDIT_RESPONSE.md).

## What the evidence shows (and what it does not)
Simulation results come from 3,540 trained gene models per ancestry, 500 simulated 1-Mb regions, and 1.06 × 10⁸ counted null tests per setting. The simulated z-scores have an exactly known covariance built from 20,000-person LD, so no estimator is involved in the truth.

1. **The polygenic inflation slope is ancestry-specific.** The slope q of a gene is a weight-specific LD score of the GWAS population (the general slope Φ = q/M of Liang, Nyasimi & Im 2026, Eq. 3). Relative to European LD, its median ratio is 1.04 in East Asian LD (5–95%: 0.62–1.78) and 0.63 in African LD (0.37–1.08).
2. **Under a polygenic (infinitesimal) background**, variance control with European q is miscalibrated outside Europeans and genomic control fails everywhere. Counted type-I error ÷ nominal, N·h² ≈ 6k–247k:

   | GWAS ancestry | α | Standard TWAS | Genomic control | Adapted VC (EUR q, EUR LD) | Adapted VC (EUR q, matched LD) | XCAL (normal ref.) | **XCAL-t** |
   |---|---|---|---|---|---|---|---|
   | European | 10⁻⁴ | 5.9–2,400 | 1.2–5.9 | 1.10–1.40 | 1.10–1.40 | 1.10–1.42 | **1.02–1.07** |
   | European | 2.5 × 10⁻⁶ | 14–64,000 | 1.7–26 | 1.2–1.9 | 1.2–1.9 | 1.2–1.9 | **0.98–1.11** |
   | East Asian | 10⁻⁴ | 6.3–2,500 | 1.2–6.0 | 1.75–9.0 | 1.32–7.1 | 1.10–1.37 | **1.02–1.03** |
   | East Asian | 2.5 × 10⁻⁶ | 16–67,000 | 1.8–25 | 3.8–61 | 1.75–35 | 1.3–1.8 | **1.02–1.22** |
   | African | 10⁻⁴ | 3.3–1,500 | 1.1–3.8 | 1.67–2.13 | 0.39–0.63 | 1.09–1.45 | **1.02–1.09** |
   | African | 2.5 × 10⁻⁶ | 5.9–34,000 | 1.4–13 | 8.1–15.8 | 0.46–1.45 | 1.3–2.2 | **1.04–1.20** |

   XCAL-t is close to nominal but **slightly liberal** in the extreme tail (up to about 1.2× at α = 2.5 × 10⁻⁶). With a 500-person panel, plain XCAL with a normal reference is not calibrated in the tail; use XCAL-t.
3. **With a and c estimated** leave-one-chromosome-out (300 genomes per setting, including confounded GWAS with a = 1.1), XCAL-t with propagated standard errors has conditional type-I ratios of 0.99–1.03 at α = 10⁻⁴ and 0.98–1.05 at α = 2.5 × 10⁻⁶.
4. **Under sparse point-normal architectures, variance calibration with a Gaussian reference is insufficient**: even the oracle that knows the exact null variance fails in the extreme tail. With π = 0.01 (about 4 causal SNPs per Mb) the oracle is 26–61× inflated at α = 10⁻⁴. XCAL-t stays at 0.97–1.01× the oracle at α = 10⁻⁴. When a trait has strong local signals, extreme-tail TWAS hits need locus-level follow-up (conditional analysis, fine-mapping or colocalization).

**Real genotypes (1000 Genomes × MESA models; model-based predictions, no phenotypes).** 12,478 MESA models (P < 0.05) in five 1000 Genomes super-populations; the truth is an independent half-panel, whose own noise puts standard TWAS at 1.28–1.52× at α = 2.5 × 10⁻⁶ even without polygenic signal.
- Median slope ratio vs European LD: African 0.61 (5–95%: 0.30–1.27), South Asian 0.90, East Asian 0.94, admixed American 1.22.
- Predicted type-I ratios at α = 2.5 × 10⁻⁶, N·h² up to 250k: European-q variance control (EUR LD) up to 11.8× (African), 5.3× (South Asian), 19.6× (East Asian), 402× (admixed American) and 4.7× (European); XCAL up to 4.1–15.7×; **XCAL-t 1.10–1.75×** in the African, South Asian, East Asian and European panels and up to 3.0× in the admixed American panel (admixture LD). See [`realdata/README.md`](realdata/README.md).

**Biobank Japan height (N = 159,095; 12,195 MESA models).** Standard TWAS: λ = 2.75 and 445 Bonferroni-significant models. XCAL-t with leave-one-chromosome-out LDSC (â ≈ 1.16, ĉ ≈ 0.063) and propagated SEs: λ = 1.09 and 44 models (32 genes), including *UQCC* (*GDF5* region), *ZBTB38* and *GNA12*. No ground truth exists here; the comparison is descriptive.

Not validated: case-control traits, eQTL/GWAS sample overlap, and admixed GWAS (local-ancestry LD is needed).

## Method in one line
```
Z_XCAL = v'z / sqrt( a·v'R̂v + c·Q̂ ),   v = w ∘ s_GWAS,  s_GWAS = sqrt(2f(1−f))
p = 2·P(t_ν > |Z_XCAL|),               ν = 2/τ̂²,  τ̂² = [Var_jk(D̂) + Q̂²se_c² + L̂²se_a²] / D̂²
```
- R̂ comes from an LD panel of the GWAS ancestry. Q̂ estimates v′Σ²v with a moment correction for panel size and a Cauchy–Schwarz positivity floor (measured bias −0.4% to −1.0% on 500-person genotype panels).
- The delete-group jackknife holds v fixed and recomputes only the panel correlation.
- a and c come from LD-score regression, ideally leaving out the tested chromosome.
- If GWAS allele frequencies are missing, panel SDs are used and inference is conditional on them (the package warns).

Derivation, estimator, jackknife and limits: [`docs/METHODS.md`](docs/METHODS.md).

## Installation
```bash
pip install "twas-xcal @ git+https://github.com/mqar96/twas-xcal@v0.3.0"
# extras for the benchmarks:
pip install "twas-xcal[sim,real] @ git+https://github.com/mqar96/twas-xcal@v0.3.0"
```
Requires Python ≥ 3.10. Tested on Linux (Python 3.10–3.13, including the minimum core and `[real]` dependency versions) and Windows.

## Quick start
```python
from xcal import harmonize, ld_scores, estimate_c, run_xcal

# 1. align weights and GWAS to the panel's counted allele (a1); strand-ambiguous SNPs are dropped
weights, rep_w = harmonize(weights_raw, panel_alleles)                 # gene, snp, effect_allele, other_allele, weight
gwas, rep_g = harmonize(gwas_raw, panel_alleles, value_col="z")       # snp, effect_allele, other_allele, z, eaf

# 2. intercept a and slope c by LD-score regression on the OTHER chromosomes (LD scores from the GWAS ancestry)
l = ld_scores(panel_G_ldsc, panel_pos_ldsc)
fit = estimate_c(gwas_other_chrom_z, l)              # c, se_c, intercept, se_intercept, cov_intercept_c

# 3. XCAL-t for every gene on the held-out chromosome
res = run_xcal(weights, gwas, panel_G, panel_snps, panel_pos,
               c=fit["c"], intercept=fit["intercept"], se_c=fit["se_c"], se_a=fit["se_intercept"])
res[["gene", "n_snps_used", "weight_coverage", "sd_source", "z_xcal", "df", "p_xcal_t", "status"]]
```
`p_xcal_t` is the recommended p-value. `sd_source` is `gwas` when an `eaf` (or `sd`) column was supplied. The output also reports weight coverage, SNPs dropped as absent from the panel or the GWAS, and the estimator floor flag.

## Reproducing the results
```bash
bash scripts/run_all.sh work 2        # every simulation table and figure in results/ (~6 h on 2 cores)
python -m pytest tests                # unit tests (incl. the audit counterexamples)
```
Each step writes a provenance record (git revision, dirty flag, module path, source SHA-256, package versions) to `results/provenance/`. `results/MANIFEST.json` lists SHA-256 hashes of all outputs. Real-genotype and Biobank Japan steps are in [`realdata/README.md`](realdata/README.md).

| Output | Script |
|---|---|
| Counted null calibration, sparse architectures, exact GC, slope transfer, estimator check | `03_null_benchmark.py` (two passes) → `05_analyze_null.py` |
| Gaussian estimator bias | `00_estimator_bias.py` |
| Individual-level GWAS check | `04_individual_level_validation.py` |
| Panel size, power, LOCO estimated a and c | `06`, `07`, `08` → `09_analyze_extra.py` |
| Real genotypes, Biobank Japan | `realdata/real_calibration.py`, `realdata/gwas_application.py` |

The per-gene raw counts of the main benchmark are in `results/raw/null_benchmark_per_gene.csv.gz`; `05_analyze_null.py --combined` rebuilds every table from it.

## Citation
> Rao, M. Q. A. (2026). *twas-xcal: ancestry-aware null calibration for transcriptome-wide association studies* (v0.3.0). Zenodo. [doi:10.5281/zenodo.24062489](https://doi.org/10.5281/zenodo.24062489)

See [`CITATION.cff`](CITATION.cff). All versions: concept DOI [10.5281/zenodo.24028643](https://doi.org/10.5281/zenodo.24028643).

**Upstream data and methods:**
- 1000 Genomes Project Consortium (2015), *Nature* 526:68.
- MESA DivPop models: Mogil et al. (2018), *PLoS Genet* 14:e1007586.
- Biobank Japan height: Akiyama et al. (2019), *Nat Commun* 10:4393 (NBDC hum0014.v15).
- Variance control and the general slope: Liang, Nyasimi & Im (2026), *AJHG* 113:276, [doi:10.1016/j.ajhg.2025.12.014](https://doi.org/10.1016/j.ajhg.2025.12.014).
- LD-score regression: Bulik-Sullivan et al. (2015), *Nat Genet* 47:291.
- Simulation software: msprime (Baumdicker et al., 2022) and stdpopsim (Adrion et al., 2020).

The MIT license covers this code only. Data and models keep their own terms.

## Author
**Muhammad Qudduse Ahmad Rao**, Independent Researcher, Lahore, Pakistan · ORCID [0009-0001-8731-7080](https://orcid.org/0009-0001-8731-7080) · qudduse1234@gmail.com

Software development and analyses were carried out with substantial assistance from an AI coding assistant (Anthropic's Claude), under the author's direction and verification. See the manuscript's AI-use statement.

Contributions and issues are welcome: [`CONTRIBUTING.md`](CONTRIBUTING.md).

## License
[MIT](LICENSE)

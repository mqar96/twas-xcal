# twas-xcal

**Ancestry-aware null calibration for transcriptome-wide association studies (TWAS)**, plus a reproducible multi-ancestry TWAS null simulator.

> Status: v0.1.0. Validated in simulation (phase 1) and on real genomes: 12,453 published MESA expression models × 2,490 1000 Genomes genomes from 5 ancestries (phase 2).

## The problem
TWAS tests (PrediXcan / S-PrediXcan / FUSION) are inflated under the null when the trait is polygenic. A gene-specific *variance-control* correction was proposed (Liang, Nyasimi & Im, AJHG 2026). Its inflation slope Φ depends on the LD of the GWAS population, and so far it has been validated only in European samples.

We show at scale that Φ is ancestry-specific, and that applying European-derived corrections to East Asian or African GWAS produces large miscalibration:
- up to 54× excess false positives at genome-wide significance in East Asian GWAS
- over-correction and lost power in African GWAS

## XCAL in one line
```
Z_XCAL = v'z / sqrt( v'R_t v + c · Q_t(v) ),   v = w ∘ s_t
```
- `R_t`, `s_t`: LD and SNP SDs from a panel **of the GWAS ancestry**
- `Q_t`: finite-sample-unbiased estimate of v'Σ²v
- `c = N·h²/M`: estimated from the GWAS by LD-score regression

**XCAL-t** adds a jackknife Satterthwaite t-reference that absorbs LD-panel noise. It is calibrated even with 500-person (1000 Genomes-sized) panels. Theory: [`docs/METHODS.md`](docs/METHODS.md).

## Install & use
```bash
pip install -e .            # core: numpy, scipy, pandas
pip install -e ".[sim]"     # + simulator (msprime, stdpopsim, scikit-learn)
```
```python
from xcal.core import run_xcal, estimate_c
c = estimate_c(gwas_z, ld_scores)["c"]                   # LDSC slope (ancestry-matched LD scores)
res = run_xcal(weights_df, gwas_df, panel_G, panel_snps, panel_pos, c=c)
```

## Phase-1 results (simulation; see `results/`)

| GWAS ancestry | Standard TWAS, EUR LD* | Variance control, EUR Φ + EUR LD* | XCAL-t* |
|---|---|---|---|
| European | 114–2,388× | 1.2–1.5× | 0.92–0.93× |
| East Asian | 118–2,348× | 3.7–8.6× | 0.89–0.90× |
| African | 45–1,361× | 1.9–2.3× | 0.92–0.94× |

\*Type-I error / nominal at α = 1e-4, N·h² ≈ 30k–245k. Based on 3,538 gene models, 500 regions and 74,298 settings. Validated against individual-level GWAS (N = 10,000 × 3 ancestries, 80,000 tests per cell; empirical/predicted variance = 0.994).

Power: XCAL matches oracle power. In African-ancestry GWAS it is about 1.9× more powerful than EUR-Φ variance control.

## Real-data results (phase 2; see `results/real/`)
- Inflation slope Φ relative to European LD, median (5–95%), across 12,453 MESA models:
  - African 0.60 (0.29–1.35)
  - South Asian 0.89 (0.53–1.42)
  - East Asian 0.90 (0.42–1.82)
  - Admixed American 1.24 (0.69–3.66)
- At genome-wide significance with N·h² = 250k, European-Φ variance control is miscalibrated 4.5× (SAS), 12× (AFR), 15× (EAS) and 402× (AMR).
- XCAL-t is 0.63–0.85× nominal in all five ancestries.

## Reproduce
```bash
bash scripts/run_all.sh work     # about 1.5 h on 2 CPU cores
python -m pytest tests
```

## Citation
See `CITATION.cff`. Author: Muhammad Qudduse Ahmad Rao (GitHub: [@mqar96](https://github.com/mqar96)). License: MIT.

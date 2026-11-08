# Real-genotype benchmark (1000 Genomes × MESA)

This benchmark uses **real genotypes and real published expression models, but no real phenotypes**. Each 1000 Genomes super-population is split into two random halves:
- **Half A** is the LD panel the methods use.
- **Half B** stands in for the GWAS-population LD. Its LD moments are estimated with the same bias-corrected estimator.

The resulting type-I error ratios are therefore **model-based predictions**, not counted false positives. Counted false positives with an exactly known truth are in the simulation benchmark (`scripts/03_null_benchmark.py`).

| Step | Command | Notes |
|---|---|---|
| 1. Download | 1000 Genomes phase 3 PLINK files `1000G_{AFR,EUR,EAS,SAS,AMR}.{bed,bim,fam}` from Zenodo [10.5281/zenodo.6614170](https://doi.org/10.5281/zenodo.6614170) (CC BY 4.0) | |
| 2. Verify | `python realdata/verify_inputs.py /path/to/geno_dir` | MD5 list from the upstream record, all 15 files |
| 3. Models | `git clone https://github.com/WheelerLab/DivPop && git -C DivPop checkout 85100c2a41e2022e1e684e165980c3e9868db925`, then `python realdata/read_divpop.py DivPop models/` | Checks the git blob SHA-1 of each model database |
| 4. Run | `python realdata/real_calibration.py --geno /path/to/geno_dir/1000G_{} --models models --out out/real` | About 2 h on 2 cores; peak RAM about 3 GB. Resumes per chromosome only when the configuration hash matches. |
| 5. Analyse | `python realdata/analyze_real.py out/real.csv.gz results/real` | |

The MESA models are distributed by their authors (Mogil et al., *PLoS Genetics* 2018) under their own terms. This repository's MIT license covers only the code.

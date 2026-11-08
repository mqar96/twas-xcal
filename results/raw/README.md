# Per-gene raw outputs

- `null_benchmark_per_gene.csv.gz`: one row per gene model, GWAS ancestry, architecture (π) and c from `scripts/03_null_benchmark.py`. It holds:
  - counts `n_<method>_<alpha>` out of `reps` replicates;
  - analytic predictions `pred_<method>_<alpha>`;
  - panel moments, the Satterthwaite df and floor flags.

  Exact genomic-control counts from the second pass are merged in (`n_GC_<alpha>`). `scripts/05_analyze_null.py --combined` rebuilds every null-benchmark table from this file alone.
- `real_genotype_per_model.csv.gz`: one row per MESA model and 1000 Genomes super-population from `realdata/real_calibration.py`. It holds split-half moments, jackknife replicates and slope ratios, and is the input to `realdata/analyze_real.py`.

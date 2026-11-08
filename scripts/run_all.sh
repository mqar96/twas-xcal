#!/usr/bin/env bash
# Regenerates every simulation table and figure in results/ from an empty working directory.
# Usage: bash scripts/run_all.sh WORKDIR [JOBS]
# Time on 2 CPU cores: about 3.5 h. Real-genotype steps: see realdata/README.md.
set -euo pipefail
W=${1:-work}; J=${2:-2}
S=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$S/.." && pwd)/results
mkdir -p "$W" "$R"
python "$S/00_estimator_bias.py" "$R/table_estimator_bias_gaussian.csv"
python "$S/01_simulate_genotypes.py" "$W/geno" --n-regions 500 --jobs "$J"
python "$S/02_train_weights.py" "$W/geno" "$W/weights" --genes-per-region 8 --jobs "$J"
python "$S/03_null_benchmark.py" "$W/geno" "$W/weights" "$W/bench" --reps 30000 --jobs "$J"
python "$S/05_analyze_null.py" "$W/bench" "$R"
python "$S/04_individual_level_validation.py" "$W/indiv.csv" --regions 40 --n-gwas 10000 --reps 1000
python "$S/06_panel_size_experiment.py" "$W/panel.csv" --regions 40
python "$S/07_power.py" "$W/power.csv" --regions 60
python "$S/08_estimated_c_benchmark.py" "$W/geno" "$W/weights" "$W/estc.csv" --genomes 100
python "$S/09_analyze_extra.py" --indiv "$W/indiv.csv" --panel "$W/panel.csv" --power "$W/power.csv" --estc "$W/estc.csv" "$R"
cp "$W/indiv.csv" "$R/individual_level_validation_raw.csv"
cp "$W/power.csv" "$R/power_raw.csv"
cp "$W/panel.csv" "$R/panel_size_raw.csv"
python "$S/10_manifest.py" "$R"

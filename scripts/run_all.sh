#!/usr/bin/env bash
# Regenerates every SIMULATION table and figure in results/ from an empty working directory.
# Real-genotype results: see realdata/README.md (separate, needs 1000 Genomes + MESA inputs).
# Each step is wrapped by scripts/provenance.py (commit, module path, source hashes, environment, timings).
# Usage: bash scripts/run_all.sh WORKDIR [JOBS]        Time on 2 CPU cores: about 4.5 h.
set -euo pipefail
export PYTHONUTF8=1
W=${1:-work}; J=${2:-2}
S=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$S/.." && pwd)/results
mkdir -p "$W" "$R"
P() { local step=$1; shift; python "$S/provenance.py" "$W" "$step" -- python "$@"; }
P 00_estimator_bias "$S/00_estimator_bias.py" "$R/table_estimator_bias_gaussian.csv"
P 01_simulate "$S/01_simulate_genotypes.py" "$W/geno" --n-regions 500 --jobs "$J"
P 02_train "$S/02_train_weights.py" "$W/geno" "$W/weights" --genes-per-region 8 --jobs "$J"
P 03_null_benchmark "$S/03_null_benchmark.py" "$W/geno" "$W/weights" "$W/bench" --reps 30000 --jobs "$J"
P 05_gc_lambda "$S/05_analyze_null.py" "$R" --benchdir "$W/bench" --lambda-only
P 03_gc_exact_counts "$S/03_null_benchmark.py" "$W/geno" "$W/weights" "$W/bench" --reps 30000 --jobs "$J" --gc-lambda "$R/table_gc_lambda.csv"
P 05_analyze_null "$S/05_analyze_null.py" "$R" --benchdir "$W/bench"
P 04_individual_level "$S/04_individual_level_validation.py" "$W/indiv.csv" --regions 40 --n-gwas 10000 --reps 1000
P 06_panel_size "$S/06_panel_size_experiment.py" "$W/panel.csv" --regions 40
P 07_power "$S/07_power.py" "$W/power.csv" --regions 60
P 08_estimated_c "$S/08_estimated_c_benchmark.py" "$W/geno" "$W/weights" "$W/estc.csv" --genomes 300
P 09_analyze_extra "$S/09_analyze_extra.py" --indiv "$W/indiv.csv" --panel "$W/panel.csv" --power "$W/power.csv" --estc "$W/estc.csv" "$R"
cp "$W/indiv.csv" "$R/individual_level_validation_raw.csv"
cp "$W/power.csv" "$R/power_raw.csv"
cp "$W/panel.csv" "$R/panel_size_raw.csv"
mkdir -p "$R/provenance" && cp "$W"/provenance/*.json "$R/provenance/"
python "$S/10_manifest.py" "$R"

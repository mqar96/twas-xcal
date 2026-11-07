#!/usr/bin/env bash
set -euo pipefail
W=${1:-work}
mkdir -p "$W" results
python scripts/01_simulate_genotypes.py "$W/geno" 500
python scripts/02_train_weights.py "$W/geno" "$W/weights" 8
python scripts/03_null_calibration.py "$W/geno" "$W/weights" results/null_calibration_raw.csv
python scripts/05_analyze_null.py results/null_calibration_raw.csv results
python scripts/04_individual_level_validation.py results/individual_level_validation_N10000.csv 20 10000 500

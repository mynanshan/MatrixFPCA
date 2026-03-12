#!/usr/bin/env bash

set -euo pipefail

# Usage:
#   ./run_simu_simple.sh
#   ./run_simu_simple.sh 50
#
# Positional argument:
#   n = number of subjects (default: 30)

N="${1:-30}"
SCRIPT="experiment/simu-simple.R"

for seed in $(seq 1 100); do
#   for p in 10 30 50; do
  for p in 30; do
    for sig in 0 0.1 0.5; do
      sbatch "$SCRIPT" \
        --n_subjects "$N" \
        --d_row "$p" \
        --d_col "$p" \
        --noise_lvl "$sig" \
        --seed "$seed"
    done
  done
done
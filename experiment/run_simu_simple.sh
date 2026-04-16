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

# check dimension
sig=0.1
distr="normal"
for seed in $(seq 1 100); do
  for p in 10 30 50; do
      sbatch \
        --job-name="simu_n${N}_p${p}_sig${sig}_d${distr}_sd${seed}" \
        --output="logs/simu_n${N}_p${p}_sig${sig}_d${distr}_sd${seed}_%j.out" \
        "$SCRIPT" \
        --n_subjects "$N" \
        --d_row "$p" \
        --d_col "$p" \
        --noise_lvl "$sig" \
        --scores_distr "$distr" \
        --seed "$seed"
    done
done

# check noise level
p=30
distr="normal"
for seed in $(seq 1 100); do
    for sig in 0 0.1 0.5; do
      sbatch \
        --job-name="simu_n${N}_p${p}_sig${sig}_sd${seed}" \
        --output="logs/simu_n${N}_p${p}_sig${sig}_sd${seed}_%j.out" \
        "$SCRIPT" \
        --n_subjects "$N" \
        --d_row "$p" \
        --d_col "$p" \
        --noise_lvl "$sig" \
        --scores_distr "$distr" \
        --seed "$seed"
    done
done

# check distr
p=30
sig=0.1
for seed in $(seq 1 100); do
    for distr in "normal" "uniform" "beta32" "t5"; do
      sbatch \
        --job-name="simu_n${N}_p${p}_sig${sig}_sd${seed}" \
        --output="logs/simu_n${N}_p${p}_sig${sig}_sd${seed}_%j.out" \
        "$SCRIPT" \
        --n_subjects "$N" \
        --d_row "$p" \
        --d_col "$p" \
        --noise_lvl "$sig" \
        --scores_distr "$distr" \
        --seed "$seed"
    done
done
#!/usr/bin/env bash

set -euo pipefail

# Usage:
#   ./run_simu_simple.sh
#   ./run_simu_simple.sh 50
#
# Positional argument:
#   n = number of subjects (default: 30)

reps="${1:-100}"
N="${2:-30}"
SCRIPT="experiment/simu-simple.R"


# check dimension
sig=0.1
distr="normal"
for ((i=1; i<=reps; i++)); do
  for p in 10 30 50; do
      sbatch \
        --job-name="simu_n${N}_p${p}_sig${sig}_d${distr}_sd${i}" \
        --output="logs/simu_n${N}_p${p}_sig${sig}_d${distr}_sd${i}_%j.out" \
        --time="1:00:00" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=2 \
        --mem="8G" \
        "$SCRIPT" \
        --n_subjects "$N" \
        --d_row "$p" \
        --d_col "$p" \
        --noise_lvl "$sig" \
        --scores_distr "$distr" \
        --seed "$i"
    done
done

check noise level
p=30
distr="normal"
for ((i=1; i<=reps; i++)); do
    for sig in 0 0.1 0.5; do
      sbatch \
        --job-name="simu_n${N}_p${p}_sig${sig}_d${distr}_sd${i}" \
        --output="logs/simu_n${N}_p${p}_sig${sig}_d${distr}_sd${i}_%j.out" \
        --time="1:00:00" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=2 \
        --mem="8G" \
        "$SCRIPT" \
        --n_subjects "$N" \
        --d_row "$p" \
        --d_col "$p" \
        --noise_lvl "$sig" \
        --scores_distr "$distr" \
        --seed "$i"
    done
done

# check distr
p=30
sig=0.1
for ((i=1; i<=reps; i++)); do
    for distr in "normal" "uniform" "beta32" "t5"; do
      sbatch \
        --job-name="simu_n${N}_p${p}_sig${sig}_d${distr}_sd${i}" \
        --output="logs/simu_n${N}_p${p}_sig${sig}_d${distr}_sd${i}_%j.out" \
        --time="1:00:00" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=2 \
        --mem="8G" \
        "$SCRIPT" \
        --n_subjects "$N" \
        --d_row "$p" \
        --d_col "$p" \
        --noise_lvl "$sig" \
        --scores_distr "$distr" \
        --seed "$i"
    done
done
#!/bin/bash

module load r/4.5.0

rm test/results/*

n_subjects=${1:-50}
reps=${2:-100}

for ((i=1; i<=reps; i++)); do
    echo "Iteration $i"
    sbatch "test/test-cov-phasechange.R" \
        --n_subjects $n_subjects \
        --seed $i
done
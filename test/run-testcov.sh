#!/bin/bash

module load r/4.5.0

n_subjects=${1:-50}
reps=${2:-100}

for ((i=1; i<=reps; i++)); do
    echo "Iteration $i"
    sbatch \
    --job-name="test-cov" \
    --output="output_files/test_cov_%j.out" \
    --nodes=1 \
    --ntasks=1 \
    --wrap="Rscript test/test-cov-phasechange.R \
        --n_subjects $n_subjects \
        --seed $i"
done
#!/bin/bash
#SBATCH --job-name=simu-simple        # Job name
#SBATCH --output=logs/simu_simple_%j.out         # Standard output file (%j expands to jobID)
#SBATCH --time=6:00:00                   # Maximum runtime (hh:mm:ss)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G                          # Memory allocation

module load r/4.5.0

# Usage check: at least one argument for the R script and one for the number of repetitions.
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <number_of_repetitions> [n_subjects] [drow] [dcol] [npc] [nfea] [init] [opt] [ls]"
  exit 1
fi

n_subjects=${2:-50}
rt=${3:-3}
reps=${4:-10}

for p in 10 20 50 100 200 400 600 800 1000; do
  # Loop to run the R script multiple times.
  for ((i=1; i<=reps; i++)); do
    echo "Running simulation iteration $i..."
    sbatch experiment/simu-simple-full.R --n_subjects "$n_subjects" \
                           --drow "$p" \
                           --dcol "$p" \
                           --rt "$rt" \
                           --repeat "$reps"
  done
done


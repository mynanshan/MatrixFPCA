#!/bin/bash
#
# run_simple.sh

module load r/4.5.0

# defaults: ./run_simple.sh [n_subjects] [rt] [reps]
n_subjects=${1:-30}
rt=${2:-3}
reps=${3:-1}

# where your helper lives (adjust path as needed)
GET_CORES_SCRIPT="$(dirname "$0")/get_num_cores.R"

for p in 10 20 50 100 200 400 600 800 1000 1200 1400 1600 1800 2000; do
  # ask R how many cores for this p
  num_cores=$(Rscript "$GET_CORES_SCRIPT" "$p")
  
  for ((i=1; i<=reps; i++)); do
    echo "[$(date +%T)] p=$p, iteration $i, requesting $num_cores cores..."
    sbatch \
      --cpus-per-task="$num_cores" \
      --wrap="Rscript experiment/simu-simple-full.R \
        --n_subjects $n_subjects \
        --drow $p \
        --dcol $p \
        --rt $rt \
        --repeat $reps"
  done
done

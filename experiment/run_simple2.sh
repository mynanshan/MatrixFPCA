#!/bin/bash
#
# run_simple.sh

module load r/4.5.0

# defaults: ./run_simple.sh [n_subjects] [rt] [reps]
rt=${2:-3}
reps=${3:-10}

# where your helper lives (adjust path as needed)
GET_CORES_SCRIPT="$(dirname "$0")/get_num_cores.R"

for p in 10 20 50 100 150 200 250 300; do
  for ((n=10; n<=100; n+=10)); do
    for ((s=1; s<=6; s++)); do
      # ask R how many cores for this p
      num_cores=$(Rscript "$GET_CORES_SCRIPT" "$p")
      for ((i=1; i<=reps; i++)); do
        echo "[$(date +%T)] p=$p, iteration $i, requesting $num_cores cores..."
        sbatch \
          --cpus-per-task="$num_cores" \
          --job-name="simu-simple" \
          --output="output_files/simu_simple_%j.out" \
          --time="3:00:00" \
          --nodes=1 \
          --ntasks=1 \
          --mem-per-cpu=4G \
          --wrap="Rscript experiment/simu-simple-full.R \
            --n_subjects $n \
            --s_obs_min $s \
            --d_row $p \
            --d_col $p \
            --rank_time $rt \
            --repeat $reps"
      done
    done
  done
done

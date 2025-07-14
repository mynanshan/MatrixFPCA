#!/bin/bash
#
# run_simple.sh

module load r/4.5.0

# defaults: ./run_simple.sh [n_subjects] [rt] [reps]
rt=${2:-3}
reps=${3:-10}

# where your helper lives (adjust path as needed)
GET_CORES_SCRIPT="$(dirname "$0")/get_num_cores.R"

echo "[$(date +%T)]"

for p in 10 20 50 100 150 200 250 300; do
  num_cores=$(Rscript "$GET_CORES_SCRIPT" "$p")
  echo "> Dimension $p: Requesting $num_cores cores..."
  for ((n=10; n<=100; n+=10)); do
    echo ">> Num. subjects: $n"
    for ((s=1; s<=6; s++)); do
      echo ">>> Sampling level: $s"
      # ask R how many cores for this p
      for ((i=1; i<=reps; i++)); do
        echo "- Iteration $i"
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

#!/bin/bash
#
# run_simple.sh

module load r/4.5.0

# defaults: ./run_simple.sh [n_subjects] [rt] [reps]
n_subjects=${1:-20}
rt=${2:-3}
reps=${3:-10}

# where your helper lives (adjust path as needed)
GET_CORES_SCRIPT="$(dirname "$0")/get_num_cores.R"
GET_MEM_SCRIPT="$(dirname "$0")/get_mem_per_core.R"

for p in 10 20 50 100 200 400 600 800 1000 1200 1400 1600 1800 2000; do
  # ask R how many cores for this p
  num_cores=$(Rscript "$GET_CORES_SCRIPT" "$p")
  mem_per_core=$(Rscript "$GET_MEM_SCRIPT" "$p")
  echo ">>>> [$(date +%T)] p=$p, requesting $num_cores cores, ${mem_per_cores}G per core..."

  for s in 1 2 3 4; do
    echo ">> Sampling level $s"

    for ((i=1; i<=reps; i++)); do
      echo "Iteration $i"
      sbatch \
        --cpus-per-task="$num_cores" \
        --job-name="simu-simple" \
        --output="output_files/simu_simple_%j.out" \
        --time="6:00:00" \
        --nodes=1 \
        --ntasks=1 \
        --mem-per-cpu="${mem_per_core}G" \
        --wrap="Rscript experiment/simu-simple-full.R \
          --n_subjects $n_subjects \
          --s_obs_min $s \
          --d_row $p \
          --d_col $p \
          --rank_time $rt \
          --seed0 $i \
          --repeat $reps"
    done
  done
done

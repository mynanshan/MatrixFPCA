#!/bin/bash
#SBATCH --account=def-jiguocao_cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --job-name=echonet-prep
#SBATCH --output=application/echonet/output/slurm-prep-%j.log
# Stage 1: preprocessing for all three alignment modes + go/no-go exploration.
#   mkdir -p application/echonet/output && sbatch application/echonet/run_prepare.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load python/3.11
source application/echonet/.venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/mpl"
OUT=application/echonet/output
for mode in none rigid similarity; do
  python -u application/echonet/prepare_echonet.py --align "$mode" \
    --out "$OUT/prep_$mode" --workers "$SLURM_CPUS_PER_TASK"
done
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
python -u application/echonet/explore_echonet.py \
  "$OUT/prep_none" "$OUT/prep_rigid" "$OUT/prep_similarity" \
  --report application/echonet/report

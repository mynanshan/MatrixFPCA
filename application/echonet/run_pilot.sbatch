#!/bin/bash
#SBATCH --account=def-jiguocao_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --job-name=echonet-pilot
#SBATCH --output=application/echonet/output/slurm-pilot-%j.log
# Stage 2: RDS conversion + pilot MatFPCA fits (56 x 56, n = 500 TRAIN,
# 300 VAL, all TEST; similarity alignment, gain-normalised), sparse vs dense
# comparison, and copies of tables/figures into report/pilot/.
#   sbatch application/echonet/run_pilot.sbatch   (after run_prepare.sbatch)
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load python/3.11 r/4.5.0
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 FLEXIBLAS_NUM_THREADS=8
OUT=application/echonet/output
EXPORT=$OUT/export_pilot56
REPORT=application/echonet/report/pilot
application/echonet/.venv/bin/python -u application/echonet/export_echonet.py \
  --prep $OUT/prep_similarity --out $EXPORT --pool 2 --n-train 500 --n-val 300 --n-test 0
Rscript application/echonet/echonet_to_rds.R $EXPORT
for design in sparse dense; do
  Rscript application/echonet/echonet_matfpca.R --input=$EXPORT \
    --output=$OUT/pilot56_$design --design=$design --cache=false
done
mkdir -p $REPORT
Rscript application/echonet/compare_pilot_fits.R $OUT/pilot56_sparse $OUT/pilot56_dense \
  $REPORT/sparse_vs_dense.csv
for design in sparse dense; do
  mkdir -p $REPORT/$design
  cp $OUT/pilot56_$design/{run_summary.txt,ef_prediction.csv,temporal_diagnostics.csv,reconstruction_summary.csv,score_correlations.csv} $REPORT/$design/
  head -51 $OUT/pilot56_$design/product_fpcs.csv > $REPORT/$design/product_fpcs_top50.csv
done
cp $OUT/pilot56_sparse/figures/*.png $REPORT/sparse/
cp $OUT/pilot56_dense/figures/temporal_modes.png $REPORT/dense/

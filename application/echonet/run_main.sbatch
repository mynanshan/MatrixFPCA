#!/bin/bash
#SBATCH --account=def-jiguocao_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --job-name=echonet-main
#SBATCH --output=application/echonet/output/slurm-main-%j.log
# Stage 3: main MatFPCA fits at full resolution (112 x 112), n = 500
# stratified TRAIN subjects, 500 VAL, all TEST; sparse and dense designs.
#   sbatch application/echonet/run_main.sbatch   (after run_prepare.sbatch)
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load python/3.11 r/4.5.0
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
OUT=application/echonet/output
EXPORT=$OUT/export_main112
REPORT=application/echonet/report/main
application/echonet/.venv/bin/python -u application/echonet/export_echonet.py \
  --prep $OUT/prep_similarity --out $EXPORT --pool 1 --n-train 500 --n-val 500 --n-test 0
Rscript application/echonet/echonet_to_rds.R $EXPORT
for design in sparse dense; do
  Rscript application/echonet/echonet_matfpca.R --input=$EXPORT \
    --output=$OUT/main112_$design --design=$design --cache=false
done
mkdir -p $REPORT/sparse $REPORT/dense
Rscript application/echonet/compare_pilot_fits.R $OUT/main112_sparse $OUT/main112_dense \
  $REPORT/sparse_vs_dense.csv
cp $EXPORT/{sample_check.csv,sample_strata.csv,sample_representativeness.png} $REPORT/
for design in sparse dense; do
  cp $OUT/main112_$design/{run_summary.txt,ef_prediction.csv,prediction_selected_k.csv,temporal_diagnostics.csv,reconstruction_summary.csv,score_correlations.csv,eigen_video_interpretation.csv,ef_group_change_agreement.txt} $REPORT/$design/
  head -51 $OUT/main112_$design/product_fpcs.csv > $REPORT/$design/product_fpcs_top50.csv
done
cp $OUT/main112_sparse/figures/*.png $REPORT/sparse/
cp $OUT/main112_dense/figures/temporal_modes.png $REPORT/dense/

#!/bin/bash
#SBATCH --account=def-jiguocao_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --job-name=echonet-nstudy
#SBATCH --output=application/echonet/output/slurm-nstudy-%j.log
# Is n = 500 enough? Sparse-design MatFPCA fits at 56 x 56 on nested
# stratified TRAIN samples (250, 500, 1000, 2000) and on a second, disjoint
# 500 (ranks 501-1000), all compared with the n = 2000 fit. The mean and
# covariance are fitted on n subjects; the EF regression and the PCA
# baselines always use the same 2000 TRAIN subjects, so only the MatFPCA
# estimate changes with n. Same VAL (300) and TEST.
#   sbatch application/echonet/run_sample_size.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load python/3.11 r/4.5.0
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
OUT=application/echonet/output
EXPORT=$OUT/export_n56
application/echonet/.venv/bin/python -u application/echonet/export_echonet.py \
  --prep $OUT/prep_similarity --out $EXPORT --pool 2 --n-train 2000 --n-val 300 --n-test 0
Rscript application/echonet/echonet_to_rds.R $EXPORT
run() {  # label n_fit fit_skip
  Rscript application/echonet/echonet_matfpca.R --input=$EXPORT --output=$OUT/nstudy_$1 \
    --design=sparse --n-fit=$2 --fit-skip=$3 --cache=false
}
run n2000 2000 0
run n1000 1000 0
run n500 500 0
run n500b 500 500
run n250 250 0
mkdir -p application/echonet/report/sample_size
Rscript application/echonet/compare_sample_sizes.R $OUT/nstudy_n2000 \
  application/echonet/report/sample_size/sample_size.csv \
  n250=$OUT/nstudy_n250 n500=$OUT/nstudy_n500 n500b=$OUT/nstudy_n500b \
  n1000=$OUT/nstudy_n1000 n2000=$OUT/nstudy_n2000
cp $EXPORT/sample_representativeness.png application/echonet/report/sample_size/sample_representativeness_n2000.png

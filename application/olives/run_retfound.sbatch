#!/bin/bash
#SBATCH --account=def-jiguocao_gpu
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --job-name=olives-retfound
#SBATCH --output=application/olives/output/retfound/slurm-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/olives-matplotlib"
NEURO_PYTHON="${NEURO_PYTHON:-/home/muye/projects/def-jiguocao/muye/env/neuro/bin/python}"
"$NEURO_PYTHON" -u application/olives/extract_retfound.py
"$NEURO_PYTHON" -u application/olives/check_retfound_clinical.py
"$NEURO_PYTHON" -u application/olives/prepare_retfound_demographics.py
"$NEURO_PYTHON" -u application/olives/prepare_retfound_biomarkers.py
module load r/4.5.0
Rscript application/olives/matfpca_retfound.R \
  application/olives/output/retfound/features \
  application/olives/output/retfound/matfpca
Rscript application/olives/associate_retfound_scores.R \
  application/olives/output/retfound/matfpca/model.rds \
  application/olives/output/retfound/matfpca
Rscript application/olives/continuum_retfound.R \
  application/olives/output/retfound/matfpca/model.rds \
  application/olives/output/retfound/matfpca \
  application/olives/output/retfound/demographics.csv \
  application/olives/output/retfound/visit_biomarkers.csv
Rscript application/olives/pca_retfound_scores.R \
  application/olives/output/retfound/matfpca/model.rds \
  application/olives/output/retfound/matfpca \
  application/olives/output/retfound/demographics.csv \
  application/olives/output/retfound/visit_biomarkers.csv

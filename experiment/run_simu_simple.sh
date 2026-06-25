#!/usr/bin/env bash

set -euo pipefail

# Usage:
#   experiment/run_simu_simple.sh [repetitions] [max_parallel]
#   experiment/run_simu_simple.sh 2 4 --dry-run
#
# Defaults produce 19 configurations x 50 replicates = 950 array tasks,
# throttled to at most 100 concurrently running tasks.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SIM_SCRIPT="${ROOT_DIR}/experiment/simu-simple.R"
RESULTS_DIR="${MATFPCA_RESULTS_DIR:-${ROOT_DIR}/experiment/results}"
MANIFEST_DIR="${RESULTS_DIR}/manifests"
LOG_DIR="${ROOT_DIR}/logs"
R_COMMAND="${R_COMMAND:-Rscript}"

run_array_task() {
  local manifest="$1"
  local task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"
  local line_number=$((task_id + 2))
  local row

  row="$(sed -n "${line_number}p" "${manifest}")"
  if [[ -z "${row}" ]]; then
    echo "No manifest row for array task ${task_id}" >&2
    exit 1
  fi

  local config_id replicate n p sigma distribution dependence dep_param
  IFS=$'\t' read -r \
    config_id replicate n p sigma distribution dependence dep_param <<< "${row}"

  echo "Running config=${config_id}, replicate=${replicate}"
  echo "n=${n}, p=q=${p}, sigma=${sigma}, distribution=${distribution}"
  echo "dependence=${dependence}, dependence_param=${dep_param}"

  cd "${ROOT_DIR}"
  "${R_COMMAND}" "${SIM_SCRIPT}" \
    --n_subjects "${n}" \
    --n_time_points 51 \
    --s_obs_min 2 \
    --s_obs_max 4 \
    --d_row "${p}" \
    --d_col "${p}" \
    --noise_lvl "${sigma}" \
    --rank_row 0 \
    --rank_col 0 \
    --rank_time 5 \
    --nbasis 0 \
    --eigset 2 \
    --scores_distr "${distribution}" \
    --score_dependence "${dependence}" \
    --dependence_param "${dep_param}" \
    --seed "${replicate}" \
    --save 1 \
    --results_dir "${RESULTS_DIR}" \
    --manual 0
}

if [[ "${1:-}" == "--worker" ]]; then
  run_array_task "$2"
  exit 0
fi

repetitions="${1:-50}"
max_parallel="${2:-100}"
dry_run="${3:-}"

if ! [[ "${repetitions}" =~ ^[1-9][0-9]*$ ]]; then
  echo "repetitions must be a positive integer" >&2
  exit 2
fi
if ! [[ "${max_parallel}" =~ ^[1-9][0-9]*$ ]]; then
  echo "max_parallel must be a positive integer" >&2
  exit 2
fi

mkdir -p "${MANIFEST_DIR}" "${LOG_DIR}"
timestamp="$(date +%Y%m%d_%H%M%S)"
manifest="${MANIFEST_DIR}/simu_simple_${timestamp}.tsv"
printf "config_id\treplicate\tn\tp\tsigma\tdistribution\tdependence\tdependence_param\n" \
  > "${manifest}"

declare -A seen_configs
config_count=0

add_config() {
  local n="$1"
  local p="$2"
  local sigma="$3"
  local distribution="$4"
  local dependence="$5"
  local dep_param="$6"
  local key="${n}|${p}|${sigma}|${distribution}|${dependence}|${dep_param}"

  if [[ -n "${seen_configs[${key}]:-}" ]]; then
    return
  fi
  seen_configs["${key}"]=1
  config_count=$((config_count + 1))

  local replicate
  for ((replicate = 1; replicate <= repetitions; replicate++)); do
    printf "%03d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${config_count}" "${replicate}" \
      "${n}" "${p}" "${sigma}" "${distribution}" \
      "${dependence}" "${dep_param}" >> "${manifest}"
  done
}

# Baseline: n=50, p=q=30, sigma=0.5, Gaussian independent scores.
# Each block changes one scientifically meaningful feature, except the final
# dependence-by-dimension block, which checks whether spatial averaging is lost.

# Subject-number series.
for n in 25 50 100; do
  add_config "${n}" 30 0.5 normal independent 0
done

# Spatial-dimension series.
for p in 15 30 60; do
  add_config 50 "${p}" 0.5 normal independent 0
done

# Noise series.
for sigma in 0 0.5 1; do
  add_config 50 30 "${sigma}" normal independent 0
done

# Marginal score-distribution series under independence.
for distribution in normal uniform beta32 t5; do
  add_config 50 30 0.5 "${distribution}" independent 0
done

# Dependence-by-dimension series.
# common_amplitude: finite fourth moments, dense dependence among squares.
# common_signed: one common magnitude with independent coordinate signs.
# heavy_amplitude: Pareto A^2 with index 1.5, finite variance of scores but
#                  infinite fourth moments; this intentionally violates theory.
for p in 15 30 60; do
  add_config 50 "${p}" 0.5 normal common_amplitude 1
  add_config 50 "${p}" 0.5 normal common_signed 0
  add_config 50 "${p}" 0.5 normal heavy_amplitude 1.5
done

task_count=$((config_count * repetitions))
echo "Manifest: ${manifest}"
echo "Configurations: ${config_count}"
echo "Replicates per configuration: ${repetitions}"
echo "Array tasks: ${task_count}"

if [[ "${dry_run}" == "--dry-run" ]]; then
  echo "Dry run requested; no jobs submitted."
  exit 0
fi

sbatch \
  --job-name="matfpca-simple" \
  --output="${LOG_DIR}/simu_simple_%A_%a.out" \
  --time="2:00:00" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task=2 \
  --mem="8G" \
  --array="0-$((task_count - 1))%${max_parallel}" \
  "$0" --worker "${manifest}"

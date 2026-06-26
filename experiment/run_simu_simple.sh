#!/usr/bin/env bash

set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module load r/4.5.0
fi

# Usage:
#   experiment/run_simu_simple.sh [repetitions] [max_parallel]
#   experiment/run_simu_simple.sh 2 4 --dry-run
#
# Defaults produce the configured experiment grid x 50 replicates, split into
# resource classes and throttled to at most 100 concurrently running tasks per
# class unless a smaller class cap is set below.
#
# Fir is sensitive to BLAS-related OOMs; the worker pins FlexiBLAS to IMKL and
# caps numerical backends to the CPUs requested by each resource class.

if [[ "${1:-}" == "--worker" ]]; then
  ROOT_DIR="${3:?project root was not passed to the Slurm worker}"
else
  ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
SIM_SCRIPT="${ROOT_DIR}/experiment/simu-simple.R"
RESULTS_DIR="${MATFPCA_RESULTS_DIR:-${ROOT_DIR}/experiment/results}"
LOG_DIR="${ROOT_DIR}/logs"
MANIFEST_DIR="${LOG_DIR}/manifests"
R_COMMAND="${R_COMMAND:-Rscript}"
TIME_COMMAND="${TIME_COMMAND:-/usr/bin/time}"

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

  # Keep every numerical backend within the one CPU allocated by Slurm.
  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
  export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
  export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
  export BLIS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

  # Compute Canada settings
  export FLEXIBLAS=IMKL

  "${TIME_COMMAND}" -v "${R_COMMAND}" "${SIM_SCRIPT}" \
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

declare -A seen_configs
declare -A manifest_by_class
declare -A class_task_count
config_count=0

DIMENSIONS=(20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200)
RESOURCE_CLASSES=(small medium large)

cpus_for_class() {
  case "$1" in
    small) echo 1 ;;
    medium) echo 2 ;;
    large) echo 4 ;;
    *) echo "unknown resource class: $1" >&2; exit 2 ;;
  esac
}

mem_for_class() {
  case "$1" in
    small) echo "4G" ;;
    medium) echo "12G" ;;
    large) echo "32G" ;;
    *) echo "unknown resource class: $1" >&2; exit 2 ;;
  esac
}

parallel_cap_for_class() {
  case "$1" in
    small) echo 100 ;;
    medium) echo 60 ;;
    large) echo 30 ;;
    *) echo "unknown resource class: $1" >&2; exit 2 ;;
  esac
}

class_for_p() {
  local p="$1"
  if (( p <= 60 )); then
    echo small
  elif (( p <= 120 )); then
    echo medium
  else
    echo large
  fi
}

min_int() {
  if (( $1 < $2 )); then echo "$1"; else echo "$2"; fi
}

init_manifest() {
  local class="$1"
  local manifest="${MANIFEST_DIR}/simu_simple_${timestamp}_${class}.tsv"
  manifest_by_class["${class}"]="${manifest}"
  class_task_count["${class}"]=0
  printf "config_id\treplicate\tn\tp\tsigma\tdistribution\tdependence\tdependence_param\n" \
    > "${manifest}"
}

for class in "${RESOURCE_CLASSES[@]}"; do
  init_manifest "${class}"
done

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

  local class manifest
  class="$(class_for_p "${p}")"
  manifest="${manifest_by_class[${class}]}"

  local replicate
  for ((replicate = 1; replicate <= repetitions; replicate++)); do
    printf "%03d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${config_count}" "${replicate}" \
      "${n}" "${p}" "${sigma}" "${distribution}" \
      "${dependence}" "${dep_param}" >> "${manifest}"
  done
  class_task_count["${class}"]=$((class_task_count["${class}"] + repetitions))
}

# Baseline: n=50, p=q=30, sigma=0.5, Gaussian independent scores.
# Each block changes one scientifically meaningful feature, except the final
# dependence-by-dimension block, which checks whether spatial averaging is lost.

# Subject-number series.
for n in 25 50 100; do
  add_config "${n}" 30 0.5 normal independent 0
done

# Spatial-dimension series.
for p in "${DIMENSIONS[@]}"; do
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
for p in "${DIMENSIONS[@]}"; do
  add_config 50 "${p}" 0.5 normal common_amplitude 1
  add_config 50 "${p}" 0.5 normal common_signed 0
  add_config 50 "${p}" 0.5 normal heavy_amplitude 1.5
done

task_count=$((config_count * repetitions))
echo "Configurations: ${config_count}"
echo "Replicates per configuration: ${repetitions}"
echo "Array tasks: ${task_count}"
for class in "${RESOURCE_CLASSES[@]}"; do
  echo "Manifest (${class}): ${manifest_by_class[${class}]}"
  echo "Array tasks (${class}): ${class_task_count[${class}]}"
done

if [[ "${dry_run}" == "--dry-run" ]]; then
  echo "Dry run requested; no jobs submitted."
  for class in "${RESOURCE_CLASSES[@]}"; do
    rm -f "${manifest_by_class[${class}]}"
  done
  exit 0
fi

submitted_jobs=()
for class in "${RESOURCE_CLASSES[@]}"; do
  class_tasks="${class_task_count[${class}]}"
  if (( class_tasks == 0 )); then
    continue
  fi

  class_parallel="$(min_int "${max_parallel}" "$(parallel_cap_for_class "${class}")")"
  class_cpus="$(cpus_for_class "${class}")"
  class_mem="$(mem_for_class "${class}")"
  manifest="${manifest_by_class[${class}]}"

  job_id="$(sbatch --parsable \
    --job-name="matfpca-simple-${class}" \
    --output="${LOG_DIR}/simu_simple_${class}_%A_%a.out" \
    --time="1:00:00" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${class_cpus}" \
    --mem="${class_mem}" \
    --array="0-$((class_tasks - 1))%${class_parallel}" \
    "$0" --worker "${manifest}" "${ROOT_DIR}")"

  submitted_jobs+=("${job_id}")
  echo "Submitted ${class} Slurm array ${job_id} (${class_tasks} tasks, ${class_cpus} CPUs, ${class_mem})."
done

echo "After completion, audit actual use with:"
echo "  sacct -j $(IFS=,; echo "${submitted_jobs[*]}") --units=M --format=JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS"

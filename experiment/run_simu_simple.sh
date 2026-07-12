#!/usr/bin/env bash

set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module load r/4.5.0
fi

# Usage:
#   experiment/run_simu_simple.sh [max_seed] [min_seed] [max_parallel] [options]
#   experiment/run_simu_simple.sh 100 1 100
#   experiment/run_simu_simple.sh 200 101 100 --dry-run
#   experiment/run_simu_simple.sh 100 1 100 --no-overwrite
#   experiment/run_simu_simple.sh 50 1 100 --eigensystem smooth-highrank
#   experiment/run_simu_simple.sh 50 1 100 --eigensystem ar1
#   experiment/run_simu_simple.sh 50 1 100 --eigensystem both --dry-run
#   experiment/run_simu_simple.sh 50 1 100 --nbasis 25
#
# Defaults produce the configured experiment grid for seeds 1-50, split into
# resource classes and throttled to at most 100 concurrently running tasks per
# class unless a smaller class cap is set below. Seed bounds are inclusive.
#
# Fir is sensitive to BLAS-related OOMs; the worker pins FlexiBLAS to IMKL and
# caps numerical backends to the CPUs requested by each resource class.
#
# Eigensystem choices:
#   ar1             default ordered spatial coordinates with AR(1)
#                   correlation 0.7 and low-rank smooth temporal covariance.
#   smooth-highrank smooth spatial Fourier basis with faster spatial
#                   eigendecay and a higher-rank temporal covariance.
#   both            run ar1 and smooth-highrank for each setting.
#
# By default the runner uses the theory-guided automatic temporal basis rule.
# Set MATFPCA_NBASIS to a positive integer or pass --nbasis to force a fixed
# marginal basis dimension.  Large fixed values are expensive because the
# bivariate temporal smoother has approximately nbasis^2 coefficients.

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
TIME_COMMAND="${TIME_COMMAND:-}"
if [[ -z "${TIME_COMMAND}" && -x /usr/bin/time ]]; then
  TIME_COMMAND="/usr/bin/time"
elif [[ -n "${TIME_COMMAND}" ]] && ! command -v "${TIME_COMMAND}" >/dev/null 2>&1; then
  echo "Warning: TIME_COMMAND=${TIME_COMMAND} is unavailable; running without it." >&2
  TIME_COMMAND=""
fi

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

  local config_id seed n p sigma distribution dependence dep_param nbasis weight_folds
  local eigensystem eigvec_set eigval_set rank_time
  IFS=$'\t' read -r \
    config_id seed n p sigma distribution dependence dep_param \
    nbasis weight_folds eigensystem eigvec_set eigval_set rank_time <<< "${row}"

  echo "Running config=${config_id}, seed=${seed}"
  echo "n=${n}, p=q=${p}, sigma=${sigma}, distribution=${distribution}"
  echo "dependence=${dependence}, dependence_param=${dep_param}"
  echo "nbasis=${nbasis}"
  echo "weight_folds=${weight_folds}"
  echo "eigensystem=${eigensystem}, eigvec_set=${eigvec_set}, eigval_set=${eigval_set}, rank_time=${rank_time}"

  cd "${ROOT_DIR}"

  # Keep every numerical backend within the one CPU allocated by Slurm.
  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
  export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
  export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
  export BLIS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

  # Compute Canada settings
  export FLEXIBLAS=IMKL

  local -a sim_command=(
    "${R_COMMAND}" "${SIM_SCRIPT}"
    --n_subjects "${n}"
    --n_time_points 51
    --s_obs_min 2
    --s_obs_max 4
    --d_row "${p}"
    --d_col "${p}"
    --noise_lvl "${sigma}"
    --rank_row 0
    --rank_col 0
    --rank_time "${rank_time}"
    --nbasis "${nbasis}"
    --weight_folds "${weight_folds}"
    --eigset "${eigvec_set}"
    --eigvalset "${eigval_set}"
    --eigensystem "${eigensystem}"
    --scores_distr "${distribution}"
    --score_dependence "${dependence}"
    --dependence_param "${dep_param}"
    --seed "${seed}"
    --save 1
    --results_dir "${RESULTS_DIR}"
    --manual 0
  )

  if [[ -n "${TIME_COMMAND}" ]]; then
    "${TIME_COMMAND}" -v "${sim_command[@]}"
  else
    "${sim_command[@]}"
  fi
}

if [[ "${1:-}" == "--worker" ]]; then
  run_array_task "$2"
  exit 0
fi

max_seed="${1:-50}"
min_seed="${2:-1}"
max_parallel="${3:-100}"
dry_run=0
no_overwrite=0
eigensystem_arg="${MATFPCA_EIGENSYSTEM:-ar1}"
nbasis_arg="${MATFPCA_NBASIS:-0}"
weight_folds_arg="${MATFPCA_WEIGHT_FOLDS:-10}"

shift_count=0
if (($# >= 1)); then shift_count=1; fi
if (($# >= 2)); then shift_count=2; fi
if (($# >= 3)); then shift_count=3; fi
shift "${shift_count}"

while (($#)); do
  case "$1" in
    --dry-run)
      dry_run=1
      shift
      ;;
    --no-overwrite)
      no_overwrite=1
      shift
      ;;
    --eigensystem)
      eigensystem_arg="${2:?--eigensystem requires a value}"
      shift 2
      ;;
    --eigensystem=*)
      eigensystem_arg="${1#*=}"
      shift
      ;;
    --nbasis)
      nbasis_arg="${2:?--nbasis requires a value}"
      shift 2
      ;;
    --nbasis=*)
      nbasis_arg="${1#*=}"
      shift
      ;;
    --weight-folds)
      weight_folds_arg="${2:?--weight-folds requires a value}"
      shift 2
      ;;
    --weight-folds=*)
      weight_folds_arg="${1#*=}"
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

if ! [[ "${max_seed}" =~ ^[1-9][0-9]*$ ]]; then
  echo "max_seed must be a positive integer" >&2
  exit 2
fi
if ! [[ "${min_seed}" =~ ^[1-9][0-9]*$ ]]; then
  echo "min_seed must be a positive integer" >&2
  exit 2
fi
if ! [[ "${max_parallel}" =~ ^[1-9][0-9]*$ ]]; then
  echo "max_parallel must be a positive integer" >&2
  exit 2
fi
if (( min_seed > max_seed )); then
  echo "min_seed must be less than or equal to max_seed" >&2
  exit 2
fi
if ! [[ "${nbasis_arg}" =~ ^[0-9]+$ ]]; then
  echo "nbasis must be a nonnegative integer" >&2
  exit 2
fi
if ! [[ "${weight_folds_arg}" =~ ^[1-9][0-9]*$ ]]; then
  echo "weight-folds must be a positive integer" >&2
  exit 2
fi

eigensystems_for_run() {
  case "$1" in
    smooth-highrank)
      printf "smooth-highrank\n"
      ;;
    ar1)
      printf "ar1\n"
      ;;
    both)
      printf "ar1\nsmooth-highrank\n"
      ;;
    *)
      echo "Unknown eigensystem: $1" >&2
      echo "Valid choices are: ar1, smooth-highrank, both" >&2
      exit 2
      ;;
  esac
}

eigvec_set_for_eigensystem() {
  case "$1" in
    smooth-highrank) echo 2 ;;
    ar1) echo 2 ;;
    *) echo "unknown eigensystem: $1" >&2; exit 2 ;;
  esac
}

eigval_set_for_eigensystem() {
  case "$1" in
    smooth-highrank) echo 2 ;;
    ar1) echo 1 ;;
    *) echo "unknown eigensystem: $1" >&2; exit 2 ;;
  esac
}

rank_time_for_eigensystem() {
  case "$1" in
    smooth-highrank) echo 9 ;;
    ar1) echo 5 ;;
    *) echo "unknown eigensystem: $1" >&2; exit 2 ;;
  esac
}

mapfile -t ACTIVE_EIGENSYSTEMS < <(eigensystems_for_run "${eigensystem_arg}")

seed_count=$((max_seed - min_seed + 1))

mkdir -p "${MANIFEST_DIR}" "${LOG_DIR}"
timestamp="$(date +%Y%m%d_%H%M%S)"

declare -A seen_configs
declare -A manifest_by_class
declare -A class_task_count
config_count=0
skipped_existing_count=0

DIMENSIONS=(20 50 100 150 200)
DISTRIBUTIONS=(normal uniform beta32 t5)
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

# Cross-fitted refined weights repeat the pilot/weight calculation once per
# fold.  With automatic temporal bases, stress tests took about 3 minutes for
# the slowest small/medium task and 4 minutes for the slowest large task.  The
# defaults below retain a generous buffer for seed and node variability without
# requesting the multi-hour limits needed by the old fixed-25 basis runs.
# Each default can be overridden without editing this runner.
time_for_class() {
  case "$1" in
    small) echo "${MATFPCA_TIME_SMALL:-1:00:00}" ;;
    medium) echo "${MATFPCA_TIME_MEDIUM:-1:00:00}" ;;
    large) echo "${MATFPCA_TIME_LARGE:-2:00:00}" ;;
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
  local manifest="${MANIFEST_DIR}/simu_simple_${timestamp}_s${min_seed}-${max_seed}_${class}.tsv"
  manifest_by_class["${class}"]="${manifest}"
  class_task_count["${class}"]=0
  printf "config_id\tseed\tn\tp\tsigma\tdistribution\tdependence\tdependence_param\tnbasis\tweight_folds\teigensystem\teigvec_set\teigval_set\trank_time\n" \
    > "${manifest}"
}

result_file_for_config() {
  local n="$1"
  local p="$2"
  local sigma="$3"
  local distribution="$4"
  local dependence="$5"
  local dep_param="$6"
  local eigensystem="$7"
  local eigvec_set="$8"
  local eigval_set="$9"
  local rank_time="${10}"
  local seed="${11}"
  local nbasis_label="theory"

  if (( nbasis_arg > 0 )); then
    nbasis_label="${nbasis_arg}"
    if (( nbasis_label < 3 )); then nbasis_label=3; fi
    if (( nbasis_label > 51 )); then nbasis_label=51; fi
  fi

  local sigma_token="${sigma//[^[:alnum:]]/p}"
  local dep_param_token="${dep_param//[^[:alnum:]]/p}"
  local eigensystem_token="${eigensystem//[^[:alnum:]]/p}"
  printf '%s/simu_n%s_p%s_sig%s_d%s_dep%s_dp%s_s2-4_nt51_nb%s_wf%s_cdrefined_r%s-%s-%s_eig%s_evc%s_evl%s_sd%s.csv\n' \
    "${RESULTS_DIR}" "${n}" "${p}" "${sigma_token}" "${distribution}" \
    "${dependence}" "${dep_param_token}" "${nbasis_label}" \
    "${weight_folds_arg}" "${p}" "${p}" "${rank_time}" \
    "${eigensystem_token}" "${eigvec_set}" "${eigval_set}" "${seed}"
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
  local eigensystem eigvec_set eigval_set rank_time key

  for eigensystem in "${ACTIVE_EIGENSYSTEMS[@]}"; do
    eigvec_set="$(eigvec_set_for_eigensystem "${eigensystem}")"
    eigval_set="$(eigval_set_for_eigensystem "${eigensystem}")"
    rank_time="$(rank_time_for_eigensystem "${eigensystem}")"
    key="${n}|${p}|${sigma}|${distribution}|${dependence}|${dep_param}|${nbasis_arg}|${weight_folds_arg}|${eigensystem}"

    if [[ -n "${seen_configs[${key}]:-}" ]]; then
      continue
    fi
    seen_configs["${key}"]=1
    config_count=$((config_count + 1))

    local class manifest
    class="$(class_for_p "${p}")"
    manifest="${manifest_by_class[${class}]}"

    local seed
    for ((seed = min_seed; seed <= max_seed; seed++)); do
      if (( no_overwrite )); then
        local result_file
        result_file="$(result_file_for_config \
          "${n}" "${p}" "${sigma}" "${distribution}" "${dependence}" \
          "${dep_param}" "${eigensystem}" "${eigvec_set}" "${eigval_set}" \
          "${rank_time}" "${seed}")"
        if [[ -f "${result_file}" ]]; then
          skipped_existing_count=$((skipped_existing_count + 1))
          continue
        fi
      fi

      printf "%03d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${config_count}" "${seed}" \
        "${n}" "${p}" "${sigma}" "${distribution}" \
        "${dependence}" "${dep_param}" \
        "${nbasis_arg}" "${weight_folds_arg}" \
        "${eigensystem}" "${eigvec_set}" "${eigval_set}" "${rank_time}" >> "${manifest}"
    done
  done
}

# Baseline: n=50, p=q=30, sigma=0.1, Gaussian independent scores.
# Each block changes one scientifically meaningful feature, except the final
# dependence-by-dimension block, which checks whether spatial averaging is lost.

# Subject-number series.
for n in 25 50 75 100; do
  add_config "${n}" 30 0.1 normal independent 0
done

# Primary study: fully cross spatial dimension and marginal score distribution
# under independent FPC scores.
for p in "${DIMENSIONS[@]}"; do
  for distribution in "${DISTRIBUTIONS[@]}"; do
    add_config 50 "${p}" 0.1 "${distribution}" independent 0
  done
done

# Noise series.
for sigma in 0.1 0.5 1 2; do
  add_config 50 30 "${sigma}" normal independent 0
done

# Extension: dependence-by-dimension series at the Gaussian marginal.
# common_amplitude: finite fourth moments, dense dependence among squares.
# common_signed: one common magnitude with independent coordinate signs.
# heavy_amplitude: Pareto A^2 with index 1.5, finite variance of scores but
#                  infinite fourth moments; this intentionally violates theory.
for p in "${DIMENSIONS[@]}"; do
  add_config 50 "${p}" 0.1 normal common_amplitude 1
  add_config 50 "${p}" 0.1 normal common_signed 0
  add_config 50 "${p}" 0.1 normal heavy_amplitude 1.5
done

task_count=0
for class in "${RESOURCE_CLASSES[@]}"; do
  class_task_count["${class}"]=$(($(wc -l < "${manifest_by_class[${class}]}") - 1))
  task_count=$((task_count + class_task_count["${class}"]))
done
echo "Configurations: ${config_count}"
echo "Eigensystem argument: ${eigensystem_arg}"
echo "Eigensystems: ${ACTIVE_EIGENSYSTEMS[*]}"
echo "Temporal marginal basis dimension: ${nbasis_arg} (0 means automatic theory rule)"
echo "Refined-weight cross-fitting folds: ${weight_folds_arg}"
echo "Seed range: ${min_seed}-${max_seed} (${seed_count} seeds per configuration)"
if (( no_overwrite )); then
  echo "Existing-result policy: skip (${skipped_existing_count} tasks skipped)"
else
  echo "Existing-result policy: overwrite"
fi
echo "Array tasks: ${task_count}"
for class in "${RESOURCE_CLASSES[@]}"; do
  echo "Manifest (${class}): ${manifest_by_class[${class}]}"
  echo "Array tasks (${class}): ${class_task_count[${class}]}"
  echo "Time limit (${class}): $(time_for_class "${class}")"
done

if (( dry_run )); then
  echo "Dry run requested; no jobs submitted."
  for class in "${RESOURCE_CLASSES[@]}"; do
    rm -f "${manifest_by_class[${class}]}"
  done
  exit 0
fi

if (( task_count == 0 )); then
  echo "All requested result files already exist; no jobs submitted."
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
  class_time="$(time_for_class "${class}")"
  manifest="${manifest_by_class[${class}]}"

  job_id="$(sbatch --parsable \
    --job-name="matfpca-simple-${class}" \
    --output="${LOG_DIR}/simu_simple_${class}_%A_%a.out" \
    --time="${class_time}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${class_cpus}" \
    --mem="${class_mem}" \
    --array="0-$((class_tasks - 1))%${class_parallel}" \
    "$0" --worker "${manifest}" "${ROOT_DIR}")"

  submitted_jobs+=("${job_id}")
  echo "Submitted ${class} Slurm array ${job_id} (${class_tasks} tasks, ${class_cpus} CPUs, ${class_mem}, ${class_time})."
done

echo "After completion, audit actual use with:"
echo "  sacct -j $(IFS=,; echo "${submitted_jobs[*]}") --units=M --format=JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS"

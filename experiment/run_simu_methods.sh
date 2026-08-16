#!/usr/bin/env bash

set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module load r/4.5.0
fi

# Compute Canada / Slurm launcher for the competitor study.
#
# Usage:
#   experiment/run_simu_methods.sh [max_seed] [min_seed] [max_parallel] [options]
#
# Positional arguments (defaulted when the script is run with no arguments):
#   max_seed       last simulation seed; default: 100
#   min_seed       first simulation seed; default: 1
#   max_parallel   maximum simultaneous tasks in each Slurm array; default: 40
#                  (a smaller, resource-class-specific cap may also apply)
# When using any option below, pass all three positional arguments first, as in
# the examples. (The launcher's positional parser expects options after them.)
#
# Options:
#   --profile NAME          scenario bundle to run; default: main
#   --method-groups LIST    comma-separated scheduler groups; default: all
#   --dry-run               create and print manifests without calling sbatch
#   --no-overwrite          skip tasks whose status file already contains "ok"
#
# A normal full study needs ONE invocation with the default/main profile. The
# profiles are alternative bundles, not pieces that must all be submitted:
#   pilot    small calibration run; its data settings are repeated by main
#   scaling  dimension-only subset (p=q=4,8,16,32) already included in main
#   main     complete final-comparison scenario bundle
# Thus `experiment/run_simu_methods.sh` alone schedules the main profile, seeds
# 1--100, and every method group. Use pilot/scaling separately only when that
# narrower calibration or diagnostic run is wanted.
#
# A "method group" is the unit placed in one Slurm array task. It keeps a slow
# or timed-out competitor from losing results produced by another task. Most
# groups contain one fitted method; `proposal` contains both MatFPCA variants:
#   proposal -> matfpca_refined, matfpca_identity
#   spline   -> spline_pca              product   -> sparse_product_fpca
#   spaco    -> spaco_minus             happ_pace -> happ_pace
#   happ_face -> happ_face              mfaces    -> mfaces
#   mgsfpca  -> mgsfpca                shi       -> shi_mdfpca
# Omitting --method-groups runs all groups. For example,
# `--method-groups proposal,shi,mfaces` schedules only those three groups for
# every profile scenario for which each group is eligible.
#
# Examples:
#   experiment/run_simu_methods.sh 100 1 40 --profile main
#   experiment/run_simu_methods.sh 3 1 20 --profile pilot --dry-run
#   experiment/run_simu_methods.sh 100 1 40 --profile scaling --no-overwrite
#   experiment/run_simu_methods.sh 100 1 20 --method-groups proposal,shi,mfaces
#
# Profiles:
#   pilot    baseline p=q=4 and 8, intended for resource calibration.
#   scaling  p=q=4,8,16,32 at the baseline data-generating setting.
#   main     dimension, sample size, sparsity, noise, spectrum, score-tail,
#            and score-dependence settings for the final comparison table.
#
# Each array task fits one method group. Consequently, a timeout in mFACEs or
# Shi MD-FPCA cannot erase completed results from faster methods. Eligibility
# is method-specific: only scalable methods are scheduled at p=16 or p=32.

if [[ "${1:-}" == "--worker" ]]; then
  ROOT_DIR="${3:?project root was not passed to the Slurm worker}"
else
  ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

SIM_SCRIPT="${ROOT_DIR}/experiment/simu-methods.R"
RESULTS_DIR="${MATFPCA_METHOD_RESULTS_DIR:-${ROOT_DIR}/experiment/results-methods}"
LOG_DIR="${MATFPCA_METHOD_LOG_DIR:-${ROOT_DIR}/logs/methods}"
MANIFEST_DIR="${LOG_DIR}/manifests"
STATUS_DIR="${RESULTS_DIR}/status"
R_COMMAND="${R_COMMAND:-Rscript}"
TIME_COMMAND="${TIME_COMMAND:-}"
if [[ -z "${TIME_COMMAND}" && -x /usr/bin/time ]]; then
  TIME_COMMAND="/usr/bin/time"
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

  local job_id scenario seed n p q nt s1 s2 sigma distribution dependence dep_param
  local eigensystem eigvec_set eigval_set rank_time method_group methods
  local fit_cores blas_threads timeout_limit cov_ncomp
  IFS=$'\t' read -r \
    job_id scenario seed n p q nt s1 s2 sigma distribution dependence dep_param \
    eigensystem eigvec_set eigval_set rank_time method_group methods \
    fit_cores blas_threads timeout_limit cov_ncomp <<< "${row}"

  echo "job_id=${job_id}"
  echo "scenario=${scenario}, seed=${seed}, method_group=${method_group}"
  echo "n=${n}, p=${p}, q=${q}, nt=${nt}, sampling=${s1}-${s2}, sigma=${sigma}"
  echo "distribution=${distribution}, dependence=${dependence}, dependence_param=${dep_param}"
  echo "eigensystem=${eigensystem}, rank_time=${rank_time}, methods=${methods}"
  echo "fit_cores=${fit_cores}, blas_threads=${blas_threads}, timeout=${timeout_limit}"

  cd "${ROOT_DIR}"
  mkdir -p "${STATUS_DIR}"

  # Future-based HAPP jobs use fit_cores workers and one BLAS thread per worker.
  # Other methods use one R process and may use threaded BLAS/LAPACK.
  export OMP_NUM_THREADS="${blas_threads}"
  export MKL_NUM_THREADS="${blas_threads}"
  export OPENBLAS_NUM_THREADS="${blas_threads}"
  export BLIS_NUM_THREADS="${blas_threads}"
  export VECLIB_MAXIMUM_THREADS="${blas_threads}"
  export FLEXIBLAS=IMKL

  local -a sim_command=(
    "${R_COMMAND}" "${SIM_SCRIPT}"
    --n_subjects "${n}"
    --d_row "${p}"
    --d_col "${q}"
    --n_time_points "${nt}"
    --s_obs_min "${s1}"
    --s_obs_max "${s2}"
    --noise_lvl "${sigma}"
    --rank_row 0
    --rank_col 0
    --rank_time "${rank_time}"
    --eigset "${eigvec_set}"
    --eigvalset "${eigval_set}"
    --eigensystem "${eigensystem}"
    --scores_distr "${distribution}"
    --score_dependence "${dependence}"
    --dependence_param "${dep_param}"
    --methods "${methods}"
    --eval_ncomp 5
    --cov_ncomp "${cov_ncomp}"
    --happ_npc 3
    --mfaces_max_outcomes 100
    --allow_slow_competitors 1
    --shi_ncomp "${cov_ncomp}"
    --spaco_rank "${cov_ncomp}"
    --ncores "${fit_cores}"
    --blas_threads "${blas_threads}"
    --backend multicore
    --weight_folds 5
    --seed "${seed}"
    --save 1
    --save_fits 0
    --results_dir "${RESULTS_DIR}"
    --result_tag "${job_id}"
  )

  local exit_code
  set +e
  if [[ -n "${TIME_COMMAND}" ]]; then
    "${TIME_COMMAND}" -v timeout --signal=TERM --kill-after=5m \
      "${timeout_limit}" "${sim_command[@]}"
    exit_code=$?
  else
    timeout --signal=TERM --kill-after=5m \
      "${timeout_limit}" "${sim_command[@]}"
    exit_code=$?
  fi
  set -e

  if (( exit_code == 124 || exit_code == 137 )); then
    printf "timeout\t%s\n" "${timeout_limit}" > "${STATUS_DIR}/${job_id}.status"
    echo "Task exceeded its ${timeout_limit} method-specific limit." >&2
    exit "${exit_code}"
  fi
  if (( exit_code != 0 )); then
    printf "error\t%d\n" "${exit_code}" > "${STATUS_DIR}/${job_id}.status"
    echo "Task failed with exit code ${exit_code}." >&2
    exit "${exit_code}"
  fi
  printf "ok\n" > "${STATUS_DIR}/${job_id}.status"
  echo "Completed ${job_id}"
}

if [[ "${1:-}" == "--worker" ]]; then
  run_array_task "$2"
  exit 0
fi

# The first three arguments are positional for compatibility with existing
# cluster commands. When omitted, they produce seeds 1--100 with a per-array
# concurrency request of 40.
max_seed="${1:-100}"
min_seed="${2:-1}"
max_parallel="${3:-40}"
shift_count=0
if (($# >= 1)); then shift_count=1; fi
if (($# >= 2)); then shift_count=2; fi
if (($# >= 3)); then shift_count=3; fi
shift "${shift_count}"

# These defaults mean that no options are required for the complete main run.
profile="main"
dry_run=0
no_overwrite=0
method_groups_arg="all"
while (($#)); do
  case "$1" in
    --profile)
      profile="${2:?--profile requires pilot, scaling, or main}"
      shift 2
      ;;
    --profile=*)
      profile="${1#*=}"
      shift
      ;;
    --method-groups)
      method_groups_arg="${2:?--method-groups requires a comma-separated value}"
      shift 2
      ;;
    --method-groups=*)
      method_groups_arg="${1#*=}"
      shift
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    --no-overwrite)
      no_overwrite=1
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

for value_name in max_seed min_seed max_parallel; do
  value="${!value_name}"
  if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${value_name} must be a positive integer" >&2
    exit 2
  fi
done
if (( min_seed > max_seed )); then
  echo "min_seed must be less than or equal to max_seed" >&2
  exit 2
fi
case "${profile}" in
  pilot|scaling|main) ;;
  *) echo "Unknown profile: ${profile}" >&2; exit 2 ;;
esac
if ! command -v timeout >/dev/null 2>&1; then
  echo "GNU timeout is required." >&2
  exit 2
fi

# Groups are scheduling/failure-isolation units, not experiment scenarios.
# Selecting a subset changes which methods are fit, while --profile changes
# which data-generating scenarios they are fit to.
ALL_GROUPS=(proposal spline product spaco happ_pace happ_face mfaces mgsfpca shi)
declare -A GROUP_SELECTED
if [[ "${method_groups_arg}" == "all" ]]; then
  for group in "${ALL_GROUPS[@]}"; do GROUP_SELECTED["${group}"]=1; done
else
  IFS=',' read -r -a requested_groups <<< "${method_groups_arg}"
  for group in "${requested_groups[@]}"; do
    valid=0
    for known in "${ALL_GROUPS[@]}"; do
      if [[ "${group}" == "${known}" ]]; then valid=1; break; fi
    done
    if (( ! valid )); then
      echo "Unknown method group: ${group}. Valid groups: ${ALL_GROUPS[*]}" >&2
      exit 2
    fi
    GROUP_SELECTED["${group}"]=1
  done
fi

methods_for_group() {
  # Expand each scheduler group into the --methods value understood by
  # simu-methods.R. Only `proposal` intentionally fits two methods together.
  case "$1" in
    proposal) echo "matfpca_refined,matfpca_identity" ;;
    spline) echo "spline_pca" ;;
    product) echo "sparse_product_fpca" ;;
    spaco) echo "spaco_minus" ;;
    happ_pace) echo "happ_pace" ;;
    happ_face) echo "happ_face" ;;
    mfaces) echo "mfaces" ;;
    mgsfpca) echo "mgsfpca" ;;
    shi) echo "shi_mdfpca" ;;
    *) echo "unknown method group: $1" >&2; exit 2 ;;
  esac
}

max_dimension_for_group() {
  case "$1" in
    proposal|product|spaco|shi) echo 32 ;;
    spline|happ_pace|mgsfpca) echo 16 ;;
    happ_face|mfaces) echo 8 ;;
    *) echo 0 ;;
  esac
}

resource_for_group() {
  local group="$1" dimension="$2"
  case "${group}" in
    proposal|spline)
      if (( dimension <= 8 )); then echo fast_small; else echo fast; fi
      ;;
    product|spaco|mgsfpca)
      if (( dimension <= 8 )); then echo medium_small; else echo medium; fi
      ;;
    happ_pace)
      if (( dimension <= 8 )); then echo parallel_small; else echo parallel; fi
      ;;
    happ_face)
      echo slow_small
      ;;
    shi)
      if (( dimension <= 8 )); then echo slow_small; else echo slow; fi
      ;;
    mfaces)
      echo mfaces
      ;;
  esac
}

fit_cores_for_group() {
  case "$1" in
    happ_pace|happ_face) echo 8 ;;
    *) echo 1 ;;
  esac
}

blas_threads_for_group() {
  case "$1" in
    happ_pace|happ_face) echo 1 ;;
    proposal|spline|product|spaco|mgsfpca) echo 4 ;;
    shi|mfaces) echo 8 ;;
  esac
}

timeout_for_group() {
  local group="$1" dimension="$2"
  if (( dimension <= 8 )); then
    case "${group}" in
      proposal|spline) echo "${MATFPCA_TIMEOUT_FAST_SMALL:-10m}" ;;
      product|spaco|mgsfpca) echo "${MATFPCA_TIMEOUT_MEDIUM_SMALL:-20m}" ;;
      happ_pace) echo "${MATFPCA_TIMEOUT_PARALLEL_SMALL:-20m}" ;;
      happ_face|shi) echo "${MATFPCA_TIMEOUT_SLOW_SMALL:-20m}" ;;
      mfaces) echo "${MATFPCA_TIMEOUT_MFACES:-20m}" ;;
    esac
    return
  fi
  case "${group}" in
    proposal|spline) echo "${MATFPCA_TIMEOUT_FAST:-50m}" ;;
    product) echo "${MATFPCA_TIMEOUT_PRODUCT:-90m}" ;;
    spaco|mgsfpca|happ_pace) echo "${MATFPCA_TIMEOUT_MEDIUM:-210m}" ;;
    shi) echo "${MATFPCA_TIMEOUT_SHI:-330m}" ;;
  esac
}

RESOURCE_CLASSES=(fast_small fast medium_small medium parallel_small parallel slow_small slow mfaces)
cpus_for_class() {
  case "$1" in
    fast_small|fast|medium_small|medium) echo 4 ;;
    parallel_small|parallel|slow_small|slow|mfaces) echo 8 ;;
  esac
}
mem_for_class() {
  case "$1" in
    fast_small) echo "${MATFPCA_MEM_FAST_SMALL:-2G}" ;;
    fast) echo "${MATFPCA_MEM_FAST:-16G}" ;;
    medium_small) echo "${MATFPCA_MEM_MEDIUM_SMALL:-2G}" ;;
    medium) echo "${MATFPCA_MEM_MEDIUM:-32G}" ;;
    parallel_small) echo "${MATFPCA_MEM_PARALLEL_SMALL:-4G}" ;;
    parallel) echo "${MATFPCA_MEM_PARALLEL:-32G}" ;;
    slow_small) echo "${MATFPCA_MEM_SLOW_SMALL:-6G}" ;;
    slow) echo "${MATFPCA_MEM_SLOW:-48G}" ;;
    mfaces) echo "${MATFPCA_MEM_MFACES:-4G}" ;;
  esac
}
time_for_class() {
  case "$1" in
    fast_small) echo "${MATFPCA_TIME_FAST_SMALL:-0:15:00}" ;;
    fast) echo "${MATFPCA_TIME_FAST:-1:00:00}" ;;
    medium_small) echo "${MATFPCA_TIME_MEDIUM_SMALL:-0:30:00}" ;;
    medium) echo "${MATFPCA_TIME_MEDIUM:-4:00:00}" ;;
    parallel_small) echo "${MATFPCA_TIME_PARALLEL_SMALL:-0:30:00}" ;;
    parallel) echo "${MATFPCA_TIME_PARALLEL:-4:00:00}" ;;
    slow_small) echo "${MATFPCA_TIME_SLOW_SMALL:-0:30:00}" ;;
    slow) echo "${MATFPCA_TIME_SLOW:-6:00:00}" ;;
    mfaces) echo "${MATFPCA_TIME_MFACES:-0:30:00}" ;;
  esac
}
parallel_cap_for_class() {
  case "$1" in
    fast_small|fast) echo 60 ;;
    medium_small|medium) echo 30 ;;
    parallel_small|parallel) echo 20 ;;
    slow_small|slow) echo 12 ;;
    mfaces) echo 6 ;;
  esac
}
min_int() { if (( $1 < $2 )); then echo "$1"; else echo "$2"; fi; }

mkdir -p "${MANIFEST_DIR}" "${LOG_DIR}" "${STATUS_DIR}"
timestamp="$(date +%Y%m%d_%H%M%S)"
declare -A manifest_by_class
declare -A task_count_by_class
for resource in "${RESOURCE_CLASSES[@]}"; do
  manifest="${MANIFEST_DIR}/simu_methods_${profile}_${timestamp}_${resource}.tsv"
  manifest_by_class["${resource}"]="${manifest}"
  task_count_by_class["${resource}"]=0
  printf "job_id\tscenario\tseed\tn\tp\tq\tnt\ts1\ts2\tsigma\tdistribution\tdependence\tdep_param\teigensystem\teigvec_set\teigval_set\trank_time\tmethod_group\tmethods\tfit_cores\tblas_threads\ttimeout_limit\tcov_ncomp\n" > "${manifest}"
done

declare -A SEEN_SCENARIOS
scenario_count=0
skipped_dimension=0
skipped_scenario_policy=0
skipped_existing=0

scenario_allowed_for_group() {
  local scenario="$1" group="$2"
  case "${group}" in
    # The quadratic collection of cross-covariance surfaces makes mFACEs the
    # most restrictive method. Keep only the two baseline dimensions and one
    # denser-visit sensitivity setting.
    mfaces)
      [[ "${scenario}" =~ ^(baseline_p4|baseline_p8|dim_p4|dim_p8|visits_20_40_p8)$ ]]
      ;;
    # FACE is much slower than PACE and showed only a small accuracy gain in
    # pilots, so use it for the baseline, denser visits, and noise sensitivity.
    happ_face)
      [[ "${scenario}" =~ ^(baseline_p4|baseline_p8|dim_p4|dim_p8|visits_20_40_p8|noise_05_p8)$ ]]
      ;;
    *)
      return 0
      ;;
  esac
}

add_scenario() {
  local scenario="$1" n="$2" p="$3" q="$4" nt="$5" s1="$6" s2="$7"
  local sigma="$8" distribution="$9" dependence="${10}" dep_param="${11}"
  local eigensystem="${12}" eigvec_set="${13}" eigval_set="${14}" rank_time="${15}"
  if [[ -n "${SEEN_SCENARIOS[${scenario}]:-}" ]]; then return; fi
  SEEN_SCENARIOS["${scenario}"]=1
  scenario_count=$((scenario_count + 1))

  local seed group max_p resource methods fit_cores blas_threads timeout_limit job_id
  for ((seed = min_seed; seed <= max_seed; seed++)); do
    for group in "${ALL_GROUPS[@]}"; do
      if [[ -z "${GROUP_SELECTED[${group}]:-}" ]]; then continue; fi
      if ! scenario_allowed_for_group "${scenario}" "${group}"; then
        skipped_scenario_policy=$((skipped_scenario_policy + 1))
        continue
      fi
      max_p="$(max_dimension_for_group "${group}")"
      if (( p > max_p || q > max_p )); then
        skipped_dimension=$((skipped_dimension + 1))
        continue
      fi
      job_id="${scenario}_sd${seed}_${group}"
      if (( no_overwrite )) && [[ -s "${STATUS_DIR}/${job_id}.status" ]] &&
          grep -qx "ok" "${STATUS_DIR}/${job_id}.status"; then
        skipped_existing=$((skipped_existing + 1))
        continue
      fi
      resource="$(resource_for_group "${group}" "${p}")"
      methods="$(methods_for_group "${group}")"
      fit_cores="$(fit_cores_for_group "${group}")"
      blas_threads="$(blas_threads_for_group "${group}")"
      timeout_limit="$(timeout_for_group "${group}" "${p}")"
      printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\n" \
        "${job_id}" "${scenario}" "${seed}" "${n}" "${p}" "${q}" "${nt}" \
        "${s1}" "${s2}" "${sigma}" "${distribution}" "${dependence}" \
        "${dep_param}" "${eigensystem}" "${eigvec_set}" "${eigval_set}" \
        "${rank_time}" "${group}" "${methods}" "${fit_cores}" "${blas_threads}" \
        "${timeout_limit}" 8 >> "${manifest_by_class[${resource}]}"
      task_count_by_class["${resource}"]=$((task_count_by_class["${resource}"] + 1))
    done
  done
}

add_baseline_dimensions() {
  local p
  for p in 4 8 16 32; do
    add_scenario "dim_p${p}" 50 "${p}" "${p}" 51 2 4 0.1 normal independent 0 ar1 2 1 5
  done
}

case "${profile}" in
  # Exactly one of these branches is used per invocation. `main` contains all
  # scenarios from `scaling`; `pilot` uses equivalent baseline settings under
  # calibration-oriented scenario names.
  pilot)
    add_scenario baseline_p4 50 4 4 51 2 4 0.1 normal independent 0 ar1 2 1 5
    add_scenario baseline_p8 50 8 8 51 2 4 0.1 normal independent 0 ar1 2 1 5
    ;;
  scaling)
    add_baseline_dimensions
    ;;
  main)
    add_baseline_dimensions
    add_scenario n25_p8 25 8 8 51 2 4 0.1 normal independent 0 ar1 2 1 5
    add_scenario n100_p8 100 8 8 51 2 4 0.1 normal independent 0 ar1 2 1 5
    add_scenario visits_20_40_p8 50 8 8 51 4 8 0.1 normal independent 0 ar1 2 1 5
    add_scenario visits_30_60_p8 50 8 8 51 6 12 0.1 normal independent 0 ar1 2 1 5
    add_scenario noise_05_p8 50 8 8 51 2 4 0.5 normal independent 0 ar1 2 1 5
    add_scenario spectral_p8 50 8 8 51 2 4 0.1 normal independent 0 spectral 2 2 9
    add_scenario t5_p8 50 8 8 51 2 4 0.1 t5 independent 0 ar1 2 1 5
    add_scenario dependent_p8 50 8 8 51 2 4 0.1 normal common_amplitude 1 ar1 2 1 5
    ;;
esac

task_count=0
for resource in "${RESOURCE_CLASSES[@]}"; do
  task_count=$((task_count + task_count_by_class["${resource}"]))
done
echo "Profile: ${profile}"
echo "Scenarios: ${scenario_count}"
echo "Seed range: ${min_seed}-${max_seed}"
echo "Selected method groups: ${!GROUP_SELECTED[*]}"
echo "Array tasks: ${task_count}"
echo "Skipped by method-specific dimension limits: ${skipped_dimension}"
echo "Skipped by method-specific scenario policies: ${skipped_scenario_policy}"
echo "Skipped completed tasks: ${skipped_existing}"
for resource in "${RESOURCE_CLASSES[@]}"; do
  echo "${resource}: ${task_count_by_class[${resource}]} tasks; manifest=${manifest_by_class[${resource}]}"
done

if (( dry_run )); then
  echo "Dry run requested; manifests were retained for inspection and no jobs were submitted."
  exit 0
fi
if (( task_count == 0 )); then
  echo "No pending tasks."
  exit 0
fi
if ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch is unavailable. Use --dry-run outside a Slurm login node." >&2
  exit 2
fi

submitted_jobs=()
for resource in "${RESOURCE_CLASSES[@]}"; do
  class_tasks="${task_count_by_class[${resource}]}"
  if (( class_tasks == 0 )); then continue; fi
  class_parallel="$(min_int "${max_parallel}" "$(parallel_cap_for_class "${resource}")")"
  class_cpus="$(cpus_for_class "${resource}")"
  class_mem="$(mem_for_class "${resource}")"
  class_time="$(time_for_class "${resource}")"
  manifest="${manifest_by_class[${resource}]}"
  job_id="$(sbatch --parsable \
    --job-name="matfpca-methods-${resource}" \
    --output="${LOG_DIR}/simu_methods_${resource}_%A_%a.out" \
    --time="${class_time}" \
    --nodes=1 --ntasks=1 \
    --cpus-per-task="${class_cpus}" \
    --mem="${class_mem}" \
    --array="0-$((class_tasks - 1))%${class_parallel}" \
    "$0" --worker "${manifest}" "${ROOT_DIR}")"
  submitted_jobs+=("${job_id}")
  echo "Submitted ${resource} array ${job_id}: ${class_tasks} tasks, ${class_cpus} CPUs, ${class_mem}, ${class_time}."
done

echo "After completion:"
echo "  sacct -j $(IFS=,; echo "${submitted_jobs[*]}") --units=M --format=JobID%24,State,ExitCode,Elapsed,TotalCPU,AllocCPUS,ReqMem,MaxRSS"
echo "  Rscript experiment/summarize_simu_methods.R --results_dir ${RESULTS_DIR}"

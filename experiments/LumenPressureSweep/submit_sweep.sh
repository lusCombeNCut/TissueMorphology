#!/bin/bash
# =============================================================================
# submit_sweep.sh — Lumen Pressure Sweep
#
# Two-phase Slurm submission for BluePebble HPC (ACRC, University of Bristol)
#
# Phase 1 (build):  Single job — git pull, cmake, make CryptBuddingApp
# Phase 2 (run):    Array job — 8 pressures x 10 replicates = 80 tasks
#                   Only starts after build succeeds (--dependency=afterok)
#
# Sweeps lumen pressure while keeping ECM stiffness fixed at the JSON default.
# Uses the -pressure CLI flag available in CryptBuddingApp.
#
# Usage (always run directly, never with sbatch):
#   ./submit_sweep.sh node2d
#   ./submit_sweep.sh vertex2d
#   ./submit_sweep.sh all2d
#   ./submit_sweep.sh all
#
# Prerequisites:
#   1. Container image at /user/work/$(whoami)/containers/tissuemorphology.sif
#   2. Source code at /user/work/$(whoami)/TissueMorphology (or auto-cloned)
# =============================================================================

# ---------- SBATCH defaults (used by simulation array; overridden for build) --
#SBATCH --account=semt036404
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --array=0-79%20
#SBATCH --output=/user/work/%u/logs/lumen_pressure_sweep/slurm-%A_%a.out
#SBATCH --error=/user/work/%u/logs/lumen_pressure_sweep/slurm-%A_%a.err

# ---------- Parse arguments ----------
MODEL_TYPE="${1:-node2d}"

case "$MODEL_TYPE" in
    all)    MODELS_TO_RUN=(node2d vertex2d node3d vertex3d) ;;
    all2d)  MODELS_TO_RUN=(node2d vertex2d) ;;
    all3d)  MODELS_TO_RUN=(node3d vertex3d) ;;
    node2d|vertex2d|node3d|vertex3d)  MODELS_TO_RUN=("$MODEL_TYPE") ;;
    *)
        echo "ERROR: Unknown model type '$MODEL_TYPE'"
        echo "Usage: ./submit_sweep.sh {node2d|vertex2d|node3d|vertex3d|all|all2d|all3d}"
        exit 1
        ;;
esac

# =============================================================================
# SUBMISSION ORCHESTRATION (runs on login node, not inside Slurm)
# =============================================================================
if [ -z "${SLURM_JOB_ID}" ]; then
    BASE_LOG_DIR="/user/work/$(whoami)/logs/lumen_pressure_sweep"
    mkdir -p "${BASE_LOG_DIR}"

    SWEEP_TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

    echo ""
    echo "=== Lumen Pressure Sweep: ${MODELS_TO_RUN[*]} ==="
    echo "  Timestamp: ${SWEEP_TIMESTAMP}"
    echo ""

    # Phase 1: Build
    BUILD_JOB_ID=$(sbatch --parsable \
        --job-name="Build_LumenSweep" \
        --array=0 \
        --cpus-per-task=4 \
        --time=01:00:00 \
        --export=ALL,SWEEP_PHASE=build,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
        --output="${BASE_LOG_DIR}/build_%j.out" \
        --error="${BASE_LOG_DIR}/build_%j.err" \
        "$0" "node2d")

    echo "  Build job:  ${BUILD_JOB_ID}"

    # Phase 2: Simulation arrays
    for model in "${MODELS_TO_RUN[@]}"; do
        SIM_JOB_ID=$(sbatch --parsable \
            --job-name="Lumen_${model}" \
            --dependency=afterok:${BUILD_JOB_ID} \
            --export=ALL,SWEEP_PHASE=run,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
            "$0" "${model}")

        echo "  ${model} sim:   ${SIM_JOB_ID}  (depends on build ${BUILD_JOB_ID})"
        echo "    Output: /user/work/$(whoami)/sim_output/${SIM_JOB_ID}_${model}_${SWEEP_TIMESTAMP}/"
    done

    echo ""
    exit 0
fi

# =============================================================================
# RUNNING INSIDE SLURM
# =============================================================================
SWEEP_PHASE="${SWEEP_PHASE:-run}"

module load apptainer
module load git

export APPTAINER_CACHEDIR="/user/work/$(whoami)/.apptainer"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

SIF_IMAGE="/user/work/$(whoami)/containers/tissuemorphology.sif"
SOURCE_DIR="/user/work/$(whoami)/TissueMorphology"
BUILD_DIR="/user/work/$(whoami)/chaste_build"
BASE_LOG_DIR="/user/work/$(whoami)/logs/lumen_pressure_sweep"
SWEEP_TIMESTAMP="${SWEEP_TIMESTAMP:-$(date +"%Y-%m-%d_%H-%M-%S")}"

mkdir -p "${BASE_LOG_DIR}"

if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container image not found at ${SIF_IMAGE}"
    exit 1
fi

# #############################################################################
# PHASE: BUILD
# #############################################################################
if [ "${SWEEP_PHASE}" = "build" ]; then

    BUILD_LOG_DIR="${BASE_LOG_DIR}/${SLURM_JOB_ID}_build_${SWEEP_TIMESTAMP}"
    mkdir -p "${BUILD_LOG_DIR}"
    LOG_FILE="${BUILD_LOG_DIR}/build.log"
    exec > >(tee -a "${LOG_FILE}")
    exec 2>&1

    echo "============================================"
    echo "  CryptBuddingApp Build (Lumen Pressure Sweep)"
    echo "  Job ID:       ${SLURM_JOB_ID}"
    echo "  Start Time:   $(date)"
    echo "============================================"

    if [ ! -d "${BUILD_DIR}" ]; then
        echo "Seeding build directory from container..."
        mkdir -p "${BUILD_DIR}"
        apptainer exec --bind "${BUILD_DIR}:/mnt" "${SIF_IMAGE}" \
            bash -c "cp -a /home/chaste/build/. /mnt"
    fi

    echo "Updating source..."
    if [ -d "${SOURCE_DIR}/.git" ]; then
        cd "${SOURCE_DIR}" && git fetch origin && git reset --hard origin/main
        echo "  Updated to: $(git rev-parse --short HEAD)"
    else
        git clone https://github.com/lusCombeNCut/TissueMorphology.git "${SOURCE_DIR}"
        cd "${SOURCE_DIR}"
    fi

    rm -rf "${BUILD_DIR}/_deps/cellml_repo-subbuild"

    echo "Building CryptBuddingApp..."
    apptainer exec \
        --bind "${BUILD_DIR}:/home/chaste/build" \
        --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
        --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
        "${SIF_IMAGE}" \
        bash -c "cd /home/chaste/build && \
                 cmake /home/chaste/src -DChaste_ERROR_ON_WARNING=OFF && \
                 make -j${SLURM_CPUS_PER_TASK} CryptBuddingApp"

    BUILD_RC=$?
    if [ ${BUILD_RC} -ne 0 ]; then
        echo "ERROR: Build failed (exit code ${BUILD_RC})"
        exit ${BUILD_RC}
    fi

    HOST_APP_PATH="${BUILD_DIR}/projects/TissueMorphology/apps/CryptBuddingApp"
    if [ ! -f "${HOST_APP_PATH}" ]; then
        echo "ERROR: Binary not found at ${HOST_APP_PATH}"
        exit 1
    fi

    echo "Build successful: ${HOST_APP_PATH}"
    exit 0
fi

# #############################################################################
# PHASE: RUN
# #############################################################################

# ---------- Lumen pressure sweep parameters ----------
PRESSURE_VALUES=(0.0 0.5 1.0 2.0 3.0 5.0 7.0 10.0)
NUM_PRESSURES=${#PRESSURE_VALUES[@]}
NUM_REPLICATES=10

PRESSURE_INDEX=$((SLURM_ARRAY_TASK_ID / NUM_REPLICATES))
RUN_NUMBER=$((SLURM_ARRAY_TASK_ID % NUM_REPLICATES))

if [ "$PRESSURE_INDEX" -ge "$NUM_PRESSURES" ]; then
    echo "Array task $SLURM_ARRAY_TASK_ID exceeds parameter space. Exiting."
    exit 0
fi

LUMEN_PRESSURE=${PRESSURE_VALUES[$PRESSURE_INDEX]}

# ---------- Dynamic job name ----------
case "$MODEL_TYPE" in
    node2d)   JOB_PREFIX="N2" ;;
    vertex2d) JOB_PREFIX="V2" ;;
    node3d)   JOB_PREFIX="N3" ;;
    vertex3d) JOB_PREFIX="V3" ;;
esac
JOB_NAME="${JOB_PREFIX}P${LUMEN_PRESSURE}R${RUN_NUMBER}"
scontrol update JobId="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" JobName="${JOB_NAME}"

# ---------- Paths ----------
RUN_TAG="${SLURM_ARRAY_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}"
LOG_DIR="${BASE_LOG_DIR}/${RUN_TAG}"
OUTPUT_BASE="/user/work/$(whoami)/sim_output/${RUN_TAG}"
OUTPUT_DIR="${OUTPUT_BASE}/p${LUMEN_PRESSURE}_r${RUN_NUMBER}"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/p${LUMEN_PRESSURE}_r${RUN_NUMBER}.log"
exec > >(tee -a "${LOG_FILE}")
exec 2>&1

# ---------- Resolve config ----------
case "$MODEL_TYPE" in
    node2d)   CONFIG_FILENAME="params-Node2d.json"   ;;
    vertex2d) CONFIG_FILENAME="params-Vertex2d.json" ;;
    node3d)   CONFIG_FILENAME="params-Node3d.json"   ;;
    vertex3d) CONFIG_FILENAME="params-Vertex3d.json" ;;
esac
HOST_CONFIG_PATH="${SOURCE_DIR}/apps/${CONFIG_FILENAME}"
CONTAINER_CONFIG_PATH="/home/chaste/src/projects/TissueMorphology/apps/${CONFIG_FILENAME}"

if [ ! -f "${HOST_CONFIG_PATH}" ]; then
    echo "ERROR: Config file not found at ${HOST_CONFIG_PATH}"
    exit 1
fi

echo "============================================"
echo "  Lumen Pressure Sweep"
echo "  Model:           ${MODEL_TYPE}"
echo "  Lumen Pressure:  ${LUMEN_PRESSURE}"
echo "  Replicate:       ${RUN_NUMBER}"
echo "  Start Time:      $(date)"
echo "============================================"

APP_PATH="/home/chaste/build/projects/TissueMorphology/apps/CryptBuddingApp"
HOST_APP_PATH="${BUILD_DIR}/projects/TissueMorphology/apps/CryptBuddingApp"

if [ ! -f "${HOST_APP_PATH}" ]; then
    echo "ERROR: Binary not found at ${HOST_APP_PATH}"
    exit 1
fi

# ---------- Run simulation ----------
apptainer exec \
    --bind "${BUILD_DIR}:/home/chaste/build" \
    --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
    --bind "${OUTPUT_DIR}:/home/chaste/output" \
    --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
    --env CHASTE_TEST_OUTPUT=/home/chaste/output \
    "${SIF_IMAGE}" \
    "${APP_PATH}" \
        -config "${CONTAINER_CONFIG_PATH}" \
        -pressure "${LUMEN_PRESSURE}" \
        -run "${RUN_NUMBER}"

EXIT_CODE=$?

# ---------- Archive ----------
ARCHIVE_DIR="${OUTPUT_BASE}/archives"
mkdir -p "${ARCHIVE_DIR}"
ARCHIVE_PATH="${ARCHIVE_DIR}/p${LUMEN_PRESSURE}_r${RUN_NUMBER}.tar.gz"

if [ -d "${OUTPUT_DIR}" ] && [ "$(ls -A ${OUTPUT_DIR})" ]; then
    tar -czf "${ARCHIVE_PATH}" -C "$(dirname ${OUTPUT_DIR})" "$(basename ${OUTPUT_DIR})"
    echo "Archive: ${ARCHIVE_PATH}"
fi

echo "============================================"
echo "  Job Complete (exit code: ${EXIT_CODE})"
echo "  Lumen Pressure:  ${LUMEN_PRESSURE}"
echo "  Replicate:       ${RUN_NUMBER}"
echo "  End Time:        $(date)"
echo "============================================"

exit ${EXIT_CODE}

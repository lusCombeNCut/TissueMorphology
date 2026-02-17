#!/bin/bash
# =============================================================================
# submit_sweep.sh
#
# Two-phase Slurm submission for BluePebble HPC (ACRC, University of Bristol)
#
# Phase 1 (build):  Single job — git pull, cmake, make CryptBuddingApp
# Phase 2 (run):    Array job — 7 stiffness x 10 replicates = 70 tasks
#                   Only starts after build succeeds (--dependency=afterok)
#
# Supports four model types: node2d, vertex2d, node3d, vertex3d
#
# Usage (always run directly, never with sbatch):
#   ./submit_sweep.sh node2d
#   ./submit_sweep.sh vertex2d
#   ./submit_sweep.sh node3d
#   ./submit_sweep.sh vertex3d
#   ./submit_sweep.sh all            # submits build+array for all 4 models
#   ./submit_sweep.sh all2d          # submits build+array for node2d + vertex2d
#   ./submit_sweep.sh all3d          # submits build+array for node3d + vertex3d
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
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --array=0-69%20
#SBATCH --output=/user/work/%u/logs/crypt_budding/slurm-%A_%a.out
#SBATCH --error=/user/work/%u/logs/crypt_budding/slurm-%A_%a.err

# ---------- Parse arguments ----------
MODEL_TYPE="${1:-node2d}"

case "$MODEL_TYPE" in
    all)
        echo "Submitting all four model types..."
        "$0" node2d
        "$0" vertex2d
        "$0" node3d
        "$0" vertex3d
        exit 0
        ;;
    all2d)
        echo "Submitting both 2D models..."
        "$0" node2d
        "$0" vertex2d
        exit 0
        ;;
    all3d)
        echo "Submitting both 3D models..."
        "$0" node3d
        "$0" vertex3d
        exit 0
        ;;
    node2d|vertex2d|node3d|vertex3d)
        ;;
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
    BASE_LOG_DIR="/user/work/$(whoami)/logs/crypt_budding"
    mkdir -p "${BASE_LOG_DIR}"

    # Generate a shared timestamp for this submission
    SWEEP_TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

    echo ""
    echo "=== Submitting ${MODEL_TYPE} sweep ==="
    echo "  Timestamp: ${SWEEP_TIMESTAMP}"

    # Phase 1: single build job (1h, no array)
    BUILD_JOB_ID=$(sbatch --parsable \
        --job-name="Build_${MODEL_TYPE}" \
        --array=0 \
        --time=01:00:00 \
        --export=ALL,SWEEP_PHASE=build,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
        --output="${BASE_LOG_DIR}/build_${MODEL_TYPE}_%j.out" \
        --error="${BASE_LOG_DIR}/build_${MODEL_TYPE}_%j.err" \
        "$0" "${MODEL_TYPE}")

    echo "  Build job:  ${BUILD_JOB_ID}"

    # Phase 2: simulation array (starts only after build succeeds)
    SIM_JOB_ID=$(sbatch --parsable \
        --job-name="Sim_${MODEL_TYPE}" \
        --dependency=afterok:${BUILD_JOB_ID} \
        --export=ALL,SWEEP_PHASE=run,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
        "$0" "${MODEL_TYPE}")

    echo "  Sim array:  ${SIM_JOB_ID}  (depends on build ${BUILD_JOB_ID})"
    echo "  Logs dir:   ${BASE_LOG_DIR}/${SIM_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}/"
    echo "  Output dir: /user/work/$(whoami)/sim_output/${SIM_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}/"
    echo ""
    exit 0
fi

# =============================================================================
# RUNNING INSIDE SLURM — determine phase
# =============================================================================
SWEEP_PHASE="${SWEEP_PHASE:-run}"

# ---------- Load modules ----------
module load apptainer
module load git

export APPTAINER_CACHEDIR="/user/work/$(whoami)/.apptainer"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# ---------- Common paths ----------
SIF_IMAGE="/user/work/$(whoami)/containers/tissuemorphology.sif"
SOURCE_DIR="/user/work/$(whoami)/TissueMorphology"
BUILD_DIR="/user/work/$(whoami)/chaste_build"
BASE_LOG_DIR="/user/work/$(whoami)/logs/crypt_budding"

# Fallback timestamp if not passed from orchestration
SWEEP_TIMESTAMP="${SWEEP_TIMESTAMP:-$(date +"%Y-%m-%d_%H-%M-%S")}"

mkdir -p "${BASE_LOG_DIR}"

# ---------- Verify container ----------
if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container image not found at ${SIF_IMAGE}"
    echo "  Pull with: apptainer pull ${SIF_IMAGE} docker://ghcr.io/luscombencut/tissuemorphology:latest"
    exit 1
fi

# #############################################################################
# PHASE: BUILD  (single job, one node, no concurrency)
# #############################################################################
if [ "${SWEEP_PHASE}" = "build" ]; then

    BUILD_LOG_DIR="${BASE_LOG_DIR}/${SLURM_JOB_ID}_build_${MODEL_TYPE}_${SWEEP_TIMESTAMP}"
    mkdir -p "${BUILD_LOG_DIR}"

    LOG_FILE="${BUILD_LOG_DIR}/build.log"
    exec > >(tee -a "${LOG_FILE}")
    exec 2>&1

    echo "============================================"
    echo "  CryptBuddingApp Build"
    echo "============================================"
    echo "  Job ID:       ${SLURM_JOB_ID}"
    echo "  Model Sweep:  ${MODEL_TYPE}"
    echo "  Node:         $(hostname)"
    echo "  Start Time:   $(date)"
    echo "============================================"

    # ---- Seed build directory from container (first time only) ----
    if [ ! -d "${BUILD_DIR}" ]; then
        echo ""
        echo "Seeding persistent build directory from container..."
        mkdir -p "${BUILD_DIR}"
        apptainer exec \
            --bind "${BUILD_DIR}:/mnt" \
            "${SIF_IMAGE}" \
            bash -c "cp -a /home/chaste/build/. /mnt"
        echo "  Done."
    fi

    # ---- Update source from GitHub ----
    echo ""
    echo "Updating TissueMorphology source..."
    if [ -d "${SOURCE_DIR}/.git" ]; then
        cd "${SOURCE_DIR}"
        git fetch origin
        git reset --hard origin/main
        echo "  Updated to: $(git rev-parse --short HEAD)"
    else
        echo "  Cloning repository..."
        git clone https://github.com/lusCombeNCut/TissueMorphology.git "${SOURCE_DIR}"
        cd "${SOURCE_DIR}"
        echo "  Cloned at: $(git rev-parse --short HEAD)"
    fi

    # ---- Clean stale FetchContent state ----
    rm -rf "${BUILD_DIR}/_deps/cellml_repo-subbuild"

    # ---- Build ----
    echo ""
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
        echo "ERROR: Build failed with exit code ${BUILD_RC}"
        exit ${BUILD_RC}
    fi

    # ---- Verify binary was produced ----
    HOST_APP_PATH="${BUILD_DIR}/projects/TissueMorphology/apps/CryptBuddingApp"
    if [ ! -f "${HOST_APP_PATH}" ]; then
        echo "ERROR: Binary not found at ${HOST_APP_PATH} despite make returning 0"
        ls -la "${BUILD_DIR}/projects/TissueMorphology/apps/" 2>&1 || true
        exit 1
    fi

    echo ""
    echo "============================================"
    echo "  Build Successful"
    echo "  Binary: ${HOST_APP_PATH}"
    echo "  End Time: $(date)"
    echo "============================================"
    exit 0
fi

# #############################################################################
# PHASE: RUN  (array job — one task per (stiffness, replicate) combo)
# #############################################################################

# ---------- Stiffness sweep parameters ----------
STIFFNESS_VALUES=(0.5 1.0 2.0 5.0 10.0 20.0 50.0)
NUM_STIFFNESS=${#STIFFNESS_VALUES[@]}
NUM_REPLICATES=1

# Decode array task ID → (stiffness_index, replicate)
STIFFNESS_INDEX=$((SLURM_ARRAY_TASK_ID / NUM_REPLICATES))
RUN_NUMBER=$((SLURM_ARRAY_TASK_ID % NUM_REPLICATES))

if [ "$STIFFNESS_INDEX" -ge "$NUM_STIFFNESS" ]; then
    echo "Array task $SLURM_ARRAY_TASK_ID exceeds parameter space. Exiting."
    exit 0
fi

ECM_STIFFNESS=${STIFFNESS_VALUES[$STIFFNESS_INDEX]}

# ---------- Dynamic job name ----------
case "$MODEL_TYPE" in
    node2d)   JOB_PREFIX="N2" ;;
    vertex2d) JOB_PREFIX="V2" ;;
    node3d)   JOB_PREFIX="N3" ;;
    vertex3d) JOB_PREFIX="V3" ;;
esac
JOB_NAME="${JOB_PREFIX}S${ECM_STIFFNESS}R${RUN_NUMBER}"
scontrol update JobId="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" JobName="${JOB_NAME}"

# ---------- Paths ----------
# Group logs and output under <jobid>_<model>_<timestamp>/
RUN_TAG="${SLURM_ARRAY_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}"
LOG_DIR="${BASE_LOG_DIR}/${RUN_TAG}"
OUTPUT_BASE="/user/work/$(whoami)/sim_output/${RUN_TAG}"
OUTPUT_DIR="${OUTPUT_BASE}/s${ECM_STIFFNESS}_r${RUN_NUMBER}"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/s${ECM_STIFFNESS}_r${RUN_NUMBER}.log"
exec > >(tee -a "${LOG_FILE}")
exec 2>&1

# ---------- Job info ----------
echo "============================================"
echo "  Crypt Budding Stiffness Sweep (App)"
echo "============================================"
echo "  Job ID:          ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Model:           ${MODEL_TYPE}"
echo "  ECM Stiffness:   ${ECM_STIFFNESS}"
echo "  Replicate:       ${RUN_NUMBER}"
echo "  Node:            $(hostname)"
echo "  Start Time:      $(date)"
echo "  Output Dir:      ${OUTPUT_DIR}"
echo "============================================"

# ---------- Verify binary exists (built by Phase 1) ----------
APP_PATH="/home/chaste/build/projects/TissueMorphology/apps/CryptBuddingApp"
HOST_APP_PATH="${BUILD_DIR}/projects/TissueMorphology/apps/CryptBuddingApp"

if [ ! -f "${HOST_APP_PATH}" ]; then
    echo "ERROR: CryptBuddingApp binary not found at ${HOST_APP_PATH}"
    echo "  The build job may have failed. Check build logs."
    ls -la "${BUILD_DIR}/projects/TissueMorphology/apps/" 2>&1 || echo "  (directory does not exist)"
    exit 1
fi

# ---------- Run simulation ----------
echo ""
echo "Running CryptBuddingApp (model=${MODEL_TYPE}, stiffness=${ECM_STIFFNESS}, run=${RUN_NUMBER})..."

apptainer exec \
    --bind "${BUILD_DIR}:/home/chaste/build" \
    --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
    --bind "${OUTPUT_DIR}:/home/chaste/output" \
    --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
    --env CHASTE_TEST_OUTPUT=/home/chaste/output \
    "${SIF_IMAGE}" \
    "${APP_PATH}" \
        -model "${MODEL_TYPE}" \
        -stiffness "${ECM_STIFFNESS}" \
        -run "${RUN_NUMBER}"

EXIT_CODE=$?

# ---------- Archive output ----------
ARCHIVE_DIR="${OUTPUT_BASE}/archives"
mkdir -p "${ARCHIVE_DIR}"
ARCHIVE_NAME="s${ECM_STIFFNESS}_r${RUN_NUMBER}.tar.gz"
ARCHIVE_PATH="${ARCHIVE_DIR}/${ARCHIVE_NAME}"

if [ -d "${OUTPUT_DIR}" ] && [ "$(ls -A ${OUTPUT_DIR})" ]; then
    tar -czf "${ARCHIVE_PATH}" -C "$(dirname ${OUTPUT_DIR})" "$(basename ${OUTPUT_DIR})"
    echo "Archive: ${ARCHIVE_PATH}"
fi

# ---------- Summary ----------
echo ""
echo "============================================"
echo "  Job Complete"
echo "  Exit Code:       ${EXIT_CODE}"
echo "  End Time:        $(date)"
echo "  Model:           ${MODEL_TYPE}"
echo "  ECM Stiffness:   ${ECM_STIFFNESS}"
echo "  Replicate:       ${RUN_NUMBER}"
echo "  Output:          ${OUTPUT_DIR}"
echo "  Archive:         ${ARCHIVE_PATH}"
echo "  Log:             ${LOG_FILE}"
echo "============================================"

exit ${EXIT_CODE}

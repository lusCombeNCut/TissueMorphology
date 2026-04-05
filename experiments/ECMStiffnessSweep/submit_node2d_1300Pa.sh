#!/bin/bash
# =============================================================================
# submit_node2d_1300Pa.sh
#
# Submit 5 replicate runs of the vertex2d model at a fixed shear modulus
# of 1300 Pa on BluePebble HPC (ACRC, University of Bristol).
#
# Two phases:
#   Phase 1 (build):  Single job — git pull, cmake, make CryptBuddingApp
#   Phase 2 (run):    Array job — 5 replicates (array indices 0-4)
#                     Starts only after build succeeds (--dependency=afterok)
#
# Usage (run directly on login node, never with sbatch):
#   ./submit_node2d_1300Pa.sh
#
# Prerequisites:
#   1. Container image at /user/work/$(whoami)/containers/tissuemorphology.sif
#   2. Source code at /user/work/$(whoami)/TissueMorphology (or auto-cloned)
# =============================================================================

#SBATCH --account=semt036404
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00
#SBATCH --array=0-4
#SBATCH --output=/user/work/%u/logs/node2d_1300Pa/slurm-%A_%a.out
#SBATCH --error=/user/work/%u/logs/node2d_1300Pa/slurm-%A_%a.err

MODEL_TYPE="node2d"
ECM_SHEAR_MODULUS=1300
NUM_REPLICATES=5

# =============================================================================
# SUBMISSION ORCHESTRATION (runs on login node)
# =============================================================================
if [ -z "${SLURM_JOB_ID}" ]; then
    BASE_LOG_DIR="/user/work/$(whoami)/logs/node2d_1300Pa"
    mkdir -p "${BASE_LOG_DIR}"

    SWEEP_TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

    echo ""
    echo "=== Submitting vertex2d at 1300 Pa — 5 replicates ==="
    echo "  Timestamp: ${SWEEP_TIMESTAMP}"
    echo ""

    # Phase 1: Build
    BUILD_JOB_ID=$(sbatch --parsable \
        --job-name="Build_node2d_1300Pa" \
        --array=0 \
        --cpus-per-task=4 \
        --time=01:00:00 \
        --export=ALL,SWEEP_PHASE=build,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
        --output="${BASE_LOG_DIR}/build_%j.out" \
        --error="${BASE_LOG_DIR}/build_%j.err" \
        "$0")

    echo "  Build job: ${BUILD_JOB_ID}"

    # Phase 2: Simulation array (5 replicates)
    SIM_JOB_ID=$(sbatch --parsable \
        --job-name="Sim_node2d_1300Pa" \
        --dependency=afterok:${BUILD_JOB_ID} \
        --export=ALL,SWEEP_PHASE=run,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
        "$0")

    echo "  Sim array: ${SIM_JOB_ID}  (depends on build ${BUILD_JOB_ID})"
    echo "  Output:    /user/work/$(whoami)/sim_output/${SIM_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}/"
    echo ""
    echo "  To download results locally after completion:"
    echo "    rsync -avz $(whoami)@bp1-login.acrc.bris.ac.uk:/user/work/$(whoami)/sim_output/${SIM_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}/archives/ ./results/"
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
BASE_LOG_DIR="/user/work/$(whoami)/logs/node2d_1300Pa"

SWEEP_TIMESTAMP="${SWEEP_TIMESTAMP:-$(date +"%Y-%m-%d_%H-%M-%S")}"

mkdir -p "${BASE_LOG_DIR}"

if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container image not found at ${SIF_IMAGE}"
    echo "  Pull with: apptainer pull ${SIF_IMAGE} docker://ghcr.io/luscombencut/tissuemorphology:latest"
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
    echo "  CryptBuddingApp Build"
    echo "============================================"
    echo "  Job ID:     ${SLURM_JOB_ID}"
    echo "  Node:       $(hostname)"
    echo "  Start Time: $(date)"
    echo "============================================"

    # Seed build directory from container if not present
    if [ ! -d "${BUILD_DIR}" ]; then
        echo "Seeding persistent build directory from container..."
        mkdir -p "${BUILD_DIR}"
        apptainer exec \
            --bind "${BUILD_DIR}:/mnt" \
            "${SIF_IMAGE}" \
            bash -c "cp -a /home/chaste/build/. /mnt"
    fi

    # Update source
    echo "Updating TissueMorphology source..."
    if [ -d "${SOURCE_DIR}/.git" ]; then
        cd "${SOURCE_DIR}"
        git fetch origin
        git reset --hard origin/main
        echo "  Updated to: $(git rev-parse --short HEAD)"
    else
        git clone https://github.com/lusCombeNCut/TissueMorphology.git "${SOURCE_DIR}"
        cd "${SOURCE_DIR}"
        echo "  Cloned at: $(git rev-parse --short HEAD)"
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
        echo "ERROR: Build failed with exit code ${BUILD_RC}"
        exit ${BUILD_RC}
    fi

    HOST_APP_PATH="${BUILD_DIR}/projects/TissueMorphology/apps/CryptBuddingApp"
    if [ ! -f "${HOST_APP_PATH}" ]; then
        echo "ERROR: Binary not found at ${HOST_APP_PATH} despite make returning 0"
        exit 1
    fi

    echo "============================================"
    echo "  Build Successful — $(date)"
    echo "============================================"
    exit 0
fi

# #############################################################################
# PHASE: RUN  (array 0-4, one replicate per task)
# #############################################################################
RUN_NUMBER=${SLURM_ARRAY_TASK_ID}

RUN_TAG="${SLURM_ARRAY_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}"
LOG_DIR="${BASE_LOG_DIR}/${RUN_TAG}"
OUTPUT_BASE="/user/work/$(whoami)/sim_output/${RUN_TAG}"
OUTPUT_DIR="${OUTPUT_BASE}/g${ECM_SHEAR_MODULUS}_r${RUN_NUMBER}"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/g${ECM_SHEAR_MODULUS}_r${RUN_NUMBER}.log"
exec > >(tee -a "${LOG_FILE}")
exec 2>&1

scontrol update JobId="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" \
    JobName="N2G${ECM_SHEAR_MODULUS}R${RUN_NUMBER}"

HOST_APP_PATH="${BUILD_DIR}/projects/TissueMorphology/apps/CryptBuddingApp"
APP_PATH="/home/chaste/build/projects/TissueMorphology/apps/CryptBuddingApp"
HOST_CONFIG_PATH="${SOURCE_DIR}/apps/params-Node2d.json"
CONTAINER_CONFIG_PATH="/home/chaste/src/projects/TissueMorphology/apps/params-Node2d.json"

echo "============================================"
echo "  vertex2d — 1300 Pa — Replicate ${RUN_NUMBER}"
echo "============================================"
echo "  Job ID:     ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Node:       $(hostname)"
echo "  Start Time: $(date)"
echo "  Output Dir: ${OUTPUT_DIR}"
echo "============================================"

if [ ! -f "${HOST_APP_PATH}" ]; then
    echo "ERROR: CryptBuddingApp binary not found at ${HOST_APP_PATH}"
    exit 1
fi

if [ ! -f "${HOST_CONFIG_PATH}" ]; then
    echo "ERROR: Config file not found at ${HOST_CONFIG_PATH}"
    exit 1
fi

apptainer exec \
    --bind "${BUILD_DIR}:/home/chaste/build" \
    --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
    --bind "${OUTPUT_DIR}:/home/chaste/output" \
    --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
    --env CHASTE_TEST_OUTPUT=/home/chaste/output \
    "${SIF_IMAGE}" \
    "${APP_PATH}" \
        -config "${CONTAINER_CONFIG_PATH}" \
        -model "${MODEL_TYPE}" \
        -shear-modulus "${ECM_SHEAR_MODULUS}" \
        -run "${RUN_NUMBER}"

EXIT_CODE=$?

# Archive output
ARCHIVE_DIR="${OUTPUT_BASE}/archives"
mkdir -p "${ARCHIVE_DIR}"
ARCHIVE_NAME="g${ECM_SHEAR_MODULUS}_r${RUN_NUMBER}.tar.gz"
ARCHIVE_PATH="${ARCHIVE_DIR}/${ARCHIVE_NAME}"

if [ -d "${OUTPUT_DIR}" ] && [ "$(ls -A ${OUTPUT_DIR})" ]; then
    tar -czf "${ARCHIVE_PATH}" -C "$(dirname ${OUTPUT_DIR})" "$(basename ${OUTPUT_DIR})"
    echo "Archive: ${ARCHIVE_PATH}"
fi

echo ""
echo "============================================"
echo "  Job Complete"
echo "  Exit Code:      ${EXIT_CODE}"
echo "  End Time:       $(date)"
echo "  Shear Modulus:  ${ECM_SHEAR_MODULUS} Pa"
echo "  Replicate:      ${RUN_NUMBER}"
echo "  Archive:        ${ARCHIVE_PATH}"
echo "============================================"

exit ${EXIT_CODE}

#!/bin/bash
# =============================================================================
# submit_crypt_budding_sweep.sh
#
# Slurm array job script for BluePebble HPC (ACRC, University of Bristol)
# Runs a stiffness sweep of crypt budding simulations:
#   - 7 stiffness levels x 10 replicates = 70 jobs per model type
#   - Each job runs one (stiffness, replicate) combination
#   - Supports both node-based and vertex-based models
#
# Usage:
#   # Node-based model:
#   sbatch submit_crypt_budding_sweep.sh node
#
#   # Vertex-based model:
#   sbatch submit_crypt_budding_sweep.sh vertex
#
#   # Both models (submits two array jobs):
#   ./submit_crypt_budding_sweep.sh both
#
# Prerequisites:
#   1. Container image at /user/work/$(whoami)/containers/tissuemorphology.sif
#   2. Source code at /user/work/$(whoami)/TissueMorphology (or auto-cloned)
# =============================================================================

#SBATCH --job-name=CryptBud
#SBATCH --account=semt036404
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=0-69

# ---------- Parse arguments ----------
MODEL_TYPE="${1:-node}"

if [ "$MODEL_TYPE" = "both" ]; then
    echo "Submitting both node-based and vertex-based sweeps..."
    sbatch --export=ALL "$0" node
    sbatch --export=ALL "$0" vertex
    exit 0
fi

if [ "$MODEL_TYPE" = "node" ]; then
    TEST_NAME="Test2dCryptBuddingNodeBased"
    TEST_METHOD="TestCryptBuddingStiffnessSweep"
elif [ "$MODEL_TYPE" = "vertex" ]; then
    TEST_NAME="Test2dCryptBuddingVertexBased"
    TEST_METHOD="TestVertexCryptBuddingStiffnessSweep"
else
    echo "ERROR: Unknown model type '$MODEL_TYPE'"
    echo "Usage: sbatch submit_crypt_budding_sweep.sh {node|vertex|both}"
    exit 1
fi

# ---------- Load modules ----------
module load apptainer
module load git

export APPTAINER_CACHEDIR=/user/work/$(whoami)/.apptainer

# ---------- Stiffness levels ----------
# 7 stiffness levels, 10 replicates each = 70 array tasks (0-69)
STIFFNESS_VALUES=(0.5 1.0 2.0 5.0 10.0 20.0 50.0)
NUM_STIFFNESS=${#STIFFNESS_VALUES[@]}
NUM_REPLICATES=10

# Decode array task ID → (stiffness_index, replicate)
STIFFNESS_INDEX=$((SLURM_ARRAY_TASK_ID / NUM_REPLICATES))
RUN_NUMBER=$((SLURM_ARRAY_TASK_ID % NUM_REPLICATES))

# Bounds check
if [ "$STIFFNESS_INDEX" -ge "$NUM_STIFFNESS" ]; then
    echo "Array task $SLURM_ARRAY_TASK_ID exceeds parameter space. Exiting."
    exit 0
fi

ECM_STIFFNESS=${STIFFNESS_VALUES[$STIFFNESS_INDEX]}

# ---------- Configuration ----------
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
SIF_IMAGE="/user/work/$(whoami)/containers/tissuemorphology.sif"
SOURCE_DIR="/user/work/$(whoami)/TissueMorphology"
BUILD_DIR="/user/work/$(whoami)/chaste_build"

# Seed build directory if needed
if [ ! -d "${BUILD_DIR}" ]; then
    echo "Seeding persistent build directory from container image..."
    mkdir -p "${BUILD_DIR}"
    apptainer exec \
        --bind "${BUILD_DIR}:/mnt" \
        "${SIF_IMAGE}" \
        bash -c "cp -a /home/chaste/build/. /mnt"
    echo "  Done."
fi

# Output directory
OUTPUT_DIR="/user/work/$(whoami)/chaste_output/CryptBudding_${MODEL_TYPE}_s${ECM_STIFFNESS}_r${RUN_NUMBER}_job${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${OUTPUT_DIR}"

# Logging
LOG_DIR="/user/work/$(whoami)/logs/crypt_budding"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/${MODEL_TYPE}_s${ECM_STIFFNESS}_r${RUN_NUMBER}_${TIMESTAMP}.log"

exec > >(tee -a "${LOG_FILE}")
exec 2>&1

# ---------- Job info ----------
echo "============================================"
echo "  Crypt Budding Stiffness Sweep"
echo "============================================"
echo "  Job ID:          ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Model:           ${MODEL_TYPE}"
echo "  Test:            ${TEST_NAME}"
echo "  ECM Stiffness:   ${ECM_STIFFNESS}"
echo "  Replicate:       ${RUN_NUMBER}"
echo "  Node:            $(hostname)"
echo "  Start Time:      $(date)"
echo "  Output Dir:      ${OUTPUT_DIR}"
echo "============================================"

# Verify container
if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container image not found at ${SIF_IMAGE}"
    exit 1
fi

# ---------- Update source from GitHub ----------
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

# ---------- Build and run ----------
echo ""
echo "Building and running ${TEST_NAME} (stiffness=${ECM_STIFFNESS}, run=${RUN_NUMBER})..."

apptainer exec \
    --bind "${BUILD_DIR}:/home/chaste/build" \
    --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
    --bind "${OUTPUT_DIR}:/home/chaste/output" \
    --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
    --env CHASTE_TEST_OUTPUT=/home/chaste/output \
    --env ECM_STIFFNESS=${ECM_STIFFNESS} \
    --env RUN_NUMBER=${RUN_NUMBER} \
    "${SIF_IMAGE}" \
    bash -c "cd /home/chaste/build && \
             cmake /home/chaste/src -DChaste_ERROR_ON_WARNING=OFF 2>&1 | tail -3 && \
             make -j${SLURM_CPUS_PER_TASK} ${TEST_NAME} && \
             ctest -R ${TEST_NAME} -V --output-on-failure"

EXIT_CODE=$?

# ---------- Archive output ----------
ARCHIVE_DIR="/user/work/$(whoami)/chaste_output/CryptBudding_${MODEL_TYPE}_results"
mkdir -p "${ARCHIVE_DIR}"
ARCHIVE_NAME="s${ECM_STIFFNESS}_r${RUN_NUMBER}.zip"
ARCHIVE_PATH="${ARCHIVE_DIR}/${ARCHIVE_NAME}"

if [ -d "${OUTPUT_DIR}" ]; then
    cd "$(dirname ${OUTPUT_DIR})"
    zip -rq "${ARCHIVE_PATH}" "$(basename ${OUTPUT_DIR})"
    echo "Archive: ${ARCHIVE_PATH}"
fi

# ---------- Summary ----------
echo ""
echo "============================================"
echo "  Job Complete"
echo "  Exit Code:       ${EXIT_CODE}"
echo "  End Time:        $(date)"
echo "  ECM Stiffness:   ${ECM_STIFFNESS}"
echo "  Replicate:       ${RUN_NUMBER}"
echo "  Output:          ${OUTPUT_DIR}"
echo "  Archive:         ${ARCHIVE_PATH}"
echo "  Log:             ${LOG_FILE}"
echo "============================================"

exit ${EXIT_CODE}

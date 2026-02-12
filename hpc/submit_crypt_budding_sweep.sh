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
MODEL_TYPE="${1:-node}"

if [ "$MODEL_TYPE" = "both" ]; then
    echo "Submitting both node-based and vertex-based sweeps..."
    # Ensure log directory exists before sbatch (Slurm needs it for --output/--error)
    mkdir -p "/user/work/$(whoami)/logs/crypt_budding"
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

# OpenMP thread count — must match cpus-per-task (see HPC docs)
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

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

# ---------- Build (serialised across array tasks) ----------
# Multiple array tasks share the same build directory. Use flock to ensure
# only one task runs cmake/make at a time — prevents build corruption.
BUILD_LOCK="${BUILD_DIR}/.build.lock"
touch "${BUILD_LOCK}"

echo ""
echo "Acquiring build lock (may wait for other tasks)..."
(
    flock -x 200   # exclusive lock on fd 200

    echo "  Lock acquired by task ${SLURM_ARRAY_TASK_ID} at $(date)"

    # Clean stale FetchContent state that breaks cmake reconfigure
    # (cellml_repo-subbuild references paths from the original container build)
    rm -rf "${BUILD_DIR}/_deps/cellml_repo-subbuild"

    apptainer exec \
        --bind "${BUILD_DIR}:/home/chaste/build" \
        --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
        --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
        "${SIF_IMAGE}" \
        bash -c "cd /home/chaste/build && \
                 cmake /home/chaste/src -DChaste_ERROR_ON_WARNING=OFF && \
                 make -j${SLURM_CPUS_PER_TASK} ${TEST_NAME}"

    echo "  Build completed by task ${SLURM_ARRAY_TASK_ID} at $(date)"
) 200>"${BUILD_LOCK}"

BUILD_EXIT=$?
if [ ${BUILD_EXIT} -ne 0 ]; then
    echo "ERROR: Build failed with exit code ${BUILD_EXIT}"
    exit ${BUILD_EXIT}
fi

# ---------- Run simulation ----------
echo ""
echo "Running ${TEST_NAME} (stiffness=${ECM_STIFFNESS}, run=${RUN_NUMBER})..."

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
             ctest -R ${TEST_NAME} -V --output-on-failure"

EXIT_CODE=$?

# ---------- Cancel all sibling jobs on failure ----------
if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "*** JOB FAILED (exit code: ${EXIT_CODE}) — cancelling all remaining array tasks ***"
    scancel ${SLURM_ARRAY_JOB_ID}
fi

# ---------- Copy output to shared results directory ----------
# All array tasks write into one shared results tree, structured as:
#   CryptBudding_<model>_results/stiffness_<s>/run_<r>/
# A single zip archive is created by the LAST finishing task.
RESULTS_DIR="/user/work/$(whoami)/chaste_output/CryptBudding_${MODEL_TYPE}_results"
TASK_RESULTS="${RESULTS_DIR}/stiffness_${ECM_STIFFNESS}/run_${RUN_NUMBER}"
mkdir -p "${TASK_RESULTS}"

if [ -d "${OUTPUT_DIR}" ]; then
    cp -a "${OUTPUT_DIR}/." "${TASK_RESULTS}/"
    echo "Results copied to: ${TASK_RESULTS}"
fi

# Mark this task as done
DONE_DIR="${RESULTS_DIR}/.done"
mkdir -p "${DONE_DIR}"
touch "${DONE_DIR}/task_${SLURM_ARRAY_TASK_ID}"

# Count completed tasks — if all 70 are done, create the final archive
NUM_DONE=$(ls "${DONE_DIR}/" 2>/dev/null | wc -l)
TOTAL_TASKS=70
echo "Completed tasks: ${NUM_DONE}/${TOTAL_TASKS}"

if [ "${NUM_DONE}" -ge "${TOTAL_TASKS}" ]; then
    echo ""
    echo "All tasks complete — creating final archive..."
    ARCHIVE_PATH="/user/work/$(whoami)/chaste_output/CryptBudding_${MODEL_TYPE}_$(date +%Y%m%d_%H%M%S).zip"
    cd "$(dirname ${RESULTS_DIR})"
    zip -rq "${ARCHIVE_PATH}" "$(basename ${RESULTS_DIR})" \
        -x "$(basename ${RESULTS_DIR})/.done/*"
    echo "Final archive: ${ARCHIVE_PATH}"
    echo "Size: $(du -h "${ARCHIVE_PATH}" | cut -f1)"
    echo ""
    echo "To download:"
    echo "  scp sv22482@bp1-login.acrc.bris.ac.uk:${ARCHIVE_PATH} ./"
fi

# ---------- Summary ----------
echo ""
echo "============================================"
echo "  Job Complete"
echo "  Exit Code:       ${EXIT_CODE}"
echo "  End Time:        $(date)"
echo "  ECM Stiffness:   ${ECM_STIFFNESS}"
echo "  Replicate:       ${RUN_NUMBER}"
echo "  Output:          ${TASK_RESULTS}"
echo "  Log:             ${LOG_FILE}"
echo "============================================"

exit ${EXIT_CODE}

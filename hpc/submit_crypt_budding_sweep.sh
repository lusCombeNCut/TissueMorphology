#!/bin/bash
# =============================================================================
# submit_crypt_budding_sweep.sh
#
# Two-phase HPC sweep for BluePebble (ACRC, University of Bristol):
#   Phase 1: A single BUILD job compiles the test binary
#   Phase 2: An ARRAY job (70 tasks) runs simulations, depending on build
#
# This avoids shared-library corruption from concurrent builds.
#
# Usage (run from login node, NOT via sbatch directly):
#   ./submit_crypt_budding_sweep.sh node     # node-based model
#   ./submit_crypt_budding_sweep.sh vertex   # vertex-based model
#   ./submit_crypt_budding_sweep.sh both     # both models
#
# Prerequisites:
#   1. Container image at /user/work/$(whoami)/containers/tissuemorphology.sif
#   2. Source code at /user/work/$(whoami)/TissueMorphology (or auto-cloned)
# =============================================================================

# ---------- Determine execution phase ----------
# SWEEP_PHASE is set by sbatch --export; when run from the command line it is
# empty, meaning we are in "submit" mode.
SWEEP_PHASE="${SWEEP_PHASE:-submit}"

# =====================================================================
# SUBMIT PHASE — orchestrate build → run dependency
# =====================================================================
if [ "$SWEEP_PHASE" = "submit" ]; then

    MODEL_TYPE="${1:-node}"

    if [ "$MODEL_TYPE" = "both" ]; then
        echo "Submitting both node-based and vertex-based sweeps..."
        "$0" node
        "$0" vertex
        exit 0
    fi

    if [ "$MODEL_TYPE" != "node" ] && [ "$MODEL_TYPE" != "vertex" ]; then
        echo "ERROR: Unknown model type '$MODEL_TYPE'"
        echo "Usage: ./submit_crypt_budding_sweep.sh {node|vertex|both}"
        exit 1
    fi

    # Ensure log directory exists (Slurm needs it at submission time)
    mkdir -p "/user/work/$(whoami)/logs/crypt_budding"

    SCRIPT_PATH="$(realpath "$0")"

    # Clean stale .done markers from any previous run so archive logic works
    RESULTS_DIR="/user/work/$(whoami)/chaste_output/CryptBudding_${MODEL_TYPE}_results"
    if [ -d "${RESULTS_DIR}/.done" ]; then
        echo "Cleaning stale .done markers from previous run..."
        rm -rf "${RESULTS_DIR}/.done"
    fi

    echo "=== Phase 1: Submitting BUILD job for ${MODEL_TYPE} ==="
    BUILD_JOB_ID=$(sbatch --parsable \
        --export=ALL,SWEEP_PHASE=build,MODEL_TYPE=${MODEL_TYPE} \
        --job-name="CB_build_${MODEL_TYPE}" \
        --account=semt036404 \
        --partition=compute \
        --nodes=1 \
        --ntasks-per-node=1 \
        --cpus-per-task=4 \
        --mem-per-cpu=4G \
        --time=01:00:00 \
        --output="/user/work/%u/logs/crypt_budding/build_${MODEL_TYPE}_%j.out" \
        --error="/user/work/%u/logs/crypt_budding/build_${MODEL_TYPE}_%j.err" \
        "${SCRIPT_PATH}")

    echo "  Build job submitted: ${BUILD_JOB_ID}"

    echo "=== Phase 2: Submitting RUN array job for ${MODEL_TYPE} (depends on build) ==="
    RUN_JOB_ID=$(sbatch --parsable \
        --dependency=afterok:${BUILD_JOB_ID} \
        --export=ALL,SWEEP_PHASE=run,MODEL_TYPE=${MODEL_TYPE} \
        --job-name="CB_run_${MODEL_TYPE}" \
        --account=semt036404 \
        --partition=compute \
        --nodes=1 \
        --ntasks-per-node=1 \
        --cpus-per-task=4 \
        --mem-per-cpu=4G \
        --time=12:00:00 \
        --array=0-69%20 \
        --output="/user/work/%u/logs/crypt_budding/slurm-%A_%a.out" \
        --error="/user/work/%u/logs/crypt_budding/slurm-%A_%a.err" \
        "${SCRIPT_PATH}")

    echo "  Run array job submitted: ${RUN_JOB_ID} (depends on ${BUILD_JOB_ID})"
    echo ""
    echo "Monitor with:  squeue -u $(whoami)"
    echo "Cancel both:   scancel ${BUILD_JOB_ID} ${RUN_JOB_ID}"
    exit 0
fi

# =====================================================================
# COMMON SETUP (for both build & run phases)
# =====================================================================
module load apptainer
module load git

export APPTAINER_CACHEDIR=/user/work/$(whoami)/.apptainer
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

SIF_IMAGE="/user/work/$(whoami)/containers/tissuemorphology.sif"
SOURCE_DIR="/user/work/$(whoami)/TissueMorphology"
BUILD_DIR="/user/work/$(whoami)/chaste_build"

if [ "$MODEL_TYPE" = "node" ]; then
    TEST_NAME="Test2dCryptBuddingNodeBased"
elif [ "$MODEL_TYPE" = "vertex" ]; then
    TEST_NAME="Test2dCryptBuddingVertexBased"
else
    echo "ERROR: MODEL_TYPE not set"; exit 1
fi

# Verify container
if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container image not found at ${SIF_IMAGE}"
    exit 1
fi

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

# Update source from GitHub
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

# =====================================================================
# BUILD PHASE — single job, compile test binary
# =====================================================================
if [ "$SWEEP_PHASE" = "build" ]; then

    LOG_DIR="/user/work/$(whoami)/logs/crypt_budding"
    mkdir -p "${LOG_DIR}"
    LOG_FILE="${LOG_DIR}/build_${MODEL_TYPE}_$(date +%Y-%m-%d_%H-%M-%S).log"
    exec > >(tee -a "${LOG_FILE}")
    exec 2>&1

    echo "============================================"
    echo "  BUILD PHASE"
    echo "============================================"
    echo "  Job ID:    ${SLURM_JOB_ID}"
    echo "  Model:     ${MODEL_TYPE}"
    echo "  Test:      ${TEST_NAME}"
    echo "  Node:      $(hostname)"
    echo "  Start:     $(date)"
    echo "============================================"

    # Clean stale FetchContent state from container-seeded build
    rm -rf "${BUILD_DIR}/_deps/cellml_repo-subbuild"

    echo ""
    echo "Running cmake + make..."
    apptainer exec \
        --bind "${BUILD_DIR}:/home/chaste/build" \
        --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
        --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
        "${SIF_IMAGE}" \
        bash -c "cd /home/chaste/build && \
                 cmake /home/chaste/src -DChaste_ERROR_ON_WARNING=OFF && \
                 make -j${SLURM_CPUS_PER_TASK} ${TEST_NAME}"

    BUILD_EXIT=$?

    echo ""
    echo "Build exit code: ${BUILD_EXIT}"
    echo "Build finished:  $(date)"

    if [ ${BUILD_EXIT} -ne 0 ]; then
        echo "*** BUILD FAILED ***"
    else
        echo "Build succeeded — array job will start automatically."
    fi

    exit ${BUILD_EXIT}
fi

# =====================================================================
# RUN PHASE — array job, each task runs one simulation
# =====================================================================
if [ "$SWEEP_PHASE" = "run" ]; then

    # ---------- Stiffness levels ----------
    STIFFNESS_VALUES=(0.5 1.0 2.0 5.0 10.0 20.0 50.0)
    NUM_STIFFNESS=${#STIFFNESS_VALUES[@]}
    NUM_REPLICATES=10

    STIFFNESS_INDEX=$((SLURM_ARRAY_TASK_ID / NUM_REPLICATES))
    RUN_NUMBER=$((SLURM_ARRAY_TASK_ID % NUM_REPLICATES))

    if [ "$STIFFNESS_INDEX" -ge "$NUM_STIFFNESS" ]; then
        echo "Array task $SLURM_ARRAY_TASK_ID exceeds parameter space. Exiting."
        exit 0
    fi

    ECM_STIFFNESS=${STIFFNESS_VALUES[$STIFFNESS_INDEX]}
    TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

    # Output directory
    OUTPUT_DIR="/user/work/$(whoami)/chaste_output/CryptBudding_${MODEL_TYPE}_s${ECM_STIFFNESS}_r${RUN_NUMBER}_job${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "${OUTPUT_DIR}"

    # Logging
    LOG_DIR="/user/work/$(whoami)/logs/crypt_budding"
    mkdir -p "${LOG_DIR}"
    LOG_FILE="${LOG_DIR}/${MODEL_TYPE}_s${ECM_STIFFNESS}_r${RUN_NUMBER}_${TIMESTAMP}.log"
    exec > >(tee -a "${LOG_FILE}")
    exec 2>&1

    echo "============================================"
    echo "  Crypt Budding Stiffness Sweep — RUN"
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

    # ---------- Run simulation (NO build — binary already compiled) ----------
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

    # ---------- Cancel all sibling tasks on failure ----------
    if [ $EXIT_CODE -ne 0 ]; then
        echo ""
        echo "*** JOB FAILED (exit code: ${EXIT_CODE}) — cancelling all remaining array tasks ***"
        scancel ${SLURM_ARRAY_JOB_ID}
    fi

    # ---------- Copy output to shared results directory ----------
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
fi

echo "ERROR: Unknown SWEEP_PHASE='${SWEEP_PHASE}'"
exit 1

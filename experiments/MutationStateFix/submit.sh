#!/bin/bash
# =============================================================================
# submit.sh — MutationStateFix
#
# Re-run vertex2d and node2d at 300 Pa and 1300 Pa (5 replicates each)
# after fixing CellMutationStatesWriter CellLabel colour override bug.
#
# Phase 1 (build):   Single job — git pull, cmake, make CryptBuddingApp
# Phase 2 (run):     Array job — 2 stiffness × 5 replicates = 10 tasks/model
# Phase 3 (archive): Single job — tar.gz the entire output directory
#
# No analysis phase — download the archive and analyse locally.
#
# Usage (run directly on login node, never with sbatch):
#   ./submit.sh vertex2d
#   ./submit.sh node2d
#   ./submit.sh all2d          # submits both vertex2d and node2d
#
# =============================================================================

# ---------- SBATCH defaults (simulation array) ----------
#SBATCH --account=semt036404
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00
#SBATCH --array=0-9%10
#SBATCH --output=/user/work/%u/logs/mutation_state_fix/slurm-%A_%a.out
#SBATCH --error=/user/work/%u/logs/mutation_state_fix/slurm-%A_%a.err

# ---------- Parse arguments ----------
MODEL_TYPE="${1:-vertex2d}"

case "$MODEL_TYPE" in
    all2d)
        MODELS_TO_RUN=(vertex2d node2d)
        ;;
    vertex2d|node2d)
        MODELS_TO_RUN=("$MODEL_TYPE")
        ;;
    *)
        echo "ERROR: Unknown model type '$MODEL_TYPE'"
        echo "Usage: ./submit.sh {vertex2d|node2d|all2d}"
        exit 1
        ;;
esac

# =============================================================================
# SUBMISSION ORCHESTRATION (login node — not inside Slurm)
# =============================================================================
if [ -z "${SLURM_JOB_ID}" ]; then
    BASE_LOG_DIR="/user/work/$(whoami)/logs/mutation_state_fix"
    mkdir -p "${BASE_LOG_DIR}"

    SWEEP_TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

    echo ""
    echo "=== MutationStateFix: submitting for ${MODELS_TO_RUN[*]} ==="
    echo "  Timestamp:       ${SWEEP_TIMESTAMP}"
    echo "  Stiffness vals:  300, 1300 Pa"
    echo "  Replicates:      5"
    echo ""

    # Phase 1: Build job (shared by all models)
    BUILD_JOB_ID=$(sbatch --parsable \
        --job-name="Build_MutFix" \
        --array=0 \
        --cpus-per-task=4 \
        --time=01:00:00 \
        --export=ALL,SWEEP_PHASE=build,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
        --output="${BASE_LOG_DIR}/build_%j.out" \
        --error="${BASE_LOG_DIR}/build_%j.err" \
        "$0" "vertex2d")

    echo "  Build job:  ${BUILD_JOB_ID}"

    # Phase 2: Simulation arrays
    SIM_JOB_IDS=()
    for model in "${MODELS_TO_RUN[@]}"; do
        SIM_JOB_ID=$(sbatch --parsable \
            --job-name="MutFix_${model}" \
            --dependency=afterok:${BUILD_JOB_ID} \
            --export=ALL,SWEEP_PHASE=run,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
            "$0" "${model}")

        SIM_JOB_IDS+=("${SIM_JOB_ID}_${model}")
        echo "  ${model} sim:   ${SIM_JOB_ID}  (depends on build ${BUILD_JOB_ID})"
        echo "    Output: /user/work/$(whoami)/sim_output/${SIM_JOB_ID}_${model}_${SWEEP_TIMESTAMP}/"
    done

    # Phase 3: Archive job (runs after all sims)
    DEP_IDS=()
    RUN_TAGS_CSV=""
    for entry in "${SIM_JOB_IDS[@]}"; do
        jid="${entry%%_*}"
        DEP_IDS+=("$jid")
        RUN_TAGS_CSV="${RUN_TAGS_CSV:+${RUN_TAGS_CSV}|}${entry}_${SWEEP_TIMESTAMP}"
    done
    DEP_STR=$(IFS=:; echo "${DEP_IDS[*]}")

    ARCHIVE_JOB_ID=$(sbatch --parsable \
        --job-name="Archive_MutFix" \
        --dependency=afterany:${DEP_STR} \
        --array=0 \
        --cpus-per-task=1 \
        --mem-per-cpu=8G \
        --time=01:00:00 \
        --export="ALL,SWEEP_PHASE=archive,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP},ARCHIVE_RUN_TAGS=${RUN_TAGS_CSV}" \
        --output="${BASE_LOG_DIR}/archive_%j.out" \
        --error="${BASE_LOG_DIR}/archive_%j.err" \
        "$0" "vertex2d")

    ARCHIVE_OUTPUT="/user/work/$(whoami)/sim_output/mutation_state_fix_${SWEEP_TIMESTAMP}"

    echo ""
    echo "  Archive:   ${ARCHIVE_JOB_ID}  (depends on sims: ${DEP_STR})"
    echo ""
    echo "  After completion, download with:"
    echo "    scp $(whoami)@bp1-login.acrc.bris.ac.uk:${ARCHIVE_OUTPUT}.tar.gz ."
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
BASE_LOG_DIR="/user/work/$(whoami)/logs/mutation_state_fix"

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
    echo "  CryptBuddingApp Build (MutationStateFix)"
    echo "  Job ID:       ${SLURM_JOB_ID}"
    echo "  Node:         $(hostname)"
    echo "  Start Time:   $(date)"
    echo "============================================"

    # Seed build directory from container if needed
    if [ ! -d "${BUILD_DIR}" ]; then
        echo "Seeding build directory from container..."
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

    # Clean stale FetchContent state
    rm -rf "${BUILD_DIR}/_deps/cellml_repo-subbuild"

    # Build
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
        echo "ERROR: Binary not found at ${HOST_APP_PATH}"
        exit 1
    fi

    echo "Build successful: ${HOST_APP_PATH}"
    echo "End Time: $(date)"
    exit 0
fi

# #############################################################################
# PHASE: ARCHIVE (single job — tar.gz all simulation output)
# #############################################################################
if [ "${SWEEP_PHASE}" = "archive" ]; then

    ARCHIVE_OUTPUT="/user/work/$(whoami)/sim_output/mutation_state_fix_${SWEEP_TIMESTAMP}"
    mkdir -p "${ARCHIVE_OUTPUT}"

    LOG_FILE="${ARCHIVE_OUTPUT}/archive.log"
    exec > >(tee -a "${LOG_FILE}")
    exec 2>&1

    echo "============================================"
    echo "  MutationStateFix — Archive Phase"
    echo "  Start Time: $(date)"
    echo "============================================"

    IFS='|' read -ra RUN_TAGS <<< "${ARCHIVE_RUN_TAGS}"
    SIM_OUTPUT_BASE="/user/work/$(whoami)/sim_output"

    # Copy each model's run dirs into a per-model subdirectory to avoid
    # vertex2d and node2d g{stiff}_r{rep} folders overwriting each other.
    for tag in "${RUN_TAGS[@]}"; do
        # tag format: <jobid>_<model>_<timestamp>
        # extract model name from the middle segment
        MODEL_NAME=$(echo "${tag}" | sed 's/^[0-9]*_//;s/_[0-9][0-9][0-9][0-9]-.*$//')
        SEARCH_DIR="${SIM_OUTPUT_BASE}/${tag}"
        DEST_DIR="${ARCHIVE_OUTPUT}/${MODEL_NAME}"
        if [ -d "${SEARCH_DIR}" ]; then
            echo "  Copying ${tag} → ${MODEL_NAME}/..."
            mkdir -p "${DEST_DIR}"
            cp -a "${SEARCH_DIR}"/g*_r* "${DEST_DIR}/" 2>/dev/null || true
        else
            echo "  WARNING: ${SEARCH_DIR} not found"
        fi
    done

    # Create the archive
    ARCHIVE_PATH="${ARCHIVE_OUTPUT}.tar.gz"
    echo ""
    echo "Creating archive: ${ARCHIVE_PATH}"
    tar -czf "${ARCHIVE_PATH}" -C "$(dirname ${ARCHIVE_OUTPUT})" "$(basename ${ARCHIVE_OUTPUT})"
    ARCHIVE_SIZE=$(du -sh "${ARCHIVE_PATH}" | cut -f1)

    echo ""
    echo "============================================"
    echo "  Archive Complete — $(date)"
    echo "  Archive: ${ARCHIVE_PATH} (${ARCHIVE_SIZE})"
    echo "  Download with:"
    echo "    scp $(whoami)@bp1-login.acrc.bris.ac.uk:${ARCHIVE_PATH} ."
    echo "============================================"
    exit 0
fi

# #############################################################################
# PHASE: RUN (array job — one task per (stiffness, replicate))
# #############################################################################

SHEAR_MODULUS_VALUES=(300 1300)
NUM_SHEAR_MODULUS=${#SHEAR_MODULUS_VALUES[@]}
NUM_REPLICATES=5

# Array 0-9: index / 5 = stiffness (0→300, 1→1300), index % 5 = replicate
SHEAR_MODULUS_INDEX=$((SLURM_ARRAY_TASK_ID / NUM_REPLICATES))
RUN_NUMBER=$((SLURM_ARRAY_TASK_ID % NUM_REPLICATES))

if [ "$SHEAR_MODULUS_INDEX" -ge "$NUM_SHEAR_MODULUS" ]; then
    echo "Array task $SLURM_ARRAY_TASK_ID exceeds parameter space. Exiting."
    exit 0
fi

ECM_SHEAR_MODULUS=${SHEAR_MODULUS_VALUES[$SHEAR_MODULUS_INDEX]}

# Dynamic job name
case "$MODEL_TYPE" in
    node2d)   JOB_PREFIX="N2" ;;
    vertex2d) JOB_PREFIX="V2" ;;
esac
JOB_NAME="${JOB_PREFIX}G${ECM_SHEAR_MODULUS}R${RUN_NUMBER}"
scontrol update JobId="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" JobName="${JOB_NAME}"

# Paths
RUN_TAG="${SLURM_ARRAY_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}"
LOG_DIR="${BASE_LOG_DIR}/${RUN_TAG}"
OUTPUT_BASE="/user/work/$(whoami)/sim_output/${RUN_TAG}"
OUTPUT_DIR="${OUTPUT_BASE}/g${ECM_SHEAR_MODULUS}_r${RUN_NUMBER}"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/g${ECM_SHEAR_MODULUS}_r${RUN_NUMBER}.log"
exec > >(tee -a "${LOG_FILE}")
exec 2>&1

# Config file
case "$MODEL_TYPE" in
    node2d)   CONFIG_FILENAME="params-Node2d.json"   ;;
    vertex2d) CONFIG_FILENAME="params-Vertex2d.json" ;;
esac
HOST_CONFIG_PATH="${SOURCE_DIR}/apps/${CONFIG_FILENAME}"
CONTAINER_CONFIG_PATH="/home/chaste/src/projects/TissueMorphology/apps/${CONFIG_FILENAME}"

if [ ! -f "${HOST_CONFIG_PATH}" ]; then
    echo "ERROR: Config file not found at ${HOST_CONFIG_PATH}"
    exit 1
fi

echo "============================================"
echo "  MutationStateFix — Simulation"
echo "============================================"
echo "  Job ID:          ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Model:           ${MODEL_TYPE}"
echo "  Config:          ${CONFIG_FILENAME}"
echo "  Shear Modulus:   ${ECM_SHEAR_MODULUS} Pa"
echo "  Replicate:       ${RUN_NUMBER}"
echo "  Node:            $(hostname)"
echo "  Start Time:      $(date)"
echo "  Output Dir:      ${OUTPUT_DIR}"
echo "============================================"

# Verify binary
APP_PATH="/home/chaste/build/projects/TissueMorphology/apps/CryptBuddingApp"
HOST_APP_PATH="${BUILD_DIR}/projects/TissueMorphology/apps/CryptBuddingApp"

if [ ! -f "${HOST_APP_PATH}" ]; then
    echo "ERROR: CryptBuddingApp binary not found at ${HOST_APP_PATH}"
    exit 1
fi

# Run simulation
echo ""
echo "Running CryptBuddingApp..."
apptainer exec \
    --bind "${BUILD_DIR}:/home/chaste/build" \
    --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
    --bind "${OUTPUT_DIR}:/home/chaste/output" \
    --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
    --env CHASTE_TEST_OUTPUT=/home/chaste/output \
    "${SIF_IMAGE}" \
    "${APP_PATH}" \
        -config "${CONTAINER_CONFIG_PATH}" \
        -shear-modulus "${ECM_SHEAR_MODULUS}" \
        -run "${RUN_NUMBER}"

EXIT_CODE=$?

echo ""
echo "============================================"
echo "  Job Complete"
echo "  Exit Code:       ${EXIT_CODE}"
echo "  End Time:        $(date)"
echo "  Shear Modulus:   ${ECM_SHEAR_MODULUS} Pa"
echo "  Replicate:       ${RUN_NUMBER}"
echo "  Output:          ${OUTPUT_DIR}"
echo "============================================"

exit ${EXIT_CODE}

#!/bin/bash
# =============================================================================
# submit_sweep.sh
#
# Slurm array job script for BluePebble HPC (ACRC, University of Bristol)
# Runs a stiffness sweep using the CryptBuddingApp executable:
#   - 7 stiffness levels x 10 replicates = 70 jobs per model type
#   - Each job runs one (stiffness, replicate) combination
#   - Supports all four model types: node2d, vertex2d, node3d, vertex3d
#
# Usage:
#   sbatch submit_sweep.sh node2d
#   sbatch submit_sweep.sh vertex2d
#   sbatch submit_sweep.sh node3d
#   sbatch submit_sweep.sh vertex3d
#   ./submit_sweep.sh all            # submits 4 array jobs (one per model)
#   ./submit_sweep.sh all2d          # submits 2 array jobs (node2d + vertex2d)
#   ./submit_sweep.sh all3d          # submits 2 array jobs (node3d + vertex3d)
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
MODEL_TYPE="${1:-node2d}"

case "$MODEL_TYPE" in
    all)
        echo "Submitting all four model types..."
        mkdir -p "/user/work/$(whoami)/logs/crypt_budding"
        sbatch --export=ALL "$0" node2d
        sbatch --export=ALL "$0" vertex2d
        sbatch --export=ALL "$0" node3d
        sbatch --export=ALL "$0" vertex3d
        exit 0
        ;;
    all2d)
        echo "Submitting both 2D models..."
        mkdir -p "/user/work/$(whoami)/logs/crypt_budding"
        sbatch --export=ALL "$0" node2d
        sbatch --export=ALL "$0" vertex2d
        exit 0
        ;;
    all3d)
        echo "Submitting both 3D models..."
        mkdir -p "/user/work/$(whoami)/logs/crypt_budding"
        sbatch --export=ALL "$0" node3d
        sbatch --export=ALL "$0" vertex3d
        exit 0
        ;;
    node2d|vertex2d|node3d|vertex3d)
        ;;
    *)
        echo "ERROR: Unknown model type '$MODEL_TYPE'"
        echo "Usage: sbatch submit_sweep.sh {node2d|vertex2d|node3d|vertex3d|all|all2d|all3d}"
        exit 1
        ;;
esac

# ---------- Load modules ----------
module load apptainer
module load git

export APPTAINER_CACHEDIR="/user/work/$(whoami)/.apptainer"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# ---------- Stiffness sweep parameters ----------
STIFFNESS_VALUES=(0.5 1.0 2.0 5.0 10.0 20.0 50.0)
NUM_STIFFNESS=${#STIFFNESS_VALUES[@]}
NUM_REPLICATES=10

# Decode array task ID → (stiffness_index, replicate)
STIFFNESS_INDEX=$((SLURM_ARRAY_TASK_ID / NUM_REPLICATES))
RUN_NUMBER=$((SLURM_ARRAY_TASK_ID % NUM_REPLICATES))

if [ "$STIFFNESS_INDEX" -ge "$NUM_STIFFNESS" ]; then
    echo "Array task $SLURM_ARRAY_TASK_ID exceeds parameter space. Exiting."
    exit 0
fi

ECM_STIFFNESS=${STIFFNESS_VALUES[$STIFFNESS_INDEX]}

# ---------- Dynamic job name ----------
# e.g. N2S5.0R3 (Node2d, Stiffness=5.0, Replicate=3)
case "$MODEL_TYPE" in
    node2d)   JOB_PREFIX="N2" ;;
    vertex2d) JOB_PREFIX="V2" ;;
    node3d)   JOB_PREFIX="N3" ;;
    vertex3d) JOB_PREFIX="V3" ;;
esac
JOB_NAME="${JOB_PREFIX}S${ECM_STIFFNESS}R${RUN_NUMBER}"
scontrol update JobId="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" JobName="${JOB_NAME}"

# ---------- Paths ----------
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
SIF_IMAGE="/user/work/$(whoami)/containers/tissuemorphology.sif"
SOURCE_DIR="/user/work/$(whoami)/TissueMorphology"
BUILD_DIR="/user/work/$(whoami)/chaste_build"
LOG_DIR="/user/work/$(whoami)/logs/crypt_budding"
OUTPUT_DIR="/user/work/$(whoami)/sim_output/CryptBudding_${MODEL_TYPE}_s${ECM_STIFFNESS}_r${RUN_NUMBER}"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/${MODEL_TYPE}_s${ECM_STIFFNESS}_r${RUN_NUMBER}_${TIMESTAMP}.log"
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

# ---------- Verify container ----------
if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container image not found at ${SIF_IMAGE}"
    echo "  Pull with: apptainer pull ${SIF_IMAGE} docker://ghcr.io/luscombencut/tissuemorphology:latest"
    exit 1
fi

# ---------- Seed build directory ----------
if [ ! -d "${BUILD_DIR}" ]; then
    echo "Seeding persistent build directory from container..."
    mkdir -p "${BUILD_DIR}"
    apptainer exec \
        --bind "${BUILD_DIR}:/mnt" \
        "${SIF_IMAGE}" \
        bash -c "cp -a /home/chaste/build/. /mnt"
    echo "  Done."
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
BUILD_LOCK="${BUILD_DIR}/.build.lock"
touch "${BUILD_LOCK}"

echo ""
echo "Acquiring build lock (may wait for other tasks)..."
(
    flock -x 200

    echo "  Lock acquired by task ${SLURM_ARRAY_TASK_ID} at $(date)"

    # Clean stale FetchContent state
    rm -rf "${BUILD_DIR}/_deps/cellml_repo-subbuild"

    apptainer exec \
        --bind "${BUILD_DIR}:/home/chaste/build" \
        --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
        --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
        "${SIF_IMAGE}" \
        bash -c "cd /home/chaste/build && \
                 cmake /home/chaste/src -DChaste_ERROR_ON_WARNING=OFF && \
                 make -j${SLURM_CPUS_PER_TASK} CryptBuddingApp"

    echo "  Build completed by task ${SLURM_ARRAY_TASK_ID} at $(date)"
) 200>"${BUILD_LOCK}"

BUILD_EXIT=$?
if [ ${BUILD_EXIT} -ne 0 ]; then
    echo "ERROR: Build failed with exit code ${BUILD_EXIT}"
    exit ${BUILD_EXIT}
fi

# ---------- Run simulation (app — no ctest!) ----------
echo ""
echo "Running CryptBuddingApp (model=${MODEL_TYPE}, stiffness=${ECM_STIFFNESS}, run=${RUN_NUMBER})..."

APP_PATH="/home/chaste/build/projects/TissueMorphology/apps/CryptBuddingApp"

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
ARCHIVE_DIR="/user/work/$(whoami)/sim_output/CryptBudding_${MODEL_TYPE}_results"
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

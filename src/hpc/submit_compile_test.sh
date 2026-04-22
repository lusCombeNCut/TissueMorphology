#!/bin/bash
# =============================================================================
# submit_compile_test.sh
#
# Quick compile-only test for a Chaste test on the HPC.
# Builds the test binary without running it — useful for catching
# compilation errors before submitting a full parameter sweep.
#
# Usage:
#   sbatch hpc/submit_compile_test.sh Test2dCryptBuddingNodeBased
#   sbatch hpc/submit_compile_test.sh Test2dCryptBuddingVertexBased
# =============================================================================

#SBATCH --job-name=CompileTest
#SBATCH --account=semt036404
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00

# ---------- Parse arguments ----------
if [ -z "$1" ]; then
    echo "ERROR: Test name not provided"
    echo "Usage: sbatch submit_compile_test.sh <TestName>"
    exit 1
fi

TEST_NAME="$1"

# ---------- Load modules ----------
module load apptainer
module load git

export APPTAINER_CACHEDIR=/user/work/$(whoami)/.apptainer

# ---------- Configuration ----------
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

LOG_DIR="/user/work/$(whoami)/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/compile_${TEST_NAME}_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "${LOG_FILE}")
exec 2>&1

echo "============================================"
echo "  Compile-only test: ${TEST_NAME}"
echo "  Job ID:  ${SLURM_JOB_ID}"
echo "  Node:    $(hostname)"
echo "  Start:   $(date)"
echo "============================================"

# Verify container
if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container image not found at ${SIF_IMAGE}"
    exit 1
fi

# Update source
echo ""
echo "Updating source..."
if [ -d "${SOURCE_DIR}/.git" ]; then
    cd "${SOURCE_DIR}"
    git fetch origin && git reset --hard origin/main
    echo "  Commit: $(git rev-parse --short HEAD)"
else
    git clone https://github.com/lusCombeNCut/TissueMorphology.git "${SOURCE_DIR}"
    cd "${SOURCE_DIR}"
    echo "  Cloned: $(git rev-parse --short HEAD)"
fi

# Clean stale FetchContent state
rm -rf "${BUILD_DIR}/_deps/cellml_repo-subbuild"

# Configure + compile only (no ctest)
echo ""
echo "Configuring and compiling ${TEST_NAME}..."

apptainer exec \
    --bind "${BUILD_DIR}:/home/chaste/build" \
    --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
    --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
    "${SIF_IMAGE}" \
    bash -c "cd /home/chaste/build && \
             cmake /home/chaste/src -DChaste_ERROR_ON_WARNING=OFF && \
             make -j${SLURM_CPUS_PER_TASK} ${TEST_NAME}"

EXIT_CODE=$?

echo ""
echo "============================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "  COMPILATION SUCCESSFUL"
else
    echo "  COMPILATION FAILED (exit code: ${EXIT_CODE})"
fi
echo "  End:     $(date)"
echo "  Log:     ${LOG_FILE}"
echo "============================================"

exit ${EXIT_CODE}

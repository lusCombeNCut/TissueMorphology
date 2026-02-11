#!/bin/bash
# =============================================================================
# submit_test.sh
#
# Slurm job script for BluePebble HPC (ACRC, University of Bristol)
# Runs a specified Chaste test via Apptainer (Singularity) container.
#
# Usage:
#   sbatch submit_test.sh <TestName>
#
# Example:
#   sbatch submit_test.sh Test3dVertexCryptOrganoid
#
# Prerequisites:
#   1. Pull the container image on the login node first:
#      apptainer pull tissuemorphology.sif docker://ghcr.io/luscombencut/tissuemorphology:latest
#   2. Place the .sif file in your home directory or project space
# =============================================================================

#SBATCH --job-name=ChasteTest
#SBATCH --account=semt036404
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# ---------- Parse arguments ----------
if [ -z "$1" ]; then
    echo "ERROR: Test name not provided"
    echo "Usage: sbatch submit_test.sh <TestName>"
    echo "Example: sbatch submit_test.sh Test3dVertexCryptOrganoid"
    exit 1
fi

TEST_NAME="$1"

# ---------- Load modules ----------
module load apptainer

# Set cache directory to scratch space
export APPTAINER_CACHEDIR=/user/work/$(whoami)/.apptainer

# ---------- Configuration ----------
# Generate timestamp
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

# Path to the Apptainer/Singularity image file (in scratch space)
SIF_IMAGE="/user/work/$(whoami)/containers/tissuemorphology.sif"

# Output directory on the host (bind-mounted into the container)
OUTPUT_DIR="/user/work/$(whoami)/chaste_output/${TEST_NAME}_${TIMESTAMP}_job${SLURM_JOB_ID}"
mkdir -p "${OUTPUT_DIR}"

# Create temporary directory for CTest
TEMP_DIR="/user/work/$(whoami)/chaste_output/temp_${TIMESTAMP}_job${SLURM_JOB_ID}"
mkdir -p "${TEMP_DIR}"

# Log directory and files with timestamp
LOG_DIR="/user/work/$(whoami)/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/${TEST_NAME}_${TIMESTAMP}_job${SLURM_JOB_ID}.log"

# Redirect all output to log file while also showing on screen
exec > >(tee -a "${LOG_FILE}")
exec 2>&1

# ---------- Environment ----------
echo "============================================"
echo "  Job ID:        ${SLURM_JOB_ID}"
echo "  Job Name:      ${SLURM_JOB_NAME}"
echo "  Test Name:     ${TEST_NAME}"
echo "  Node:          $(hostname)"
echo "  CPUs:          ${SLURM_CPUS_PER_TASK}"
echo "  Memory:        ${SLURM_MEM_PER_NODE}MB"
echo "  Account:       semt036404"
echo "  Timestamp:     ${TIMESTAMP}"
echo "  Start Time:    $(date)"
echo "  SIF Image:     ${SIF_IMAGE}"
echo "  Output Dir:    ${OUTPUT_DIR}"
echo "  Log File:      ${LOG_FILE}"
echo "============================================"

# Verify the container image exists
if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container image not found at ${SIF_IMAGE}"
    echo "Run on login node first:"
    echo "  apptainer pull tissuemorphology.sif docker://ghcr.io/luscombencut/tissuemorphology:latest"
    exit 1
fi

# ---------- Run the simulation ----------
# --bind: mount host directories into the container
#   - Output directory for test results
#   - Temp directory for CTest temporary files (writable)
# --env: pass OpenMP thread count to match allocated CPUs
apptainer exec \
    --bind "${OUTPUT_DIR}:/home/chaste/output" \
    --bind "${TEMP_DIR}:/home/chaste/build/Testing/Temporary" \
    --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
    --env CHASTE_TEST_OUTPUT=/home/chaste/output \
    "${SIF_IMAGE}" \
    bash -c "cd /home/chaste/build && ctest -R ${TEST_NAME} -V --output-on-failure"

EXIT_CODE=$?

# Prepare archive name (zip format with timestamp)
ARCHIVE_NAME="${TIMESTAMP}.zip"
ARCHIVE_PATH="/user/work/$(whoami)/chaste_output/${ARCHIVE_NAME}"

# ---------- Summary ----------
echo ""
echo "============================================"
echo "  Job Complete"
echo "  Exit Code:     ${EXIT_CODE}"
echo "  End Time:      $(date)"
echo "  Output Dir:    ${OUTPUT_DIR}"
echo "  Log File:      ${LOG_FILE}"
echo "  Archive:       ${ARCHIVE_PATH}"
echo "============================================"

# List output files
if [ -d "${OUTPUT_DIR}" ]; then
    echo ""
    echo "Output files:"
    find "${OUTPUT_DIR}" -type f | head -50
    
    # Compress output directory as zip archive
    echo ""
    echo "Compressing output directory as zip..."
    
    # Change to parent directory and zip from there
    cd "$(dirname ${OUTPUT_DIR})"
    zip -r "${ARCHIVE_PATH}" "$(basename ${OUTPUT_DIR})"
    
    ARCHIVE_SIZE=$(du -h "${ARCHIVE_PATH}" | cut -f1)
    echo "  Created: ${ARCHIVE_PATH}"
    echo "  Size: ${ARCHIVE_SIZE}"
    echo ""
    echo "  To copy to local machine:"
    echo "    scp sv22482@bp1-login.acrc.bris.ac.uk:${ARCHIVE_PATH} ./"
fi

# Clean up temporary directory
echo ""
echo "Cleaning up temporary files..."
rm -rf "${TEMP_DIR}"

exit ${EXIT_CODE}

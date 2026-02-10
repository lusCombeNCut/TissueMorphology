#!/bin/bash
# =============================================================================
# submit_vertex_organoid.sh
#
# Slurm job script for BluePebble HPC (ACRC, University of Bristol)
# Runs Test3dVertexCryptOrganoid via Apptainer (Singularity) container.
#
# Usage:
#   sbatch submit_vertex_organoid.sh
#
# Prerequisites:
#   1. Pull the container image on the login node first:
#      apptainer pull tissuemorphology.sif docker://ghcr.io/luscombencut/tissuemorphology:latest
#   2. Place the .sif file in your home directory or project space
# =============================================================================

#SBATCH --job-name=VertexOrganoid
#SBATCH --account=semt036404
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=organoid_%j.out
#SBATCH --error=organoid_%j.err

# ---------- Load modules ----------
module load apptainer/1.1.9

# Set cache directory to scratch space
export APPTAINER_CACHEDIR=/user/work/$(whoami)/.apptainer

# ---------- Configuration ----------
# Path to the Apptainer/Singularity image file (in scratch space)
SIF_IMAGE="/user/work/$(whoami)/containers/tissuemorphology.sif"

# Output directory on the host (bind-mounted into the container)
OUTPUT_DIR="/user/work/$(whoami)/chaste_output/vertex_organoid_${SLURM_JOB_ID}"
mkdir -p "${OUTPUT_DIR}"

# ---------- Environment ----------
echo "============================================"
echo "  Job ID:        ${SLURM_JOB_ID}"
echo "  Job Name:      ${SLURM_JOB_NAME}"
echo "  Node:          $(hostname)"
echo "  CPUs:          ${SLURM_CPUS_PER_TASK}"
echo "  Memory:        ${SLURM_MEM_PER_NODE}MB"
echo "  Account:       semt036404"
echo "  Start Time:    $(date)"
echo "  SIF Image:     ${SIF_IMAGE}"
echo "  Output Dir:    ${OUTPUT_DIR}"
echo "============================================"

# Verify the container image exists
if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container image not found at ${SIF_IMAGE}"
    echo "Run on login node first:"
    echo "  apptainer pull tissuemorphology.sif docker://ghcr.io/luscombencut/tissuemorphology:latest"
    exit 1
fi

# ---------- Run the simulation ----------
# --bind: mount host output directory into the container's test output path
# --cleanenv: start with a clean environment (avoids host env contamination)
# --env: pass OpenMP thread count to match allocated CPUs
apptainer exec \
    --cleanenv \
    --bind "${OUTPUT_DIR}:/home/chaste/output" \
    --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
    "${SIF_IMAGE}" \
    bash -c "cd /home/chaste/build && ctest -R Test3dVertexCryptOrganoid -V --output-on-failure"

EXIT_CODE=$?

# ---------- Summary ----------
echo ""
echo "============================================"
echo "  Job Complete"
echo "  Exit Code:     ${EXIT_CODE}"
echo "  End Time:      $(date)"
echo "  Output:        ${OUTPUT_DIR}"
echo "============================================"

# List output files
if [ -d "${OUTPUT_DIR}" ]; then
    echo ""
    echo "Output files:"
    find "${OUTPUT_DIR}" -type f | head -50
fi

exit ${EXIT_CODE}

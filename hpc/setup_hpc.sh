#!/bin/bash
# =============================================================================
# setup_hpc.sh
#
# Run this script ONCE on the BluePebble login node to pull the container
# image and set up the working directory.
#
# Usage (on BluePebble login node):
#   bash setup_hpc.sh
# =============================================================================

set -euo pipefail

IMAGE="docker://ghcr.io/luscombencut/tissuemorphology:latest"
CONTAINER_DIR="/user/work/$(whoami)/containers"
SIF_FILE="${CONTAINER_DIR}/tissuemorphology.sif"

echo "============================================"
echo "  BluePebble HPC — Container Setup"
echo "  User: $(whoami)"
echo "  Host: $(hostname)"
echo "============================================"

# Load Apptainer module (required on BluePebble)
echo ""
echo "Loading Apptainer module..."
module load apptainer

# Set cache directory to scratch space (avoids filling home directory)
export APPTAINER_CACHEDIR=/user/work/$(whoami)/.apptainer
mkdir -p "$APPTAINER_CACHEDIR"
mkdir -p "$CONTAINER_DIR"
echo "Cache directory: $APPTAINER_CACHEDIR"
echo "Container directory: $CONTAINER_DIR"

# ---------- Check authentication ----------
echo ""
echo "[1/4] Checking GHCR authentication..."
if ! apptainer remote status docker://ghcr.io &>/dev/null; then
    echo ""
    echo "⚠️  Authentication required for GitHub Container Registry"
    echo ""
    echo "Please run this command first:"
    echo "  apptainer remote login -u lusCombeNCut docker://ghcr.io"
    echo ""
    echo "When prompted for password, use a GitHub Personal Access Token"
    echo "with 'read:packages' permission."
    echo ""
    echo "Create token at: https://github.com/settings/tokens"
    echo ""
    exit 1
fi
echo "✓ Authenticated with GHCR"

# ---------- Pull container image ----------
echo ""
echo "[2/4] Pulling container image..."
echo "  Source: ${IMAGE}"
echo "  Target: ${SIF_FILE}"
echo ""

# If the image exists, ask to replace
if [ -f "${SIF_FILE}" ]; then
    echo "WARNING: ${SIF_FILE} already exists."
    read -p "Replace it? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Skipping pull. Using existing image."
    else
        rm -f "${SIF_FILE}"
        apptainer pull "${SIF_FILE}" "${IMAGE}"
    fi
else
    apptainer pull "${SIF_FILE}" "${IMAGE}"
fi

echo ""
echo "[3/4] Verifying container..."
apptainer exec "${SIF_FILE}" cat /etc/os-release | head -2
echo ""
apptainer exec "${SIF_FILE}" bash -c 'cd /home/chaste/build && ctest -N -R Test3dVertexCryptOrganoid' 2>/dev/null || \
    echo "  (Test listing will work after full build)"

# ---------- Create output and log directories ----------
echo ""
echo "[4/4] Creating output and log directories..."
mkdir -p "/user/work/$(whoami)/chaste_output"
mkdir -p "/user/work/$(whoami)/logs"

echo ""
echo "============================================"
echo "  Setup complete!"
echo ""
echo "  Container: ${SIF_FILE}"
echo "  Output:    /user/work/$(whoami)/chaste_output"
echo "  Logs:      /user/work/$(whoami)/logs"
echo ""
echo "  To submit a job:"
echo "    sbatch hpc/submit_vertex_organoid.sh"
echo ""
echo "  To run interactively (for testing):"
echo "    srun --account=semt036404 --cpus-per-task=4 --mem=16G --time=01:00:00 --pty bash"
echo "    module load apptainer"
echo "    export WORK=/user/work/\$(whoami)"
echo "    apptainer exec --bind \${WORK}/chaste_output:/home/chaste/output \${WORK}/containers/tissuemorphology.sif \\"
echo "        bash -c 'cd /home/chaste/build && ctest -R Test3dVertexCryptOrganoid -V'"
echo "============================================"

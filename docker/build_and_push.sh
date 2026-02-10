#!/bin/bash
# =============================================================================
# build_and_push.sh
#
# Build the Docker image locally and push to GitHub Container Registry.
# Then provides instructions for pulling onto BluePebble HPC.
#
# Usage (from TissueMorphology/ directory):
#   ./docker/build_and_push.sh
#
# Or from docker/ directory:
#   ./build_and_push.sh
#
# Prerequisites:
#   - Docker Desktop running
#   - GitHub Personal Access Token with write:packages permission
#
# To authenticate (run ONCE before this script):
#   In WSL/bash:
#     echo YOUR_TOKEN_HERE | docker login ghcr.io -u lusCombeNCut --password-stdin
#
#   In PowerShell:
#     "YOUR_TOKEN_HERE" | docker login ghcr.io -u lusCombeNCut --password-stdin
#
# Get token at: https://github.com/settings/tokens/new
# =============================================================================

set -euo pipefail

IMAGE_NAME="ghcr.io/luscombencut/tissuemorphology"
IMAGE_TAG="latest"
FULL_IMAGE="${IMAGE_NAME}:${IMAGE_TAG}"

echo "============================================"
echo "  Building Chaste + TissueMorphology image"
echo "  Image: ${FULL_IMAGE}"
echo "============================================"

# ---------- Check authentication ----------
echo ""
echo "[1/4] Checking authentication..."
if docker pull ghcr.io/luscombencut/tissuemorphology:latest >/dev/null 2>&1 || \
   grep -q "ghcr.io" ~/.docker/config.json 2>/dev/null; then
    echo "✓ Already authenticated to ghcr.io"
else
    echo ""
    echo "❌ Not authenticated to GitHub Container Registry"
    echo ""
    echo "Please authenticate first:"
    echo ""
    echo "1. Get your token at: https://github.com/settings/tokens/new"
    echo "   - Name: Docker/Chaste push"
    echo "   - Scopes: ✓ write:packages"
    echo ""
    echo "2. Run this command and paste your token when prompted:"
    echo ""
    if [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "win32" ]]; then
        echo "   # PowerShell:"
        echo "   \"YOUR_TOKEN\" | docker login ghcr.io -u lusCombeNCut --password-stdin"
    else
        echo "   # WSL/bash:"
        echo "   docker login ghcr.io -u lusCombeNCut"
        echo "   (then paste token when prompted)"
    fi
    echo ""
    exit 1
fi

# ---------- Build ----------
# Resolve path: build context = Chaste/projects/ (parent of TissueMorphology)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"          # TissueMorphology/
CONTEXT_DIR="$(dirname "${PROJECT_DIR}")"           # Chaste/projects/

echo ""
echo "[2/4] Building Docker image..."
echo "  Dockerfile: ${SCRIPT_DIR}/Dockerfile"
echo "  Context:    ${CONTEXT_DIR}"
docker build \
    --progress=plain \
    -f "${SCRIPT_DIR}/Dockerfile" \
    -t "${FULL_IMAGE}" \
    "${CONTEXT_DIR}"

echo ""
echo "[3/4] Build complete. Image size:"
docker images "${IMAGE_NAME}" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}"

# ---------- Push ----------
echo ""
echo "[4/4] Pushing to GitHub Container Registry..."
docker push "${FULL_IMAGE}"

echo ""
echo "============================================"
echo "  Push complete!"
echo "============================================"
echo ""
echo "Next steps on BluePebble HPC (ssh sv22482@bp1-login.acrc.bris.ac.uk):"
echo ""
echo "  # 1. Pull the container image (login node only, do NOT run on compute)"
echo "  apptainer pull tissuemorphology.sif docker://${FULL_IMAGE}"
echo ""
echo "  # 2. Submit the Slurm job"
echo "  sbatch hpc/submit_vertex_organoid.sh"
echo ""
echo "  # 3. Monitor the job"
echo "  squeue -u \$USER"
echo "  tail -f organoid_<JOBID>.out"
echo ""

#!/bin/bash
#
# ECM Stiffness Sweep using CompressSphere
# =========================================
# This script runs a parameter sweep over ECM stiffness using the CompressSphere
# application from OrganoidChaste. ECM stiffness is varied via:
#   1. Basal surface tension (GammaB) - cell-ECM interface stiffness
#   2. Soft potential boundary (HardPotential=0) - spring-like ECM substrate
#
# Two-phase approach:
#   Phase 1: Generate a relaxed spherical organoid mesh (MinimizeHomogeneousSphere)
#   Phase 2: Compress the sphere with varying ECM stiffness (CompressSphere)
#
# Usage:
#   chmod +x run_ecm_stiffness_sweep.sh
#   ./run_ecm_stiffness_sweep.sh
#

set -e

# ============================================================
# Configuration
# ============================================================

# Build directory (where Chaste was compiled)
BUILD_DIR="$HOME/Thesis"

# Path to executables
MINIMIZE_APP="${BUILD_DIR}/projects/OrganoidChaste/apps/MinimizeHomogeneousSphere"
COMPRESS_APP="${BUILD_DIR}/projects/OrganoidChaste/apps/CompressSphere"

# Number of cells in the organoid
NUM_CELLS=100

# Fixed surface tension parameters
GAMMA_A=0.9           # Apical tension (fixed)
GAMMA_LATERAL=1.0     # Lateral tension (fixed, normalized to 1)

# ECM stiffness sweep: varying basal tension (GammaB)
# Lower GammaB = softer ECM interface, Higher GammaB = stiffer ECM interface
GAMMA_B_VALUES=(0.3 0.5 0.7 0.9 1.1 1.3 1.5)

# Compression parameters
INDENTATION_PERCENT=0.7      # Compress to 70% of original diameter
BOUNDARY_VELOCITY=0.4        # Compression speed
USE_HARD_POTENTIAL=0          # 0 = soft potential (spring-like ECM)
DAMPING_CONSTANT_STEP=1.0    # Damping during compression
DAMPING_CONSTANT_RELAX=1.0   # Damping during relaxation

# Simulation timing
MAX_SIM_TIME=5000             # Max simulation time for convergence loop
TIME_INITIAL_RELAXATION=100   # Relaxation time after compression

# Output settings
PRINT_EVERY_N_COMPRESS=100    # Save frequency during compression
PRINT_EVERY_N_RELAX=200       # Save frequency during relaxation

# How many simulations to run in parallel (adjust to your CPU cores)
MAX_PARALLEL=4

# ============================================================
# Phase 1: Generate relaxed sphere mesh
# ============================================================

echo "=============================================="
echo " Phase 1: Generating relaxed organoid mesh"
echo "=============================================="
echo " NumCells: ${NUM_CELLS}"
echo " GammaA: ${GAMMA_A}, GammaB: 0.9 (default for mesh generation)"
echo ""

# Build if needed
cd "${BUILD_DIR}"
make MinimizeHomogeneousSphere 2>&1 | tail -3

# Run mesh generation with default parameters
# This creates a relaxed sphere that CompressSphere can read
${MINIMIZE_APP} \
    -NumCells ${NUM_CELLS} \
    -ApicalTension ${GAMMA_A} \
    -BasalTension 0.9 \
    -LateralTension ${GAMMA_LATERAL} \
    -MaxSimTime 2000 \
    -PrintEveryNSteps 500 \
    -AreaDiffCutoff 1e-4 \
    -FinalDT 0.1 \
    -Comment "_mesh_for_compression"

echo ""
echo "Phase 1 complete. Finding output mesh..."

# Find the most recently created mesh directory
MESH_DIR=$(ls -td ${HOME}/Thesis/testoutput/RandomizedHomogeneousSphere_*_mesh_for_compression_* 2>/dev/null | head -1)

if [ -z "${MESH_DIR}" ]; then
    echo "ERROR: Could not find mesh output directory"
    exit 1
fi

echo "Mesh directory: ${MESH_DIR}"

# ============================================================
# Phase 2: Compression sweep over ECM stiffness
# ============================================================

echo ""
echo "=============================================="
echo " Phase 2: ECM stiffness sweep (CompressSphere)"
echo "=============================================="
echo " Sweeping GammaB values: ${GAMMA_B_VALUES[*]}"
echo " Using soft potential (spring-like ECM boundary)"
echo " Indentation: ${INDENTATION_PERCENT}"
echo ""

# Build CompressSphere
make CompressSphere 2>&1 | tail -3

# Track running processes
PIDS=()
RUNNING=0

for GB in "${GAMMA_B_VALUES[@]}"; do
    echo "[$(date +%H:%M:%S)] Starting simulation: GammaB=${GB}"

    ${COMPRESS_APP} \
        -InputMeshPath "${MESH_DIR}/" \
        -GammaA ${GAMMA_A} \
        -GammaB ${GB} \
        -Indent ${INDENTATION_PERCENT} \
        -BoundaryVel ${BOUNDARY_VELOCITY} \
        -HardPotential ${USE_HARD_POTENTIAL} \
        -DCStep ${DAMPING_CONSTANT_STEP} \
        -DCRelax ${DAMPING_CONSTANT_RELAX} \
        -MaxSimTime ${MAX_SIM_TIME} \
        -TimeForInitialRelaxation ${TIME_INITIAL_RELAXATION} \
        -PrintEveryNStepsC ${PRINT_EVERY_N_COMPRESS} \
        -PrintEveryNStepsR ${PRINT_EVERY_N_RELAX} \
        -Comment "_ECM_sweep_GB${GB}" &

    PIDS+=($!)
    RUNNING=$((RUNNING + 1))

    # Limit parallel jobs
    if [ ${RUNNING} -ge ${MAX_PARALLEL} ]; then
        # Wait for any one process to finish
        wait -n 2>/dev/null || true
        RUNNING=$((RUNNING - 1))
    fi
done

# Wait for all remaining processes
echo ""
echo "Waiting for all simulations to complete..."
wait

echo ""
echo "=============================================="
echo " All simulations complete!"
echo "=============================================="
echo ""
echo "Output directories:"
ls -td ${HOME}/Thesis/testoutput/*Compress_Sphere*ECM_sweep* 2>/dev/null | head -20
echo ""
echo "Each directory contains .vtu files that can be"
echo "visualized in ParaView. Open the results_*.vtu"
echo "files as a time series to see the compression."
echo ""
echo "Parameters for each run are saved in"
echo "CompressSphereParameters.txt in each directory."

#!/bin/bash

# Run Painter Replication with OpenMP parallelization
# Uses all available CPU cores for cell force calculations

echo "=== Painter Replication with OpenMP Parallelization ==="
echo ""

# Detect number of cores
NUM_CORES=$(nproc)
echo "Available CPU cores: $NUM_CORES"

# Set OpenMP threads (use all cores)
export OMP_NUM_THREADS=$NUM_CORES

# Optional: Enable dynamic threading
export OMP_DYNAMIC=TRUE

# Optional: Set thread affinity for better performance
export OMP_PROC_BIND=true

echo "OpenMP threads: $OMP_NUM_THREADS"
echo ""

# Change to Chaste directory
cd /home/orlando/Thesis/Chaste

# Run the test
echo "Starting simulation..."
echo ""

/home/orlando/Thesis/projects/TissueMorphology/test/Test2dPainterReplication

echo ""
echo "=== Simulation Complete ==="
echo "Results: ~/Thesis/testoutput/PainterReplication2d/"

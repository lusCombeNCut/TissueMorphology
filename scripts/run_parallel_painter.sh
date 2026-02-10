#!/bin/bash

# Parallel Simulation Runner for Painter Replication
# Runs multiple simulations concurrently to utilize all CPU cores

echo "=== Parallel Painter Replication Simulation Runner ==="
echo "Available CPUs: $(nproc)"
echo ""

# Configuration
NUM_PARALLEL_JOBS=6   # Run 6 simulations concurrently (adjust based on your cores)
FAST_MODE=${FAST_MODE:-0}  # Set to 1 for reduced output

# Export for child processes
if [ "$FAST_MODE" -eq 1 ]; then
    export FAST_MODE=1
    echo "[FAST MODE ENABLED: Reduced I/O for faster execution]"
fi

# Test executable
TEST_EXE="/home/orlando/Thesis/projects/TissueMorphology/test/Test2dPainterReplication"

if [ ! -f "$TEST_EXE" ]; then
    echo "ERROR: Test executable not found at $TEST_EXE"
    echo "Please build first: cd /home/orlando/Thesis && make Test2dPainterReplication"
    exit 1
fi

echo ""
echo "=== Running simulations in parallel ===" 
echo "Jobs: $NUM_PARALLEL_JOBS concurrent processes"
echo ""

# Function to run a single simulation
run_simulation() {
    local run_num=$1
    local test_name=$2
    
    export RUN_NUMBER=$run_num
    echo "[Job $run_num] Starting $test_name..."
    
    cd /home/orlando/Thesis/Chaste
    "$TEST_EXE" --test $test_name > "/tmp/painter_run_${run_num}_${test_name}.log" 2>&1
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "[Job $run_num] ✓ $test_name completed successfully"
    else
        echo "[Job $run_num] ✗ $test_name failed (exit code: $exit_code)"
    fi
    
    return $exit_code
}

# Run RandomECM simulations in parallel
echo "--- Random ECM Simulations ---"
for i in $(seq 0 $((NUM_PARALLEL_JOBS-1))); do
    run_simulation $i "TestPainterRandomECM" &
    
    # Stagger starts to avoid I/O contention
    sleep 2
done

# Wait for all background jobs
wait

echo ""
echo "=== All simulations completed ==="
echo ""
echo "Results saved to: ~/Thesis/testoutput/PainterReplication2d/"
echo "Logs saved to: /tmp/painter_run_*.log"
echo ""
echo "To view logs: ls -lh /tmp/painter_run_*.log"
echo "To analyze results, run:"
echo "  python3 export_ecm_snapshots.py --data-dir ~/Thesis/testoutput/PainterReplication2d/random"

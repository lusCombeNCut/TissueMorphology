#!/bin/bash
# Run simulations in parallel using multiple CPU cores

# Configuration
NUM_CORES=10  # Use 10 cores (leave 2 for system)
NUM_RUNS=20   # 20 runs per ECM type
TEST_BINARY="./projects/TissueMorphology/test/Test2dDynamicECMInvasion"
LOG_DIR="$HOME/Thesis/simulation_logs"

# Clean up old outputs
echo "Cleaning old simulation outputs..."
rm -rf ~/Thesis/testoutput/DynamicECMInvasion2d
mkdir -p "$LOG_DIR"
rm -f "$LOG_DIR"/*.log

cd ~/Thesis

echo "========================================"
echo "Running 60 simulations in parallel"
echo "Using $NUM_CORES CPU cores"
echo "========================================"
echo ""

# Export function to run a single test
run_single_simulation() {
    ecm_type=$1
    run_num=$2
    
    log_file="$LOG_DIR/${ecm_type}_run_${run_num}.log"
    
    export RUN_NUMBER=$run_num
    
    case $ecm_type in
        perpendicular)
            $TEST_BINARY --test-name=TestDynamicPerpendicularECM > "$log_file" 2>&1
            ;;
        parallel)
            $TEST_BINARY --test-name=TestDynamicParallelECM > "$log_file" 2>&1
            ;;
        random)
            $TEST_BINARY --test-name=TestDynamicRandomECM > "$log_file" 2>&1
            ;;
    esac
    
    if [ $? -eq 0 ]; then
        echo "✓ Completed: $ecm_type run $run_num"
    else
        echo "✗ Failed: $ecm_type run $run_num (see $log_file)"
        return 1
    fi
}

export -f run_single_simulation
export TEST_BINARY
export LOG_DIR

# Create job list
job_list=$(mktemp)
for ecm_type in perpendicular parallel random; do
    for run_num in $(seq 0 19); do
        echo "$ecm_type $run_num" >> "$job_list"
    done
done

echo "Starting parallel execution at $(date)"
echo ""

# Check if GNU parallel is installed
if command -v parallel &> /dev/null; then
    # Use GNU parallel (more efficient)
    cat "$job_list" | parallel -j $NUM_CORES --colsep ' ' run_single_simulation {1} {2}
else
    # Fallback: use xargs (less efficient but still parallel)
    echo "GNU parallel not found, using xargs..."
    cat "$job_list" | xargs -n 2 -P $NUM_CORES bash -c 'run_single_simulation "$@"' _
fi

rm -f "$job_list"

echo ""
echo "========================================"
echo "All simulations completed at $(date)"
echo "========================================"
echo ""

# Summary
echo "Results summary:"
echo "----------------"
for ecm_type in perpendicular parallel random; do
    count=$(ls -d ~/Thesis/testoutput/DynamicECMInvasion2d/$ecm_type/run_* 2>/dev/null | wc -l)
    echo "$ecm_type: $count / 20 runs completed"
done

echo ""
echo "Logs saved to: $LOG_DIR"
echo ""
echo "To generate plots, run:"
echo "  cd ~/Thesis/Chaste/projects/TissueMorphology"
echo "  python3 plot_invasion_analysis_averaged.py"

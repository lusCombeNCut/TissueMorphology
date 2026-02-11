#!/bin/bash
# Check progress of multiple simulation runs

echo "=== Simulation Progress Monitor ==="
echo ""

# Check if process is running
if ps aux | grep -q "[T]est2dDynamicECMInvasion"; then
    echo "✓ Simulation process is RUNNING"
    echo ""
else
    echo "✗ Simulation process is NOT running"
    echo ""
fi

# Count completed run directories
base_dir=~/Thesis/testoutput/DynamicECMInvasion2d

if [ -d "$base_dir" ]; then
    perpendicular_count=$(ls -d $base_dir/perpendicular/run_* 2>/dev/null | wc -l)
    parallel_count=$(ls -d $base_dir/parallel/run_* 2>/dev/null | wc -l)
    random_count=$(ls -d $base_dir/random/run_* 2>/dev/null | wc -l)
    total_count=$((perpendicular_count + parallel_count + random_count))
    
    echo "Completed runs:"
    echo "  Perpendicular: $perpendicular_count / 20"
    echo "  Parallel:      $parallel_count / 20"
    echo "  Random:        $random_count / 20"
    echo "  Total:         $total_count / 60"
    echo ""
    
    # Calculate percentage
    percentage=$((total_count * 100 / 60))
    echo "Progress: $percentage% complete"
    echo ""
else
    echo "Output directory not found yet"
    echo ""
fi

# Show last few lines of log
echo "Last 10 lines of simulation log:"
echo "=================================="
tail -10 ~/Thesis/simulation_log.txt 2>/dev/null || echo "Log file not found"

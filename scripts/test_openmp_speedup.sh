#!/bin/bash

# Quick test to verify OpenMP is working

echo "=== OpenMP Parallelization Test ==="
echo ""

# Test with 1 thread
echo "Running with 1 thread..."
export OMP_NUM_THREADS=1
time -p /home/orlando/Thesis/projects/TissueMorphology/test/Test2dPainterReplication --test TestPainterRandomECM > /dev/null 2>&1
SERIAL_TIME=$?

echo ""

# Test with all threads
CORES=$(nproc)
echo "Running with $CORES threads..."
export OMP_NUM_THREADS=$CORES
time -p /home/orlando/Thesis/projects/TissueMorphology/test/Test2dPainterReplication --test TestPainterRandomECM > /dev/null 2>&1
PARALLEL_TIME=$?

echo ""
echo "=== Speedup Analysis ==="
echo "Compare the 'real' times above to see the speedup from parallelization"
echo "Expected speedup: 2-6x with $CORES cores"

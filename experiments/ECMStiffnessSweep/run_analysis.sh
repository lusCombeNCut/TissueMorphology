#!/bin/bash
# =============================================================================
# run_analysis.sh — Standalone analysis job for ECM stiffness sweep
#
# Usage:
#   sbatch run_analysis.sh <analysis_dir>
#
# Example:
#   BASE="/user/work/sv22482/sim_output/analysis_ecm_stiffness_2026-03-30_01-36-37"
#   sbatch run_analysis.sh "${BASE}"
#
# <analysis_dir> must contain a merged/ subdirectory produced by submit_sweep.sh.
# =============================================================================

#SBATCH --account=semt036404
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=Analyse_ECMStiffness

# ---------- Argument ----------
if [ -z "$1" ]; then
    echo "ERROR: No analysis directory provided."
    echo "Usage: sbatch run_analysis.sh <analysis_dir>"
    exit 1
fi

BASE="$1"
MERGED_DIR="${BASE}/merged"
ANALYSIS_OUT="${BASE}/crypt_analysis_output"
TIMESTEP_OUT="${BASE}/timestep_analysis_output"
PLOTS_OUT="${BASE}/plots"
SCRIPT_DIR="/user/work/$(whoami)/TissueMorphology/experiments/ECMStiffnessSweep"
LOG_DIR="${BASE}"

# Set log file paths (can't use $1 in #SBATCH directives, so rename after)
# Slurm will write to default location; we tee everything ourselves below.

# ---------- Setup ----------
module load languages/python 2>/dev/null || module load python 2>/dev/null || true

mkdir -p "${ANALYSIS_OUT}" "${TIMESTEP_OUT}" "${PLOTS_OUT}"

echo "============================================"
echo "  Job ID:      ${SLURM_JOB_ID}"
echo "  Base dir:    ${BASE}"
echo "  Merged dir:  ${MERGED_DIR}"
echo "  Start:       $(date)"
echo "============================================"

if [ ! -d "${MERGED_DIR}" ]; then
    echo "ERROR: merged/ not found at ${MERGED_DIR}"
    exit 1
fi

# ---------- Crypt budding analysis ----------
echo ""
echo "=== Crypt budding analysis ==="
python3 "${SCRIPT_DIR}/analyse_crypt_budding.py" \
    --base-dir "${MERGED_DIR}" \
    -o "${ANALYSIS_OUT}" \
    --method fourier \
    --use-outline \
    --debug-plots \
    2>&1 | tee "${ANALYSIS_OUT}/crypt_analysis.log" || true

# ---------- Timestep summary ----------
echo ""
echo "=== Timestep summary ==="
python3 "${SCRIPT_DIR}/timestep_summary.py" \
    --base-dir "${MERGED_DIR}" \
    -o "${TIMESTEP_OUT}" \
    2>&1 | tee "${TIMESTEP_OUT}/timestep_analysis.log" || true

# ---------- Summary plots ----------
echo ""
echo "=== Summary plots ==="
python3 "${SCRIPT_DIR}/analyse_and_plot.py" \
    --data-dir "${MERGED_DIR}" \
    -o "${PLOTS_OUT}" \
    2>&1 | tee "${PLOTS_OUT}/plot_generation.log" || true

echo ""
echo "============================================"
echo "  Done. End: $(date)"
echo "  Plots:    ${PLOTS_OUT}/"
echo "  To download plots:"
echo "    scp -r sv22482@bp1-login.acrc.bris.ac.uk:${PLOTS_OUT}/ ./plots/"
echo "============================================"

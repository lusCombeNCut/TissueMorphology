#!/bin/bash
#
# fetch_and_analyse.sh
#
# Download simulation results from BluePebble and run analysis scripts.
#
# Usage:
#   ./fetch_and_analyse.sh --job-dir <job_folder> [--job-dir <job_folder2> ...] [options]
#
# Required (unless --skip-scp):
#   --job-dir <folder>     Job subdirectory on BluePebble (can be repeated)
#                          (e.g., 16514549_node2d_2026-03-28_16-59-42)
#                          Full path: /user/work/sv22482/sim_output/<job_folder>/archives/
#
# Options:
#   --model <type>         Only fetch/analyse specific model: node2d, vertex2d, node3d, vertex3d
#   --runs <list>          Only fetch specific runs (comma-separated, e.g., "0,1,2" or "0")
#   --remote               Run analysis on BluePebble, then download plots only
#   --skip-scp             Skip SCP download step (analyse existing data only)
#   --skip-analysis        Skip analysis step (SCP only)
#   --fix-pvd              Fix broken PVD files during analysis
#   --output-base <dir>    Local output directory (default: ../../output)
#   --help                 Show this help message
#
# Examples:
#   # Run analysis on BluePebble and just download the plots:
#   ./fetch_and_analyse.sh --remote \
#     --job-dir 16514549_node2d_2026-03-28_16-59-42 \
#     --job-dir 16514550_vertex2d_2026-03-28_16-59-42 \
#     --job-dir 16514551_node3d_2026-03-28_16-59-42 \
#     --job-dir 16514552_vertex3d_2026-03-28_16-59-42
#
#   # Download all raw data and analyse locally:
#   ./fetch_and_analyse.sh \
#     --job-dir 16514549_node2d_2026-03-28_16-59-42 \
#     --job-dir 16514550_vertex2d_2026-03-28_16-59-42 \
#     --job-dir 16514551_node3d_2026-03-28_16-59-42 \
#     --job-dir 16514552_vertex3d_2026-03-28_16-59-42
#
#   # Analyse existing local data without downloading:
#   ./fetch_and_analyse.sh --skip-scp --model vertex3d --output-base ../../output
#

set -e  # Exit on error

# =============================================================================
# Default values
# =============================================================================
REMOTE_HOST="sv22482@bp1-login.acrc.bris.ac.uk"
REMOTE_BASE="/user/work/sv22482/sim_output"
REMOTE_WORK="/user/work/sv22482"
JOB_DIRS=()
OUTPUT_BASE="../../output"
MODEL_FILTER=""
RUN_FILTER=""
REMOTE_MODE=false
SKIP_SCP=false
SKIP_ANALYSIS=false
FIX_PVD=false

# =============================================================================
# Parse command line arguments
# =============================================================================
while [[ $# -gt 0 ]]; do
  case $1 in
    --job-dir)
      JOB_DIRS+=("$2")
      shift 2
      ;;
    --model)
      MODEL_FILTER="$2"
      if [[ ! "$MODEL_FILTER" =~ ^(node2d|vertex2d|node3d|vertex3d)$ ]]; then
        echo "ERROR: Invalid model '$MODEL_FILTER'. Must be one of: node2d, vertex2d, node3d, vertex3d"
        exit 1
      fi
      shift 2
      ;;
    --runs)
      RUN_FILTER="$2"
      shift 2
      ;;
    --remote)
      REMOTE_MODE=true
      shift
      ;;
    --skip-scp)
      SKIP_SCP=true
      shift
      ;;
    --skip-analysis)
      SKIP_ANALYSIS=true
      shift
      ;;
    --fix-pvd)
      FIX_PVD=true
      shift
      ;;
    --output-base)
      OUTPUT_BASE="$2"
      shift 2
      ;;
    --help|-h)
      grep '^#' "$0" | grep -v '#!/bin/bash' | sed 's/^# //' | sed 's/^#//'
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# =============================================================================
# Validate arguments
# =============================================================================
if [[ ${#JOB_DIRS[@]} -eq 0 && "$SKIP_SCP" == false ]]; then
  echo "ERROR: At least one --job-dir is required (unless --skip-scp is used)"
  echo "Use --help for usage information"
  exit 1
fi

# =============================================================================
# Setup paths
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPERIMENTS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ "$REMOTE_MODE" == true ]]; then
  # For remote mode, output goes to a local directory for downloaded plots
  mkdir -p "$SCRIPT_DIR/$OUTPUT_BASE"
  OUTPUT_DIR="$(cd "$SCRIPT_DIR" && cd "$OUTPUT_BASE" && pwd)"
else
  OUTPUT_DIR="$(cd "$SCRIPT_DIR" && cd "$OUTPUT_BASE" && pwd)"
fi

MERGED_DIR="$OUTPUT_DIR/merged"
ANALYSIS_OUT="$OUTPUT_DIR/crypt_analysis_output"
TIMESTEP_OUT="$OUTPUT_DIR/timestep_analysis_output"

echo "======================================================================="
echo "  BluePebble Data Fetch & Analysis"
echo "======================================================================="
echo "  Script directory:  $SCRIPT_DIR"
if [[ "$REMOTE_MODE" == true ]]; then
  echo "  Mode:              REMOTE (analyse on BluePebble, download plots)"
else
  echo "  Mode:              LOCAL (download data, analyse here)"
fi
echo "  Output directory:  $OUTPUT_DIR"
if [[ "$REMOTE_MODE" == true || "$SKIP_SCP" == false ]]; then
  echo "  Remote host:       $REMOTE_HOST"
  echo "  Job directories:   ${JOB_DIRS[*]}"
fi
echo "  Model filter:      ${MODEL_FILTER:-all}"
echo "  Run filter:        ${RUN_FILTER:-all}"
echo "  Fix PVD files:     $FIX_PVD"
echo "======================================================================="
echo ""

# #############################################################################
# REMOTE MODE: run analysis on BluePebble, download plots
# #############################################################################
if [[ "$REMOTE_MODE" == true ]]; then

  REMOTE_REPO="${REMOTE_WORK}/TissueMorphology"
  REMOTE_SCRIPTS="${REMOTE_REPO}/experiments/ECMStiffnessSweep"
  REMOTE_OUTPUT="${REMOTE_BASE}/analysis_output"

  # Step 1: SSH in and run merge + analysis
  echo "======================================================================="
  echo "  Step 1: Running analysis on BluePebble"
  echo "======================================================================="
  echo "  (Make sure you have pushed locally and pulled on BluePebble first)"
  echo ""

  # Build job dirs as space-separated string for the remote script
  JOB_DIRS_STR="${JOB_DIRS[*]}"

  ssh -t "$REMOTE_HOST" bash -s <<REMOTE_SCRIPT
set -e

REMOTE_BASE="${REMOTE_BASE}"
REMOTE_SCRIPTS="${REMOTE_REPO}/experiments/ECMStiffnessSweep"
REMOTE_OUTPUT="${REMOTE_OUTPUT}"
MODEL_FILTER="${MODEL_FILTER}"
RUN_FILTER="${RUN_FILTER}"
FIX_PVD="${FIX_PVD}"
JOB_DIRS_STR="${JOB_DIRS_STR}"

read -ra JOB_DIRS <<< "\$JOB_DIRS_STR"

echo "  Working on BluePebble..."
echo "  Job directories: \${JOB_DIRS[*]}"

# --- Create merged directory ---
MERGED_DIR="\${REMOTE_OUTPUT}/merged"
rm -rf "\$MERGED_DIR"

MODELS_TO_MERGE=()
if [[ -n "\$MODEL_FILTER" ]]; then
  MODELS_TO_MERGE=("\$MODEL_FILTER")
else
  for model in node2d vertex2d node3d vertex3d; do
    for JOB_DIR in "\${JOB_DIRS[@]}"; do
      if ls -d "\${REMOTE_BASE}/\${JOB_DIR}"/s*_r*/CryptBudding/*/"\$model" 2>/dev/null | head -1 | grep -q .; then
        MODELS_TO_MERGE+=("\$model")
        break
      fi
    done
  done
fi

echo "  Models found: \${MODELS_TO_MERGE[*]}"

for model in "\${MODELS_TO_MERGE[@]}"; do
  echo "  Merging \$model..."
  MERGED_MODEL="\$MERGED_DIR/CryptBudding/\$model"

  for JOB_DIR in "\${JOB_DIRS[@]}"; do
    SEARCH_DIR="\${REMOTE_BASE}/\${JOB_DIR}"
    for sdir in "\$SEARCH_DIR"/s*_r*/CryptBudding/*/"\$model"/stiffness_*; do
      [[ -d "\$sdir" ]] || continue
      stiff=\$(basename "\$sdir")
      mkdir -p "\$MERGED_MODEL/\$stiff"

      for rdir in "\$sdir"/run_*; do
        [[ -d "\$rdir" ]] || continue
        run=\$(basename "\$rdir")
        run_num="\${run#run_}"

        if [[ -n "\$RUN_FILTER" ]]; then
          if [[ ! ",\$RUN_FILTER," =~ ,"\$run_num", ]]; then
            continue
          fi
        fi

        if [[ ! -e "\$MERGED_MODEL/\$stiff/\$run" ]]; then
          ln -s "\$rdir" "\$MERGED_MODEL/\$stiff/\$run"
        fi
      done
    done
  done
done

echo "  Merged directory: \$MERGED_DIR"

# --- Load Python environment ---
module load languages/python 2>/dev/null || module load python 2>/dev/null || true

# --- Run analysis ---
ANALYSIS_OUT="\${REMOTE_OUTPUT}/crypt_analysis_output"
TIMESTEP_OUT="\${REMOTE_OUTPUT}/timestep_analysis_output"
PLOTS_OUT="\${REMOTE_OUTPUT}/plots"
mkdir -p "\$ANALYSIS_OUT" "\$TIMESTEP_OUT" "\$PLOTS_OUT"

MODEL_ARG=""
[[ -n "\$MODEL_FILTER" ]] && MODEL_ARG="--model \$MODEL_FILTER"

echo ""
echo "  Running crypt count analysis..."
python3 "\$REMOTE_SCRIPTS/analyse_crypt_budding.py" \
  --base-dir "\$MERGED_DIR" \
  \$MODEL_ARG \
  -o "\$ANALYSIS_OUT" \
  --method fourier \
  --use-outline \
  --debug-plots \
  2>&1 | tee "\$ANALYSIS_OUT/analysis.log"

echo ""
echo "  Running timestep summary analysis..."
TIMESTEP_ARGS="--base-dir \$MERGED_DIR \$MODEL_ARG -o \$TIMESTEP_OUT"
[[ "\$FIX_PVD" == true ]] && TIMESTEP_ARGS="\$TIMESTEP_ARGS --fix-pvd"
python3 "\$REMOTE_SCRIPTS/timestep_summary.py" \
  \$TIMESTEP_ARGS \
  2>&1 | tee "\$TIMESTEP_OUT/timestep_analysis.log"

echo ""
echo "  Running summary plot generation..."
python3 "\$REMOTE_SCRIPTS/analyse_and_plot.py" \
  --data-dir "\$MERGED_DIR" \
  \$MODEL_ARG \
  -o "\$PLOTS_OUT" \
  2>&1 | tee "\$PLOTS_OUT/plot_generation.log"

echo ""
echo "  Remote analysis complete."
echo "  Plots at: \$PLOTS_OUT"
echo "  Crypt counts at: \$ANALYSIS_OUT"
REMOTE_SCRIPT

  echo ""

  # Step 2: Download plots
  echo "======================================================================="
  echo "  Step 2: Downloading results"
  echo "======================================================================="

  mkdir -p "$OUTPUT_DIR"

  rsync -avz \
    "${REMOTE_HOST}:${REMOTE_OUTPUT}/plots/" \
    "$OUTPUT_DIR/plots/"

  rsync -avz \
    "${REMOTE_HOST}:${REMOTE_OUTPUT}/crypt_analysis_output/" \
    "$OUTPUT_DIR/crypt_analysis_output/"

  rsync -avz \
    "${REMOTE_HOST}:${REMOTE_OUTPUT}/timestep_analysis_output/" \
    "$OUTPUT_DIR/timestep_analysis_output/"

  echo ""
  echo "======================================================================="
  echo "  Done! Results downloaded to:"
  echo "    Plots:           $OUTPUT_DIR/plots/"
  echo "    Crypt counts:    $OUTPUT_DIR/crypt_analysis_output/"
  echo "    Timestep summary:$OUTPUT_DIR/timestep_analysis_output/"
  echo "======================================================================="

  exit 0
fi

# #############################################################################
# LOCAL MODE: download raw data, analyse locally
# #############################################################################

# =============================================================================
# SCP: Download archives from BluePebble
# =============================================================================
if [[ "$SKIP_SCP" == false ]]; then
  echo "======================================================================="
  echo "  Step 1: Downloading archives from BluePebble"
  echo "======================================================================="

  # Build rsync source paths — one per job dir, using --relative so rsync
  # preserves the job directory structure in a single authenticated transfer.
  # The /. marker tells rsync where the relative path starts.
  RSYNC_SOURCES=()
  for JOB_DIR in "${JOB_DIRS[@]}"; do
    RSYNC_SOURCES+=("${REMOTE_HOST}:${REMOTE_BASE}/./${JOB_DIR}/archives/")
  done

  RSYNC_FILTER=()
  if [[ -n "$RUN_FILTER" ]]; then
    RSYNC_FILTER+=(--include='*/')
    IFS=',' read -ra RUN_ARRAY <<< "$RUN_FILTER"
    for run in "${RUN_ARRAY[@]}"; do
      run=$(echo "$run" | xargs)
      RSYNC_FILTER+=(--include="s*_r${run}.tar.gz")
    done
    RSYNC_FILTER+=(--exclude='*')
  fi

  echo "  Downloading ${#JOB_DIRS[@]} job directories in one transfer..."
  rsync -avz --relative "${RSYNC_FILTER[@]}" "${RSYNC_SOURCES[@]}" "$OUTPUT_DIR/" || {
    echo "ERROR: rsync transfer failed"; exit 1
  }

  # Move archives from archives/ subdirs up to job dir level
  for JOB_DIR in "${JOB_DIRS[@]}"; do
    if [[ -d "$OUTPUT_DIR/$JOB_DIR/archives" ]]; then
      mv "$OUTPUT_DIR/$JOB_DIR/archives/"*.tar.gz "$OUTPUT_DIR/$JOB_DIR/" 2>/dev/null || true
      rm -rf "$OUTPUT_DIR/$JOB_DIR/archives"
    fi
    NUM_ARCHIVES=$(ls "$OUTPUT_DIR/$JOB_DIR/"*.tar.gz 2>/dev/null | wc -l)
    echo "  $JOB_DIR: $NUM_ARCHIVES archives"
  done

  echo ""

  # Extract archives from all job directories
  echo "======================================================================="
  echo "  Step 2: Extracting archives"
  echo "======================================================================="

  for JOB_DIR in "${JOB_DIRS[@]}"; do
    JOB_OUTPUT_DIR="$OUTPUT_DIR/$JOB_DIR"
    cd "$JOB_OUTPUT_DIR"
    for archive in *.tar.gz; do
      [[ -e "$archive" ]] || continue
      echo "  Extracting: $JOB_DIR/$archive"
      tar xzf "$archive"
    done
  done
  echo "  Extraction complete"
  echo ""
fi

# =============================================================================
# Create merged directory structure
# =============================================================================
if [[ "$SKIP_ANALYSIS" == false ]]; then
  echo "======================================================================="
  echo "  Step 3: Creating merged directory structure"
  echo "======================================================================="

  # Remove old merged directory if it exists
  rm -rf "$MERGED_DIR"

  # Determine search directories
  SEARCH_DIRS=()
  if [[ ${#JOB_DIRS[@]} -gt 0 ]]; then
    for JOB_DIR in "${JOB_DIRS[@]}"; do
      dir="$OUTPUT_DIR/$JOB_DIR"
      [[ -d "$dir" ]] && SEARCH_DIRS+=("$dir")
    done
  fi
  # If no job dirs specified (--skip-scp), find all job directories
  if [[ ${#SEARCH_DIRS[@]} -eq 0 ]]; then
    SEARCH_DIRS=("$OUTPUT_DIR"/*_*_*-*-*_*-*-*/)
  fi

  # Determine which models to merge
  MODELS_TO_MERGE=()
  if [[ -n "$MODEL_FILTER" ]]; then
    MODELS_TO_MERGE=("$MODEL_FILTER")
  else
    for model in node2d vertex2d node3d vertex3d; do
      for search_dir in "${SEARCH_DIRS[@]}"; do
        if ls -d "$search_dir"/s*_r*/CryptBudding/*/"$model" 2>/dev/null | head -1 | grep -q .; then
          MODELS_TO_MERGE+=("$model")
          break
        fi
      done
    done
  fi

  echo "  Models to merge: ${MODELS_TO_MERGE[*]}"

  for model in "${MODELS_TO_MERGE[@]}"; do
    echo "  Merging $model..."
    MERGED_MODEL="$MERGED_DIR/CryptBudding/$model"

    for search_dir in "${SEARCH_DIRS[@]}"; do
    for sdir in "$search_dir"/s*_r*/CryptBudding/*/"$model"/stiffness_*; do
      [[ -d "$sdir" ]] || continue

      stiff=$(basename "$sdir")
      mkdir -p "$MERGED_MODEL/$stiff"

      for rdir in "$sdir"/run_*; do
        [[ -d "$rdir" ]] || continue

        run=$(basename "$rdir")
        run_num="${run#run_}"

        # Apply run filter if specified
        if [[ -n "$RUN_FILTER" ]]; then
          if [[ ! ",$RUN_FILTER," =~ ,"$run_num", ]]; then
            continue
          fi
        fi

        if [[ ! -e "$MERGED_MODEL/$stiff/$run" ]]; then
          ln -s "$rdir" "$MERGED_MODEL/$stiff/$run"
        fi
      done
    done
    done
  done

  echo "  Merged directory created at: $MERGED_DIR"
  echo ""

  # =============================================================================
  # Run analysis scripts
  # =============================================================================
  echo "======================================================================="
  echo "  Step 4: Running analysis scripts"
  echo "======================================================================="

  mkdir -p "$ANALYSIS_OUT"
  mkdir -p "$TIMESTEP_OUT"

  # Build model argument for analysis scripts
  MODEL_ARG=""
  if [[ -n "$MODEL_FILTER" ]]; then
    MODEL_ARG="--model $MODEL_FILTER"
  fi

  # Run crypt budding analysis
  echo "  Running crypt count analysis..."
  python3 "$SCRIPT_DIR/analyse_crypt_budding.py" \
    --base-dir "$MERGED_DIR" \
    $MODEL_ARG \
    -o "$ANALYSIS_OUT" \
    --method fourier \
    --use-outline \
    --debug-plots \
    2>&1 | tee "$ANALYSIS_OUT/analysis.log"

  echo ""

  # Run timestep summary analysis
  echo "  Running timestep summary analysis..."
  TIMESTEP_ARGS="--base-dir $MERGED_DIR $MODEL_ARG -o $TIMESTEP_OUT"
  if [[ "$FIX_PVD" == true ]]; then
    TIMESTEP_ARGS="$TIMESTEP_ARGS --fix-pvd"
  fi

  python3 "$SCRIPT_DIR/timestep_summary.py" \
    $TIMESTEP_ARGS \
    2>&1 | tee "$TIMESTEP_OUT/timestep_analysis.log"

  echo ""

  # Run comprehensive summary plots (time-series, boxplots, comparisons)
  PLOTS_OUT="$OUTPUT_DIR/plots"
  mkdir -p "$PLOTS_OUT"

  echo "  Running summary plot generation..."
  python3 "$SCRIPT_DIR/analyse_and_plot.py" \
    --data-dir "$MERGED_DIR" \
    $MODEL_ARG \
    -o "$PLOTS_OUT" \
    2>&1 | tee "$PLOTS_OUT/plot_generation.log"

  echo ""
  echo "======================================================================="
  echo "  Analysis complete!"
  echo "======================================================================="
  echo "  Crypt counts:     $ANALYSIS_OUT/"
  echo "  Timestep summary: $TIMESTEP_OUT/"
  echo "  Summary plots:    $PLOTS_OUT/"
  echo "======================================================================="
fi

echo ""
echo "Done."

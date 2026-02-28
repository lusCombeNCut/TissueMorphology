#!/bin/bash
#
# fetch_and_analyse.sh
#
# Download simulation results from BluePebble and run analysis scripts.
#
# Usage:
#   ./fetch_and_analyse.sh --job-dir <job_folder> [options]
#
# Required:
#   --job-dir <folder>     Job subdirectory on BluePebble
#                          (e.g., 15712463_vertex3d_2026-02-17_00-23-53)
#                          Full path: /user/work/sv22482/sim_output/<job_folder>/archives/
#
# Options:
#   --model <type>         Only fetch/analyse specific model: node2d, vertex2d, node3d, vertex3d
#   --runs <list>          Only fetch specific runs (comma-separated, e.g., "0,1,2" or "0")
#   --skip-scp             Skip SCP download step (analyse existing data only)
#   --skip-analysis        Skip analysis step (SCP only)
#   --fix-pvd              Fix broken PVD files during analysis
#   --output-base <dir>    Local output directory (default: ../../output)
#   --help                 Show this help message
#
# Examples:
#   # Download all archives and analyse everything:
#   ./fetch_and_analyse.sh --job-dir 15712463_vertex3d_2026-02-17_00-23-53
#
#   # Download only run 0 from each parameter, analyse vertex3d only:
#   ./fetch_and_analyse.sh --job-dir 15712463_vertex3d_2026-02-17_00-23-53 \
#                          --model vertex3d --runs 0
#
#   # Download from vertex2d job:
#   ./fetch_and_analyse.sh --job-dir 15712513_vertex2d_2026-02-17_00-50-10 --model vertex2d
#
#   # Analyse existing data without downloading:
#   ./fetch_and_analyse.sh --skip-scp --model vertex3d --output-base ../../output
#

set -e  # Exit on error

# =============================================================================
# Default values
# =============================================================================
REMOTE_HOST="sv22482@bp1-login.acrc.bris.ac.uk"
REMOTE_BASE="/user/work/sv22482/sim_output"
JOB_DIR=""
OUTPUT_BASE="../../output"
MODEL_FILTER=""
RUN_FILTER=""
SKIP_SCP=false
SKIP_ANALYSIS=false
FIX_PVD=false

# =============================================================================
# Parse command line arguments
# =============================================================================
while [[ $# -gt 0 ]]; do
  case $1 in
    --job-dir)
      JOB_DIR="$2"
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
if [[ "$SKIP_SCP" == false && -z "$JOB_DIR" ]]; then
  echo "ERROR: --job-dir is required (unless --skip-scp is used)"
  echo "Use --help for usage information"
  exit 1
fi

# =============================================================================
# Setup paths
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$(cd "$SCRIPT_DIR" && cd "$OUTPUT_BASE" && pwd)"
MERGED_DIR="$OUTPUT_DIR/merged"
ANALYSIS_OUT="$OUTPUT_DIR/crypt_analysis_output"
TIMESTEP_OUT="$OUTPUT_DIR/timestep_analysis_output"

# Construct full remote path
REMOTE_DIR=""
if [[ -n "$JOB_DIR" ]]; then
  REMOTE_DIR="$REMOTE_BASE/$JOB_DIR"
fi

echo "======================================================================="
echo "  BluePebble Data Fetch & Analysis"
echo "======================================================================="
echo "  Script directory:  $SCRIPT_DIR"
echo "  Output directory:  $OUTPUT_DIR"
if [[ "$SKIP_SCP" == false ]]; then
  echo "  Remote host:       $REMOTE_HOST"
  echo "  Job directory:     $JOB_DIR"
  echo "  Full remote path:  $REMOTE_DIR/archives/"
fi
echo "  Model filter:      ${MODEL_FILTER:-all}"
echo "  Run filter:        ${RUN_FILTER:-all}"
echo "  Skip SCP:          $SKIP_SCP"
echo "  Skip analysis:     $SKIP_ANALYSIS"
echo "  Fix PVD files:     $FIX_PVD"
echo "======================================================================="
echo ""

# =============================================================================
# SCP: Download archives from BluePebble
# =============================================================================
if [[ "$SKIP_SCP" == false ]]; then
  echo "======================================================================="
  echo "  Step 1: Downloading archives from BluePebble"
  echo "======================================================================="
  
  ARCHIVES_PATH="${REMOTE_DIR}/archives"
  
  # Build archive pattern for filtering
  if [[ -n "$RUN_FILTER" ]]; then
    # User specified specific runs
    IFS=',' read -ra RUN_ARRAY <<< "$RUN_FILTER"
    PATTERNS=()
    for run in "${RUN_ARRAY[@]}"; do
      run=$(echo "$run" | xargs)  # Trim whitespace
      PATTERNS+=("s*_r${run}.tar.gz")
    done
    
    echo "  Fetching archives for runs: ${RUN_FILTER}"
    for pattern in "${PATTERNS[@]}"; do
      echo "    Pattern: $pattern"
      scp "${REMOTE_HOST}:${ARCHIVES_PATH}/${pattern}" "$OUTPUT_DIR/" || true
    done
  else
    # Fetch all archives
    echo "  Fetching all archives from: ${REMOTE_HOST}:${ARCHIVES_PATH}/"
    scp "${REMOTE_HOST}:${ARCHIVES_PATH}/*.tar.gz" "$OUTPUT_DIR/" || {
      echo "  ERROR: Failed to download archives"
      exit 1
    }
  fi
  
  # Count downloaded archives
  NUM_ARCHIVES=$(ls "$OUTPUT_DIR"/*.tar.gz 2>/dev/null | wc -l)
  echo "  Downloaded $NUM_ARCHIVES archives"
  echo ""
  
  # Extract archives
  echo "======================================================================="
  echo "  Step 2: Extracting archives"
  echo "======================================================================="
  
  cd "$OUTPUT_DIR"
  for archive in *.tar.gz; do
    [[ -e "$archive" ]] || continue
    echo "  Extracting: $archive"
    tar xzf "$archive"
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
  
  # Determine which models to merge
  MODELS_TO_MERGE=()
  if [[ -n "$MODEL_FILTER" ]]; then
    MODELS_TO_MERGE=("$MODEL_FILTER")
  else
    # Auto-detect models from extracted directories
    # App output: s*_r*/CryptBudding/<git_hash>/<model>/stiffness_*/run_*
    for model in node2d vertex2d node3d vertex3d; do
      if ls -d "$OUTPUT_DIR"/s*_r*/CryptBudding/*/"$model" 2>/dev/null | head -1 | grep -q .; then
        MODELS_TO_MERGE+=("$model")
      fi
    done
  fi
  
  echo "  Models to merge: ${MODELS_TO_MERGE[*]}"
  
  for model in "${MODELS_TO_MERGE[@]}"; do
    echo "  Merging $model..."
    MERGED_MODEL="$MERGED_DIR/CryptBudding/$model"
    
    # Create symlinks at the run level
    # App output: s*_r*/CryptBudding/<git_hash>/<model>/stiffness_*/run_*
    for sdir in "$OUTPUT_DIR"/s*_r*/CryptBudding/*/"$model"/stiffness_*; do
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
            continue  # Skip this run
          fi
        fi
        
        if [[ ! -e "$MERGED_MODEL/$stiff/$run" ]]; then
          ln -s "$rdir" "$MERGED_MODEL/$stiff/$run"
        fi
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
  echo "======================================================================="
  echo "  Analysis complete!"
  echo "======================================================================="
  echo "  Crypt counts:     $ANALYSIS_OUT/"
  echo "  Timestep summary: $TIMESTEP_OUT/"
  echo "======================================================================="
fi

echo ""
echo "Done."

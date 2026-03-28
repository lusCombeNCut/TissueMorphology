#!/bin/bash
#
# fetch_and_analyse.sh — Lumen Pressure Sweep
#
# Download results from BluePebble and run analysis + plotting.
#
# Usage:
#   ./fetch_and_analyse.sh --job-dir <folder> [--model <type>] [--skip-scp] [--full]
#   ./fetch_and_analyse.sh --skip-scp --output-base ../../output --model node2d

set -e

REMOTE_HOST="sv22482@bp1-login.acrc.bris.ac.uk"
REMOTE_BASE="/user/work/sv22482/sim_output"
JOB_DIR=""
OUTPUT_BASE="../../output"
MODEL_FILTER=""
RUN_FILTER=""
SKIP_SCP=false
SKIP_ANALYSIS=false
FULL_ANALYSIS=false

while [[ $# -gt 0 ]]; do
  case $1 in
    --job-dir)      JOB_DIR="$2"; shift 2 ;;
    --model)        MODEL_FILTER="$2"; shift 2 ;;
    --runs)         RUN_FILTER="$2"; shift 2 ;;
    --skip-scp)     SKIP_SCP=true; shift ;;
    --skip-analysis) SKIP_ANALYSIS=true; shift ;;
    --full)         FULL_ANALYSIS=true; shift ;;
    --output-base)  OUTPUT_BASE="$2"; shift 2 ;;
    --help|-h)      grep '^#' "$0" | sed 's/^# //' | sed 's/^#//'; exit 0 ;;
    *)              echo "ERROR: Unknown option: $1"; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$(cd "$SCRIPT_DIR" && cd "$OUTPUT_BASE" && pwd 2>/dev/null || echo "$OUTPUT_BASE")"

if [[ "$SKIP_SCP" == false && -z "$JOB_DIR" ]]; then
  echo "ERROR: --job-dir required unless --skip-scp"
  exit 1
fi

echo "======================================================================="
echo "  Lumen Pressure Sweep — Fetch & Analyse"
echo "======================================================================="

# --- SCP Download ---
if [[ "$SKIP_SCP" == false ]]; then
  JOB_OUTPUT_DIR="$OUTPUT_DIR/$JOB_DIR"
  mkdir -p "$JOB_OUTPUT_DIR"

  ARCHIVES_PATH="${REMOTE_BASE}/${JOB_DIR}/archives"

  if [[ -n "$RUN_FILTER" ]]; then
    IFS=',' read -ra RUN_ARRAY <<< "$RUN_FILTER"
    for run in "${RUN_ARRAY[@]}"; do
      run=$(echo "$run" | xargs)
      scp "${REMOTE_HOST}:${ARCHIVES_PATH}/p*_r${run}.tar.gz" "$JOB_OUTPUT_DIR/" || true
    done
  else
    scp "${REMOTE_HOST}:${ARCHIVES_PATH}/*.tar.gz" "$JOB_OUTPUT_DIR/" || exit 1
  fi

  cd "$JOB_OUTPUT_DIR"
  for archive in *.tar.gz; do
    [[ -e "$archive" ]] || continue
    echo "  Extracting: $archive"
    tar xzf "$archive"
  done
fi

# --- Analysis ---
if [[ "$SKIP_ANALYSIS" == false ]]; then
  ANALYSIS_OUT="$OUTPUT_DIR/lumen_pressure_analysis"
  mkdir -p "$ANALYSIS_OUT"

  MODEL_ARG=""
  [[ -n "$MODEL_FILTER" ]] && MODEL_ARG="--model $MODEL_FILTER"

  FULL_ARG=""
  [[ "$FULL_ANALYSIS" == true ]] && FULL_ARG="--full"

  SEARCH_DIR="$OUTPUT_DIR"
  [[ -n "$JOB_DIR" ]] && SEARCH_DIR="$OUTPUT_DIR/$JOB_DIR"

  echo "  Running lumen pressure analysis..."
  python3 "$SCRIPT_DIR/analyse_and_plot.py" \
    --data-dir "$SEARCH_DIR" \
    $MODEL_ARG \
    $FULL_ARG \
    -o "$ANALYSIS_OUT/plots" \
    2>&1 | tee "$ANALYSIS_OUT/analysis.log"

  echo ""
  echo "  Analysis output: $ANALYSIS_OUT/"
fi

echo "Done."

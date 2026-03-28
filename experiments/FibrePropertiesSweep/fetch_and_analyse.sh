#!/bin/bash
# fetch_and_analyse.sh — Fibre Properties Sweep
set -e

REMOTE_HOST="sv22482@bp1-login.acrc.bris.ac.uk"
REMOTE_BASE="/user/work/sv22482/sim_output"
JOB_DIR=""
OUTPUT_BASE="../../output"
MODEL_FILTER=""
SKIP_SCP=false
SKIP_ANALYSIS=false

while [[ $# -gt 0 ]]; do
  case $1 in
    --job-dir)       JOB_DIR="$2"; shift 2 ;;
    --model)         MODEL_FILTER="$2"; shift 2 ;;
    --skip-scp)      SKIP_SCP=true; shift ;;
    --skip-analysis) SKIP_ANALYSIS=true; shift ;;
    --output-base)   OUTPUT_BASE="$2"; shift 2 ;;
    --help|-h)       grep '^#' "$0" | sed 's/^# //' | sed 's/^#//'; exit 0 ;;
    *)               echo "ERROR: Unknown option: $1"; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$(cd "$SCRIPT_DIR" && cd "$OUTPUT_BASE" && pwd 2>/dev/null || echo "$OUTPUT_BASE")"

if [[ "$SKIP_SCP" == false && -z "$JOB_DIR" ]]; then
  echo "ERROR: --job-dir required unless --skip-scp"; exit 1
fi

echo "=== Fibre Properties Sweep — Fetch & Analyse ==="

if [[ "$SKIP_SCP" == false ]]; then
  JOB_OUTPUT_DIR="$OUTPUT_DIR/$JOB_DIR"
  mkdir -p "$JOB_OUTPUT_DIR"
  scp "${REMOTE_HOST}:${REMOTE_BASE}/${JOB_DIR}/archives/*.tar.gz" "$JOB_OUTPUT_DIR/" || exit 1
  cd "$JOB_OUTPUT_DIR"
  for archive in *.tar.gz; do
    [[ -e "$archive" ]] || continue
    echo "  Extracting: $archive"
    tar xzf "$archive"
  done
fi

if [[ "$SKIP_ANALYSIS" == false ]]; then
  ANALYSIS_OUT="$OUTPUT_DIR/fibre_analysis"
  mkdir -p "$ANALYSIS_OUT"
  MODEL_ARG=""
  [[ -n "$MODEL_FILTER" ]] && MODEL_ARG="--model $MODEL_FILTER"
  SEARCH_DIR="$OUTPUT_DIR"
  [[ -n "$JOB_DIR" ]] && SEARCH_DIR="$OUTPUT_DIR/$JOB_DIR"

  python3 "$SCRIPT_DIR/analyse_and_plot.py" \
    --data-dir "$SEARCH_DIR" $MODEL_ARG \
    -o "$ANALYSIS_OUT/plots" \
    2>&1 | tee "$ANALYSIS_OUT/analysis.log"

  echo "  Output: $ANALYSIS_OUT/"
fi
echo "Done."

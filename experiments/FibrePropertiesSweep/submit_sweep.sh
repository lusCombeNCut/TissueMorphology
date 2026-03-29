#!/bin/bash
# =============================================================================
# submit_sweep.sh — Fibre Properties Sweep
#
# Two-phase Slurm submission for BluePebble HPC (ACRC, University of Bristol)
#
# Sweeps ECM fibre remodeling rate to investigate how fibre dynamics affect
# crypt budding morphology. Higher remodeling rates mean faster fibre
# alignment with cell traction forces.
#
# Since fibreRemodelingRate has no CLI flag, we generate modified JSON configs
# at runtime using Python.
#
# Sweep parameters:
#   fibreRemodelingRate: 0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0
#   10 replicates each = 70 tasks per model
#
# Phase 3 (analyse): Single job — run analysis, generate plots.
#
# Usage:
#   ./submit_sweep.sh node2d
#   ./submit_sweep.sh all2d
#   ./submit_sweep.sh all
#
# To download results locally after jobs complete:
#   scp -r sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/sim_output/analysis_fibre_properties_<TIMESTAMP>/plots/ ./plots/
#   (The exact command with timestamp is printed when you run this script)
# =============================================================================

#SBATCH --account=semt036404
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --array=0-69%20
#SBATCH --output=/user/work/%u/logs/fibre_sweep/slurm-%A_%a.out
#SBATCH --error=/user/work/%u/logs/fibre_sweep/slurm-%A_%a.err

MODEL_TYPE="${1:-node2d}"

case "$MODEL_TYPE" in
    all)    MODELS_TO_RUN=(node2d vertex2d node3d vertex3d) ;;
    all2d)  MODELS_TO_RUN=(node2d vertex2d) ;;
    all3d)  MODELS_TO_RUN=(node3d vertex3d) ;;
    node2d|vertex2d|node3d|vertex3d)  MODELS_TO_RUN=("$MODEL_TYPE") ;;
    *)
        echo "ERROR: Unknown model type '$MODEL_TYPE'"
        exit 1
        ;;
esac

# =============================================================================
# ORCHESTRATION (login node)
# =============================================================================
if [ -z "${SLURM_JOB_ID}" ]; then
    BASE_LOG_DIR="/user/work/$(whoami)/logs/fibre_sweep"
    mkdir -p "${BASE_LOG_DIR}"
    SWEEP_TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

    echo ""
    echo "=== Fibre Properties Sweep: ${MODELS_TO_RUN[*]} ==="
    echo "  Timestamp: ${SWEEP_TIMESTAMP}"
    echo ""

    BUILD_JOB_ID=$(sbatch --parsable \
        --job-name="Build_FibreSweep" \
        --array=0 \
        --cpus-per-task=4 \
        --time=01:00:00 \
        --export=ALL,SWEEP_PHASE=build,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
        --output="${BASE_LOG_DIR}/build_%j.out" \
        --error="${BASE_LOG_DIR}/build_%j.err" \
        "$0" "node2d")

    echo "  Build job: ${BUILD_JOB_ID}"

    SIM_JOB_IDS=()
    RUN_TAGS_CSV=""
    for model in "${MODELS_TO_RUN[@]}"; do
        SIM_JOB_ID=$(sbatch --parsable \
            --job-name="Fibre_${model}" \
            --dependency=afterok:${BUILD_JOB_ID} \
            --export=ALL,SWEEP_PHASE=run,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
            "$0" "${model}")
        SIM_JOB_IDS+=("$SIM_JOB_ID")
        RUN_TAGS_CSV="${RUN_TAGS_CSV:+${RUN_TAGS_CSV},}${SIM_JOB_ID}_${model}_${SWEEP_TIMESTAMP}"
        echo "  ${model} sim: ${SIM_JOB_ID}  (depends on ${BUILD_JOB_ID})"
    done

    # Phase 3: Analysis
    DEP_STR=$(IFS=:; echo "${SIM_JOB_IDS[*]}")
    ANALYSIS_OUTPUT="/user/work/$(whoami)/sim_output/analysis_fibre_properties_${SWEEP_TIMESTAMP}"

    ANALYSE_JOB_ID=$(sbatch --parsable \
        --job-name="Analyse_FibreSweep" \
        --dependency=afterany:${DEP_STR} \
        --array=0 \
        --cpus-per-task=2 \
        --mem-per-cpu=8G \
        --time=02:00:00 \
        --export=ALL,SWEEP_PHASE=analyse,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP},ANALYSIS_RUN_TAGS=${RUN_TAGS_CSV},ANALYSIS_OUTPUT=${ANALYSIS_OUTPUT} \
        --output="${BASE_LOG_DIR}/analyse_%j.out" \
        --error="${BASE_LOG_DIR}/analyse_%j.err" \
        "$0" "node2d")

    echo ""
    echo "  Analysis:  ${ANALYSE_JOB_ID}  (depends on all sims: ${DEP_STR})"
    echo "    Results: ${ANALYSIS_OUTPUT}/plots/"
    echo ""
    echo "  To download plots locally after completion:"
    echo "    scp -r $(whoami)@bp1-login.acrc.bris.ac.uk:${ANALYSIS_OUTPUT}/plots/ ./plots/"
    echo ""
    exit 0
fi

# =============================================================================
# INSIDE SLURM
# =============================================================================
SWEEP_PHASE="${SWEEP_PHASE:-run}"

module load apptainer
module load git

export APPTAINER_CACHEDIR="/user/work/$(whoami)/.apptainer"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

SIF_IMAGE="/user/work/$(whoami)/containers/tissuemorphology.sif"
SOURCE_DIR="/user/work/$(whoami)/TissueMorphology"
BUILD_DIR="/user/work/$(whoami)/chaste_build"
BASE_LOG_DIR="/user/work/$(whoami)/logs/fibre_sweep"
SWEEP_TIMESTAMP="${SWEEP_TIMESTAMP:-$(date +"%Y-%m-%d_%H-%M-%S")}"

mkdir -p "${BASE_LOG_DIR}"

if [ ! -f "${SIF_IMAGE}" ]; then
    echo "ERROR: Container not found at ${SIF_IMAGE}"
    exit 1
fi

# #############################################################################
# BUILD
# #############################################################################
if [ "${SWEEP_PHASE}" = "build" ]; then
    BUILD_LOG_DIR="${BASE_LOG_DIR}/${SLURM_JOB_ID}_build_${SWEEP_TIMESTAMP}"
    mkdir -p "${BUILD_LOG_DIR}"
    LOG_FILE="${BUILD_LOG_DIR}/build.log"
    exec > >(tee -a "${LOG_FILE}")
    exec 2>&1

    echo "============================================"
    echo "  Build (Fibre Properties Sweep)"
    echo "  Job ID: ${SLURM_JOB_ID}  Start: $(date)"
    echo "============================================"

    if [ ! -d "${BUILD_DIR}" ]; then
        mkdir -p "${BUILD_DIR}"
        apptainer exec --bind "${BUILD_DIR}:/mnt" "${SIF_IMAGE}" \
            bash -c "cp -a /home/chaste/build/. /mnt"
    fi

    if [ -d "${SOURCE_DIR}/.git" ]; then
        cd "${SOURCE_DIR}" && git fetch origin && git reset --hard origin/main
    else
        git clone https://github.com/lusCombeNCut/TissueMorphology.git "${SOURCE_DIR}"
        cd "${SOURCE_DIR}"
    fi

    rm -rf "${BUILD_DIR}/_deps/cellml_repo-subbuild"

    apptainer exec \
        --bind "${BUILD_DIR}:/home/chaste/build" \
        --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
        --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
        "${SIF_IMAGE}" \
        bash -c "cd /home/chaste/build && \
                 cmake /home/chaste/src -DChaste_ERROR_ON_WARNING=OFF && \
                 make -j${SLURM_CPUS_PER_TASK} CryptBuddingApp"

    [ $? -ne 0 ] && { echo "Build failed"; exit 1; }

    HOST_APP_PATH="${BUILD_DIR}/projects/TissueMorphology/apps/CryptBuddingApp"
    [ ! -f "${HOST_APP_PATH}" ] && { echo "Binary not found"; exit 1; }

    echo "Build successful: ${HOST_APP_PATH}"
    exit 0
fi

# #############################################################################
# PHASE: ANALYSE
# #############################################################################
if [ "${SWEEP_PHASE}" = "analyse" ]; then

    ANALYSIS_OUTPUT="${ANALYSIS_OUTPUT:-/user/work/$(whoami)/sim_output/analysis_fibre_properties_${SWEEP_TIMESTAMP}}"
    MERGED_DIR="${ANALYSIS_OUTPUT}/data"
    PLOTS_OUT="${ANALYSIS_OUTPUT}/plots"
    mkdir -p "${MERGED_DIR}" "${PLOTS_OUT}"

    LOG_FILE="${ANALYSIS_OUTPUT}/analysis.log"
    exec > >(tee -a "${LOG_FILE}")
    exec 2>&1

    echo "============================================"
    echo "  Fibre Properties Sweep — Analysis Phase"
    echo "  Start Time:   $(date)"
    echo "============================================"

    IFS=',' read -ra RUN_TAGS <<< "${ANALYSIS_RUN_TAGS}"
    SIM_OUTPUT_BASE="/user/work/$(whoami)/sim_output"
    for tag in "${RUN_TAGS[@]}"; do
        [ -d "${SIM_OUTPUT_BASE}/${tag}" ] && \
            ln -sf "${SIM_OUTPUT_BASE}/${tag}" "${MERGED_DIR}/${tag}"
    done

    module load languages/python 2>/dev/null || module load python 2>/dev/null || true

    SCRIPT_DIR="${SOURCE_DIR}/experiments/FibrePropertiesSweep"

    echo "  Running analysis..."
    python3 "${SCRIPT_DIR}/analyse_and_plot.py" \
        --data-dir "${MERGED_DIR}" \
        -o "${PLOTS_OUT}" \
        2>&1 | tee "${PLOTS_OUT}/plot_generation.log" || true

    echo ""
    echo "============================================"
    echo "  Analysis Complete — $(date)"
    echo "  Plots: ${PLOTS_OUT}/"
    echo "============================================"
    exit 0
fi

# #############################################################################
# RUN
# #############################################################################

FIBRE_REMODEL_VALUES=(0.0 0.01 0.05 0.1 0.2 0.5 1.0)
NUM_PARAMS=${#FIBRE_REMODEL_VALUES[@]}
NUM_REPLICATES=10

PARAM_INDEX=$((SLURM_ARRAY_TASK_ID / NUM_REPLICATES))
RUN_NUMBER=$((SLURM_ARRAY_TASK_ID % NUM_REPLICATES))

if [ "$PARAM_INDEX" -ge "$NUM_PARAMS" ]; then
    echo "Task $SLURM_ARRAY_TASK_ID exceeds parameter space. Exiting."
    exit 0
fi

FIBRE_RATE=${FIBRE_REMODEL_VALUES[$PARAM_INDEX]}

case "$MODEL_TYPE" in
    node2d)   JOB_PREFIX="N2"; CONFIG_FILENAME="params-Node2d.json"   ;;
    vertex2d) JOB_PREFIX="V2"; CONFIG_FILENAME="params-Vertex2d.json" ;;
    node3d)   JOB_PREFIX="N3"; CONFIG_FILENAME="params-Node3d.json"   ;;
    vertex3d) JOB_PREFIX="V3"; CONFIG_FILENAME="params-Vertex3d.json" ;;
esac

JOB_NAME="${JOB_PREFIX}F${FIBRE_RATE}R${RUN_NUMBER}"
scontrol update JobId="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" JobName="${JOB_NAME}"

RUN_TAG="${SLURM_ARRAY_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}"
LOG_DIR="${BASE_LOG_DIR}/${RUN_TAG}"
OUTPUT_BASE="/user/work/$(whoami)/sim_output/${RUN_TAG}"
OUTPUT_DIR="${OUTPUT_BASE}/f${FIBRE_RATE}_r${RUN_NUMBER}"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/f${FIBRE_RATE}_r${RUN_NUMBER}.log"
exec > >(tee -a "${LOG_FILE}")
exec 2>&1

# --- Generate modified config with fibre remodeling rate ---
HOST_CONFIG_PATH="${SOURCE_DIR}/apps/${CONFIG_FILENAME}"
MODIFIED_CONFIG="${OUTPUT_DIR}/params_modified.json"

python3 -c "
import json, sys
with open('${HOST_CONFIG_PATH}') as f:
    cfg = json.load(f)
cfg['forces']['GhostNodeECM']['fibreRemodelingRate'] = ${FIBRE_RATE}
cfg['simulation']['runNumber'] = ${RUN_NUMBER}
with open('${MODIFIED_CONFIG}', 'w') as f:
    json.dump(cfg, f, indent=2)
print('Modified config: fibreRemodelingRate=${FIBRE_RATE}')
" || { echo "ERROR: Failed to generate modified config"; exit 1; }

# Copy modified config into a path visible inside the container
CONTAINER_CONFIG="/home/chaste/output/params_modified.json"

echo "============================================"
echo "  Fibre Properties Sweep"
echo "  Model:              ${MODEL_TYPE}"
echo "  Fibre Remodel Rate: ${FIBRE_RATE}"
echo "  Replicate:          ${RUN_NUMBER}"
echo "  Start Time:         $(date)"
echo "============================================"

APP_PATH="/home/chaste/build/projects/TissueMorphology/apps/CryptBuddingApp"
HOST_APP_PATH="${BUILD_DIR}/projects/TissueMorphology/apps/CryptBuddingApp"
[ ! -f "${HOST_APP_PATH}" ] && { echo "Binary not found"; exit 1; }

apptainer exec \
    --bind "${BUILD_DIR}:/home/chaste/build" \
    --bind "${SOURCE_DIR}:/home/chaste/src/projects/TissueMorphology" \
    --bind "${OUTPUT_DIR}:/home/chaste/output" \
    --env OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
    --env CHASTE_TEST_OUTPUT=/home/chaste/output \
    "${SIF_IMAGE}" \
    "${APP_PATH}" \
        -config "${CONTAINER_CONFIG}" \
        -run "${RUN_NUMBER}"

EXIT_CODE=$?

ARCHIVE_DIR="${OUTPUT_BASE}/archives"
mkdir -p "${ARCHIVE_DIR}"
ARCHIVE_PATH="${ARCHIVE_DIR}/f${FIBRE_RATE}_r${RUN_NUMBER}.tar.gz"

if [ -d "${OUTPUT_DIR}" ] && [ "$(ls -A ${OUTPUT_DIR})" ]; then
    tar -czf "${ARCHIVE_PATH}" -C "$(dirname ${OUTPUT_DIR})" "$(basename ${OUTPUT_DIR})"
fi

echo "============================================"
echo "  Job Complete (exit: ${EXIT_CODE})  End: $(date)"
echo "============================================"

exit ${EXIT_CODE}

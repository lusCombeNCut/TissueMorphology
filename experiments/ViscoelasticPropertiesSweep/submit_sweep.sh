#!/bin/bash
# =============================================================================
# submit_sweep.sh — Viscoelastic Properties Sweep
#
# Two-phase Slurm submission for BluePebble HPC (ACRC, University of Bristol)
#
# Sweeps the viscoelastic ECM relaxation time (tau) to investigate how ECM
# stress relaxation kinetics affect crypt budding. Longer tau = stiffer
# initial response; shorter tau = more fluid-like ECM.
#
# Also sweeps ghostRelaxationModulus (E_1) in a secondary experiment to
# separate the effects of relaxation rate vs relaxation magnitude.
#
# Primary sweep: ghostRelaxationTime (tau)
#   Values: 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0
#   10 replicates = 70 tasks per model
#
# Usage:
#   ./submit_sweep.sh node2d
#   ./submit_sweep.sh all2d
#   ./submit_sweep.sh all
# =============================================================================

#SBATCH --account=semt036404
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --array=0-69%20
#SBATCH --output=/user/work/%u/logs/viscoelastic_sweep/slurm-%A_%a.out
#SBATCH --error=/user/work/%u/logs/viscoelastic_sweep/slurm-%A_%a.err

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
# ORCHESTRATION
# =============================================================================
if [ -z "${SLURM_JOB_ID}" ]; then
    BASE_LOG_DIR="/user/work/$(whoami)/logs/viscoelastic_sweep"
    mkdir -p "${BASE_LOG_DIR}"
    SWEEP_TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

    echo ""
    echo "=== Viscoelastic Properties Sweep: ${MODELS_TO_RUN[*]} ==="
    echo "  Timestamp: ${SWEEP_TIMESTAMP}"
    echo ""

    BUILD_JOB_ID=$(sbatch --parsable \
        --job-name="Build_ViscoSweep" \
        --array=0 \
        --cpus-per-task=4 \
        --time=01:00:00 \
        --export=ALL,SWEEP_PHASE=build,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
        --output="${BASE_LOG_DIR}/build_%j.out" \
        --error="${BASE_LOG_DIR}/build_%j.err" \
        "$0" "node2d")

    echo "  Build job: ${BUILD_JOB_ID}"

    for model in "${MODELS_TO_RUN[@]}"; do
        SIM_JOB_ID=$(sbatch --parsable \
            --job-name="Visco_${model}" \
            --dependency=afterok:${BUILD_JOB_ID} \
            --export=ALL,SWEEP_PHASE=run,SWEEP_TIMESTAMP=${SWEEP_TIMESTAMP} \
            "$0" "${model}")
        echo "  ${model} sim: ${SIM_JOB_ID}  (depends on ${BUILD_JOB_ID})"
    done
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
BASE_LOG_DIR="/user/work/$(whoami)/logs/viscoelastic_sweep"
SWEEP_TIMESTAMP="${SWEEP_TIMESTAMP:-$(date +"%Y-%m-%d_%H-%M-%S")}"

mkdir -p "${BASE_LOG_DIR}"
[ ! -f "${SIF_IMAGE}" ] && { echo "ERROR: Container not found"; exit 1; }

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
    echo "  Build (Viscoelastic Sweep)  Start: $(date)"
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
    echo "Build successful"
    exit 0
fi

# #############################################################################
# RUN
# #############################################################################

TAU_VALUES=(0.1 0.5 1.0 2.0 5.0 10.0 50.0)
NUM_PARAMS=${#TAU_VALUES[@]}
NUM_REPLICATES=10

PARAM_INDEX=$((SLURM_ARRAY_TASK_ID / NUM_REPLICATES))
RUN_NUMBER=$((SLURM_ARRAY_TASK_ID % NUM_REPLICATES))

if [ "$PARAM_INDEX" -ge "$NUM_PARAMS" ]; then
    echo "Task $SLURM_ARRAY_TASK_ID exceeds parameter space. Exiting."
    exit 0
fi

TAU=${TAU_VALUES[$PARAM_INDEX]}

case "$MODEL_TYPE" in
    node2d)   JOB_PREFIX="N2"; CONFIG_FILENAME="params-Node2d.json"   ;;
    vertex2d) JOB_PREFIX="V2"; CONFIG_FILENAME="params-Vertex2d.json" ;;
    node3d)   JOB_PREFIX="N3"; CONFIG_FILENAME="params-Node3d.json"   ;;
    vertex3d) JOB_PREFIX="V3"; CONFIG_FILENAME="params-Vertex3d.json" ;;
esac

JOB_NAME="${JOB_PREFIX}T${TAU}R${RUN_NUMBER}"
scontrol update JobId="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" JobName="${JOB_NAME}"

RUN_TAG="${SLURM_ARRAY_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}"
LOG_DIR="${BASE_LOG_DIR}/${RUN_TAG}"
OUTPUT_BASE="/user/work/$(whoami)/sim_output/${RUN_TAG}"
OUTPUT_DIR="${OUTPUT_BASE}/tau${TAU}_r${RUN_NUMBER}"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/tau${TAU}_r${RUN_NUMBER}.log"
exec > >(tee -a "${LOG_FILE}")
exec 2>&1

# --- Generate modified config ---
HOST_CONFIG_PATH="${SOURCE_DIR}/apps/${CONFIG_FILENAME}"
MODIFIED_CONFIG="${OUTPUT_DIR}/params_modified.json"

python3 -c "
import json
with open('${HOST_CONFIG_PATH}') as f:
    cfg = json.load(f)
# Ensure viscoelastic ECM is enabled
cfg['features']['enableViscoelasticECM'] = True
cfg['features']['enableGhostNodeECM'] = True
# Set relaxation time
cfg['forces']['GhostNodeECM']['ViscoelasticECM']['ghostRelaxationTime'] = ${TAU}
cfg['simulation']['runNumber'] = ${RUN_NUMBER}
with open('${MODIFIED_CONFIG}', 'w') as f:
    json.dump(cfg, f, indent=2)
print('Modified config: ghostRelaxationTime=${TAU}')
" || { echo "ERROR: Config generation failed"; exit 1; }

CONTAINER_CONFIG="/home/chaste/output/params_modified.json"

echo "============================================"
echo "  Viscoelastic Properties Sweep"
echo "  Model:             ${MODEL_TYPE}"
echo "  Relaxation Time:   ${TAU}"
echo "  Replicate:         ${RUN_NUMBER}"
echo "  Start Time:        $(date)"
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
ARCHIVE_PATH="${ARCHIVE_DIR}/tau${TAU}_r${RUN_NUMBER}.tar.gz"

if [ -d "${OUTPUT_DIR}" ] && [ "$(ls -A ${OUTPUT_DIR})" ]; then
    tar -czf "${ARCHIVE_PATH}" -C "$(dirname ${OUTPUT_DIR})" "$(basename ${OUTPUT_DIR})"
fi

echo "============================================"
echo "  Job Complete (exit: ${EXIT_CODE})  End: $(date)"
echo "============================================"

exit ${EXIT_CODE}

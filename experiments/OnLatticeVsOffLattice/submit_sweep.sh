#!/bin/bash
# =============================================================================
# submit_sweep.sh — On-Lattice vs Off-Lattice ECM Comparison
#
# Two-phase Slurm submission for BluePebble HPC (ACRC, University of Bristol)
#
# Compares two ECM representations at matched stiffness levels:
#   A) Grid-based ECMConfinementForce (on-lattice):
#      enableEcmConfinement=true, enableGhostNodeECM=false
#   B) Particle-based GhostNodeEcmForce (off-lattice):
#      enableGhostNodeECM=true (default)
#
# Sweep parameters:
#   ECM type:  on-lattice, off-lattice (2 configs)
#   Stiffness: 0.5, 2.0, 5.0, 10.0, 20.0, 50.0 (6 values)
#   Replicates: 10
#   Total: 2 x 6 x 10 = 120 tasks per model
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
#SBATCH --array=0-119%20
#SBATCH --output=/user/work/%u/logs/lattice_comparison/slurm-%A_%a.out
#SBATCH --error=/user/work/%u/logs/lattice_comparison/slurm-%A_%a.err

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
    BASE_LOG_DIR="/user/work/$(whoami)/logs/lattice_comparison"
    mkdir -p "${BASE_LOG_DIR}"
    SWEEP_TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

    echo ""
    echo "=== On-Lattice vs Off-Lattice Comparison: ${MODELS_TO_RUN[*]} ==="
    echo "  Timestamp: ${SWEEP_TIMESTAMP}"
    echo "  Tasks per model: 120 (2 ECM types x 6 stiffness x 10 replicates)"
    echo ""

    BUILD_JOB_ID=$(sbatch --parsable \
        --job-name="Build_Lattice" \
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
            --job-name="Lattice_${model}" \
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
BASE_LOG_DIR="/user/work/$(whoami)/logs/lattice_comparison"
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
    echo "  Build (Lattice Comparison)  Start: $(date)"
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

# Parameter space: 2 ECM types x 6 stiffness values x 10 replicates
ECM_TYPES=("on-lattice" "off-lattice")
STIFFNESS_VALUES=(0.5 2.0 5.0 10.0 20.0 50.0)
NUM_ECM_TYPES=${#ECM_TYPES[@]}
NUM_STIFFNESS=${#STIFFNESS_VALUES[@]}
NUM_REPLICATES=10

# Decode: task_id = ecm_index * (stiffness * replicates) + stiff_index * replicates + run
TASKS_PER_ECM=$((NUM_STIFFNESS * NUM_REPLICATES))
ECM_INDEX=$((SLURM_ARRAY_TASK_ID / TASKS_PER_ECM))
REMAINDER=$((SLURM_ARRAY_TASK_ID % TASKS_PER_ECM))
STIFFNESS_INDEX=$((REMAINDER / NUM_REPLICATES))
RUN_NUMBER=$((REMAINDER % NUM_REPLICATES))

if [ "$ECM_INDEX" -ge "$NUM_ECM_TYPES" ] || [ "$STIFFNESS_INDEX" -ge "$NUM_STIFFNESS" ]; then
    echo "Task $SLURM_ARRAY_TASK_ID exceeds parameter space. Exiting."
    exit 0
fi

ECM_TYPE=${ECM_TYPES[$ECM_INDEX]}
ECM_STIFFNESS=${STIFFNESS_VALUES[$STIFFNESS_INDEX]}

case "$MODEL_TYPE" in
    node2d)   JOB_PREFIX="N2"; CONFIG_FILENAME="params-Node2d.json"   ;;
    vertex2d) JOB_PREFIX="V2"; CONFIG_FILENAME="params-Vertex2d.json" ;;
    node3d)   JOB_PREFIX="N3"; CONFIG_FILENAME="params-Node3d.json"   ;;
    vertex3d) JOB_PREFIX="V3"; CONFIG_FILENAME="params-Vertex3d.json" ;;
esac

# Short ECM type label for job name
ECM_LABEL="ON"
[ "$ECM_TYPE" = "off-lattice" ] && ECM_LABEL="OFF"

JOB_NAME="${JOB_PREFIX}${ECM_LABEL}S${ECM_STIFFNESS}R${RUN_NUMBER}"
scontrol update JobId="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" JobName="${JOB_NAME}"

RUN_TAG="${SLURM_ARRAY_JOB_ID}_${MODEL_TYPE}_${SWEEP_TIMESTAMP}"
LOG_DIR="${BASE_LOG_DIR}/${RUN_TAG}"
OUTPUT_BASE="/user/work/$(whoami)/sim_output/${RUN_TAG}"
OUTPUT_DIR="${OUTPUT_BASE}/${ECM_TYPE}_s${ECM_STIFFNESS}_r${RUN_NUMBER}"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/${ECM_TYPE}_s${ECM_STIFFNESS}_r${RUN_NUMBER}.log"
exec > >(tee -a "${LOG_FILE}")
exec 2>&1

# --- Generate modified config ---
HOST_CONFIG_PATH="${SOURCE_DIR}/apps/${CONFIG_FILENAME}"
MODIFIED_CONFIG="${OUTPUT_DIR}/params_modified.json"

python3 -c "
import json
with open('${HOST_CONFIG_PATH}') as f:
    cfg = json.load(f)

ecm_type = '${ECM_TYPE}'
stiffness = ${ECM_STIFFNESS}

if ecm_type == 'on-lattice':
    # Grid-based ECM only
    cfg['features']['enableEcmConfinement'] = True
    cfg['features']['enableGhostNodeECM'] = False
    cfg['features']['enableViscoelasticECM'] = False
    cfg['forces']['ECMConfinementForce']['ecmStiffness'] = stiffness
else:
    # Particle-based ghost node ECM
    cfg['features']['enableEcmConfinement'] = False
    cfg['features']['enableGhostNodeECM'] = True
    cfg['features']['enableViscoelasticECM'] = True
    cfg['forces']['ECMConfinementForce']['ecmStiffness'] = stiffness
    cfg['forces']['GhostNodeECM']['ViscoelasticECM']['ghostRelaxedStiffness'] = stiffness

# Tag the ECM type in params for analysis identification
cfg['simulation']['ecmType'] = ecm_type
cfg['simulation']['runNumber'] = ${RUN_NUMBER}

with open('${MODIFIED_CONFIG}', 'w') as f:
    json.dump(cfg, f, indent=2)
print('Config: ecmType=${ECM_TYPE}, stiffness=${ECM_STIFFNESS}')
" || { echo "ERROR: Config generation failed"; exit 1; }

CONTAINER_CONFIG="/home/chaste/output/params_modified.json"

echo "============================================"
echo "  On-Lattice vs Off-Lattice Comparison"
echo "  Model:      ${MODEL_TYPE}"
echo "  ECM Type:   ${ECM_TYPE}"
echo "  Stiffness:  ${ECM_STIFFNESS}"
echo "  Replicate:  ${RUN_NUMBER}"
echo "  Start:      $(date)"
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
        -stiffness "${ECM_STIFFNESS}" \
        -run "${RUN_NUMBER}"

EXIT_CODE=$?

ARCHIVE_DIR="${OUTPUT_BASE}/archives"
mkdir -p "${ARCHIVE_DIR}"
ARCHIVE_PATH="${ARCHIVE_DIR}/${ECM_TYPE}_s${ECM_STIFFNESS}_r${RUN_NUMBER}.tar.gz"

if [ -d "${OUTPUT_DIR}" ] && [ "$(ls -A ${OUTPUT_DIR})" ]; then
    tar -czf "${ARCHIVE_PATH}" -C "$(dirname ${OUTPUT_DIR})" "$(basename ${OUTPUT_DIR})"
fi

echo "============================================"
echo "  Job Complete (exit: ${EXIT_CODE})  End: $(date)"
echo "============================================"

exit ${EXIT_CODE}

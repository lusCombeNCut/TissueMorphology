# =============================================================================
# Dockerfile: Chaste + OrganoidChaste + TissueMorphology
#
# Builds a self-contained image with the 3D vertex crypt organoid simulation
# (Test3dVertexCryptOrganoid) ready to run on HPC via Apptainer/Singularity.
#
# Base: chaste/base (Ubuntu + all Chaste dependencies pre-installed)
# Projects: OrganoidChaste (3D vertex model) + TissueMorphology (ECM forces)
#
# Build (from TissueMorphology/ dir):
#   docker build -f docker/Dockerfile -t ghcr.io/luscombencut/tissuemorphology:latest ../../
# Push:   docker push ghcr.io/luscombencut/tissuemorphology:latest
# HPC:    apptainer pull docker://ghcr.io/luscombencut/tissuemorphology:latest
# =============================================================================

FROM chaste/base:latest

# ---------- System dependencies ----------
USER root

# CGAL is required by OrganoidChaste (3D Voronoi sphere mesh generation)
RUN apt-get update && \
    apt-get install -y --no-install-recommends libcgal-dev && \
    rm -rf /var/lib/apt/lists/*

# ---------- Environment ----------
USER chaste

ENV CHASTE_DIR=/home/chaste
ENV CHASTE_SOURCE_DIR=${CHASTE_DIR}/src
ENV CHASTE_BUILD_DIR=${CHASTE_DIR}/build
ENV CHASTE_TEST_OUTPUT=${CHASTE_DIR}/output

WORKDIR ${CHASTE_DIR}

# ---------- Chaste source (develop branch, shallow clone) ----------
RUN git clone --depth 1 -b develop https://github.com/Chaste/Chaste.git src

# ---------- Copy custom projects ----------
# Build context = Chaste/projects/ (set via -f docker/Dockerfile ../../)
# OrganoidChaste: finite-thickness 3D vertex model (Drozdowski & Schwarz)
COPY --chown=chaste:chaste OrganoidChaste/ ${CHASTE_SOURCE_DIR}/projects/OrganoidChaste/

# TissueMorphology: ECM forces, basement membrane, dynamic ECM field
COPY --chown=chaste:chaste TissueMorphology/ ${CHASTE_SOURCE_DIR}/projects/TissueMorphology/

# ---------- Build (Release mode for HPC performance) ----------
RUN mkdir -p ${CHASTE_BUILD_DIR} ${CHASTE_TEST_OUTPUT} && \
    cd ${CHASTE_BUILD_DIR} && \
    cmake ${CHASTE_SOURCE_DIR} \
        -DCMAKE_BUILD_TYPE=Release \
        -DChaste_ERROR_ON_WARNING=OFF \
        -DChaste_UPDATE_PROVENANCE=OFF && \
    make -j$(nproc) project_OrganoidChaste && \
    make -j$(nproc) project_TissueMorphology

# ---------- Default working directory ----------
WORKDIR ${CHASTE_BUILD_DIR}

# Default command: run the vertex crypt organoid test
CMD ["ctest", "-R", "Test3dVertexCryptOrganoid", "-V", "--output-on-failure"]

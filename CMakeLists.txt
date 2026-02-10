# Copyright (c) 2005-2023, University of Oxford.
# All rights reserved.
# 
# University of Oxford means the Chancellor, Masters and Scholars of the
# University of Oxford, having an administrative office at Wellington
# Square, Oxford OX1 2JD, UK.
# 
# This file is part of Chaste.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * Neither the name of the University of Oxford nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# Example CMake build configuration file for a Chaste user project.
# See also https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/CmakeBuildGuide for full documentation.

# You can set which Chaste components (or other projects) your project depends on by editing the
# find_package() call for Chaste.  E.g.
#   find_package(Chaste COMPONENTS cell_based)
# for a project just using the cell_based component (and its dependencies), or
#   find_package(Chaste COMPONENTS heart lung)
# for a project combining heart & lung simulations.
#
# Note that the order in which components are specified does not matter.

# Here we just depend on core components (nothing application-specific).
find_package(Chaste COMPONENTS cell_based crypt)

# OrganoidChaste provides the 3D vertex model (MonolayerVertexBasedCellPopulation,
# FiniteThicknessSimulation3d, SurfaceTensionForce, etc.)
# Required for the vertex-based organoid simulation (Test3dVertexCryptOrganoid)
find_package(CGAL)

# Add OrganoidChaste include directories BEFORE chaste_do_project()
# so they are available when test targets are created
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_mesh")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/forces")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/forces/sub_forces")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/cell_writers")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/face_writers")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/population_writers")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/boundary_conditions")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/moving_boundary_conditions")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/moving_boundary_conditions/move_profiles")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/division_rule")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_population/cell_property_mutations")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_simulation")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_simulation/modifiers")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste/src/monolayer_common")
endif()

# Alternatively, to specify a Chaste installation directory use a line like that below.
# This is needed if your project is not contained in the projects folder within a Chaste source tree.
#find_package(Chaste COMPONENTS heart crypt PATHS /path/to/chaste-install NO_DEFAULT_PATH)

# Change the project name in the line below to match the folder this file is in,
# i.e. the name of your project.
chaste_do_project(TissueMorphology)

# Link against OrganoidChaste library AFTER chaste_do_project()
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../OrganoidChaste")
    if (BUILD_SHARED_LIBS)
        target_link_libraries(chaste_project_TissueMorphology LINK_PUBLIC chaste_project_OrganoidChaste)
    endif()
endif()

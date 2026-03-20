/*
 * SurfaceMeshWriter.hpp
 *
 * A simulation modifier that writes the 3D organoid surface as VTK
 * PolyData (.vtp) files with triangulated surface connectivity.
 *
 * Uses the SurfaceTopologyTracker's edge graph to extract triangles
 * (faces where all three pairs are connected in the adjacency graph).
 *
 * Each output file contains:
 *   - Node positions (Points)
 *   - Triangular faces (Polys)
 *   - Per-node cell data (cell type, polarity, generation, etc.)
 *
 * A companion .pvd collection file is written so the whole time series
 * can be loaded in ParaView with a single File -> Open.
 */
#ifndef SURFACEMESHWRITER_HPP_
#define SURFACEMESHWRITER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "SurfaceTopologyTracker.hpp"
#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include "SimProfiler.hpp"
#include "AbstractCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <algorithm>

template<unsigned DIM>
class SurfaceMeshWriter : public AbstractCellBasedSimulationModifier<DIM>
{
private:
    /** Pointer to the surface topology tracker (not owned) */
    SurfaceTopologyTracker<DIM>* mpSurfaceTracker;

    /** Absolute path to the output directory */
    std::string mOutputDirectory;

    /** Sampling timestep multiple (must match the simulation's value) */
    unsigned mSamplingTimestepMultiple;

    /** PVD collection file stream */
    out_stream mpPvdFile;

public:
    /**
     * Constructor.
     *
     * @param pSurfaceTracker  Raw pointer to a SurfaceTopologyTracker that is
     *                         already registered with the simulation.
     * @param samplingTimestepMultiple  How often to write (in timesteps).
     */
    SurfaceMeshWriter(SurfaceTopologyTracker<DIM>* pSurfaceTracker,
                      unsigned samplingTimestepMultiple = 12)
        : AbstractCellBasedSimulationModifier<DIM>(),
          mpSurfaceTracker(pSurfaceTracker),
          mOutputDirectory(""),
          mSamplingTimestepMultiple(samplingTimestepMultiple)
    {
    }

    virtual ~SurfaceMeshWriter()
    {
        if (mpPvdFile)
        {
            *mpPvdFile << "  </Collection>\n";
            *mpPvdFile << "</VTKFile>\n";
            mpPvdFile->close();
        }
    }

    /* ------------------------------------------------------------------ */
    /*  SimulationModifier interface                                       */
    /* ------------------------------------------------------------------ */

    void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                    std::string outputDirectory) override
    {
        OutputFileHandler handler(outputDirectory, false);
        mOutputDirectory = handler.GetOutputDirectoryFullPath();

        // Create the PVD collection file
        mpPvdFile = handler.OpenOutputFile("surface_mesh.pvd");
        *mpPvdFile << "<?xml version=\"1.0\"?>\n";
        *mpPvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        *mpPvdFile << "  <Collection>\n";

        // Write the initial frame (timestep 0)
        WriteFrame(rCellPopulation, 0);
    }

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation) override
    {
        ScopedTimer timer("SurfaceMeshWriter");
        unsigned ts = SimulationTime::Instance()->GetTimeStepsElapsed();
        if (ts % mSamplingTimestepMultiple == 0)
        {
            WriteFrame(rCellPopulation, ts);
        }
    }

    void OutputSimulationModifierParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<SamplingTimestepMultiple>"
                     << mSamplingTimestepMultiple
                     << "</SamplingTimestepMultiple>\n";
        AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
    }

private:
    /* ------------------------------------------------------------------ */
    /*  Core writing logic                                                 */
    /* ------------------------------------------------------------------ */

    /**
     * Write a single .vtp frame and append an entry to the PVD file.
     *
     * @param rCellPopulation  The current cell population.
     * @param timestep         The current timestep number (used in filenames
     *                         and PVD timestamps so it matches results.pvd).
     */
    void WriteFrame(AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                    unsigned timestep)
    {
        // --- 1. Collect all live cells and build index mapping --------------
        std::vector<unsigned> cellIndices;
        std::map<unsigned, unsigned> indexToLocal;  // cell index -> local VTK index
        std::map<unsigned, CellPtr> cellMap;

        for (auto it = rCellPopulation.Begin(); it != rCellPopulation.End(); ++it)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*it);
            indexToLocal[idx] = cellIndices.size();
            cellIndices.push_back(idx);
            cellMap[idx] = *it;
        }

        unsigned nPts = cellIndices.size();
        if (nPts < 3) return;

        // --- 2. Build positions and per-node data ---------------------------
        std::vector<c_vector<double, 3>> positions(nPts);
        std::vector<double> cellTypeId(nPts, 0.0);
        std::vector<double> cellGeneration(nPts, 0.0);
        std::vector<double> cellAge(nPts, 0.0);

        for (unsigned k = 0; k < nPts; k++)
        {
            unsigned idx = cellIndices[k];
            c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(cellMap[idx]);
            positions[k][0] = pos[0];
            positions[k][1] = pos[1];
            positions[k][2] = (DIM >= 3) ? pos[2] : 0.0;

            CellPtr pCell = cellMap[idx];
            if (pCell->GetCellProliferativeType()->template IsType<StemCellProliferativeType>())
                cellTypeId[k] = 0.0;
            else if (pCell->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
                cellTypeId[k] = 1.0;
            else
                cellTypeId[k] = 2.0;

            cellAge[k] = pCell->GetAge();

            // Try to get generation from cell data if available
            if (pCell->GetCellData()->HasItem("generation"))
                cellGeneration[k] = pCell->GetCellData()->GetItem("generation");
        }

        // --- 3. Extract triangles from the adjacency graph ------------------
        // A triangle exists if three cells A, B, C are all pairwise connected
        std::vector<std::tuple<unsigned, unsigned, unsigned>> triangles;
        std::set<std::tuple<unsigned, unsigned, unsigned>> triangleSet;

        for (unsigned i = 0; i < nPts; i++)
        {
            unsigned idxA = cellIndices[i];
            const std::set<unsigned>& neighborsA = mpSurfaceTracker->GetNeighbors(idxA);

            for (unsigned nbB : neighborsA)
            {
                if (indexToLocal.find(nbB) == indexToLocal.end()) continue;
                unsigned localB = indexToLocal[nbB];
                if (localB <= i) continue;  // avoid duplicates

                const std::set<unsigned>& neighborsB = mpSurfaceTracker->GetNeighbors(nbB);

                // Find common neighbors of A and B (these form triangles)
                for (unsigned nbC : neighborsA)
                {
                    if (indexToLocal.find(nbC) == indexToLocal.end()) continue;
                    unsigned localC = indexToLocal[nbC];
                    if (localC <= localB) continue;  // avoid duplicates

                    // Check if B and C are also connected
                    if (neighborsB.count(nbC) > 0)
                    {
                        // Found a triangle: i, localB, localC
                        // Sort indices to create canonical form
                        std::vector<unsigned> tri = {i, localB, localC};
                        std::sort(tri.begin(), tri.end());
                        auto triTuple = std::make_tuple(tri[0], tri[1], tri[2]);

                        if (triangleSet.find(triTuple) == triangleSet.end())
                        {
                            triangleSet.insert(triTuple);
                            triangles.push_back(std::make_tuple(i, localB, localC));
                        }
                    }
                }
            }
        }

        unsigned nTris = triangles.size();

        // --- 4. Write the .vtp file -----------------------------------------
        std::ostringstream fname;
        fname << "surface_" << timestep << ".vtp";
        std::string vtpFilename = fname.str();
        std::string fullPath = mOutputDirectory + vtpFilename;

        std::ofstream vtp(fullPath.c_str());
        if (!vtp.is_open())
        {
            EXCEPTION("SurfaceMeshWriter: cannot open " + fullPath);
        }

        vtp << "<?xml version=\"1.0\"?>\n";
        vtp << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        vtp << "  <PolyData>\n";
        vtp << "    <Piece NumberOfPoints=\"" << nPts
            << "\" NumberOfVerts=\"0\""
            << " NumberOfLines=\"0\""
            << " NumberOfStrips=\"0\""
            << " NumberOfPolys=\"" << nTris << "\">\n";

        // --- Points ---------------------------------------------------------
        vtp << "      <Points>\n";
        vtp << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (unsigned k = 0; k < nPts; k++)
        {
            vtp << "          " << positions[k][0] << " "
                                << positions[k][1] << " "
                                << positions[k][2] << "\n";
        }
        vtp << "        </DataArray>\n";
        vtp << "      </Points>\n";

        // --- Polys (triangles) - connectivity + offsets ---------------------
        vtp << "      <Polys>\n";

        // connectivity: three vertex indices per triangle
        vtp << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        vtp << "          ";
        for (const auto& tri : triangles)
        {
            vtp << std::get<0>(tri) << " "
                << std::get<1>(tri) << " "
                << std::get<2>(tri) << " ";
        }
        vtp << "\n        </DataArray>\n";

        // offsets: cumulative count of indices per polygon (each triangle has 3 vertices)
        vtp << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        vtp << "          ";
        for (unsigned k = 0; k < nTris; k++)
        {
            vtp << 3 * (k + 1) << " ";
        }
        vtp << "\n        </DataArray>\n";

        vtp << "      </Polys>\n";

        // --- PointData (per-node scalars) ------------------------------------
        vtp << "      <PointData Scalars=\"cell_type\">\n";

        vtp << "        <DataArray type=\"Float64\" Name=\"cell_type\" format=\"ascii\">\n";
        vtp << "          ";
        for (unsigned k = 0; k < nPts; k++)
        {
            vtp << cellTypeId[k] << " ";
        }
        vtp << "\n        </DataArray>\n";

        vtp << "        <DataArray type=\"Float64\" Name=\"cell_age\" format=\"ascii\">\n";
        vtp << "          ";
        for (unsigned k = 0; k < nPts; k++)
        {
            vtp << cellAge[k] << " ";
        }
        vtp << "\n        </DataArray>\n";

        vtp << "        <DataArray type=\"Float64\" Name=\"generation\" format=\"ascii\">\n";
        vtp << "          ";
        for (unsigned k = 0; k < nPts; k++)
        {
            vtp << cellGeneration[k] << " ";
        }
        vtp << "\n        </DataArray>\n";

        vtp << "      </PointData>\n";

        // Close tags
        vtp << "    </Piece>\n";
        vtp << "  </PolyData>\n";
        vtp << "</VTKFile>\n";
        vtp.close();

        // --- 5. Append to PVD -----------------------------------------------
        if (mpPvdFile)
        {
            *mpPvdFile << "    <DataSet timestep=\"" << timestep
                       << "\" group=\"\" part=\"0\" file=\""
                       << vtpFilename << "\"/>\n";
        }
    }
};

#endif // SURFACEMESHWRITER_HPP_

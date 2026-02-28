/*
 * RingOutlineWriter.hpp
 *
 * A simulation modifier that writes the organoid ring outline as VTK
 * PolyData (.vtp) files with line-segment connectivity.  ParaView can
 * render this as a closed polygon / wireframe overlay on top of the
 * ECM grid or point cloud.
 *
 * Each output file contains:
 *   - Node positions (Points)
 *   - Line segments connecting ring-adjacent cells (Lines/Polys)
 *   - Per-node cell data (cell type, polarity, etc.)
 *
 * A companion .pvd collection file is written so the whole time series
 * can be loaded in ParaView with a single File -> Open.
 */
#ifndef RINGOUTLINEWRITER_HPP_
#define RINGOUTLINEWRITER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "RingTopologyTracker.hpp"
#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include "AbstractCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>

template<unsigned DIM>
class RingOutlineWriter : public AbstractCellBasedSimulationModifier<DIM>
{
private:
    /** Pointer to the ring topology tracker (not owned) */
    RingTopologyTracker<DIM>* mpRingTracker;

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
     * @param pRingTracker  Raw pointer to a RingTopologyTracker that is
     *                      already registered with the simulation.
     * @param samplingTimestepMultiple  How often to write (in timesteps).
     */
    RingOutlineWriter(RingTopologyTracker<DIM>* pRingTracker,
                      unsigned samplingTimestepMultiple = 12)
        : AbstractCellBasedSimulationModifier<DIM>(),
          mpRingTracker(pRingTracker),
          mOutputDirectory(""),
          mSamplingTimestepMultiple(samplingTimestepMultiple)
    {
    }

    virtual ~RingOutlineWriter()
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
        mpPvdFile = handler.OpenOutputFile("outline_results.pvd");
        *mpPvdFile << "<?xml version=\"1.0\"?>\n";
        *mpPvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        *mpPvdFile << "  <Collection>\n";

        // Write the initial frame (timestep 0)
        WriteFrame(rCellPopulation, 0);
    }

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation) override
    {
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
        // --- 1. Build ordered ring from the topology tracker ----------------
        //     Walk the doubly-linked list to get a consistent ordering.

        // Collect all live cell indices
        std::set<unsigned> liveIndices;
        std::map<unsigned, CellPtr> cellMap;
        for (auto it = rCellPopulation.Begin(); it != rCellPopulation.End(); ++it)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*it);
            liveIndices.insert(idx);
            cellMap[idx] = *it;
        }

        if (liveIndices.empty()) return;

        // Walk the ring starting from any live cell
        unsigned startIdx = *liveIndices.begin();
        std::vector<unsigned> ring;          // ordered cell indices
        ring.reserve(liveIndices.size());

        unsigned current = startIdx;
        do
        {
            ring.push_back(current);
            auto nbrs = mpRingTracker->GetNeighbors(current);
            unsigned next = nbrs.second;  // walk rightward
            if (next == UINT_MAX || liveIndices.find(next) == liveIndices.end())
            {
                break;  // broken ring — shouldn't happen
            }
            current = next;
        }
        while (current != startIdx);

        unsigned nPts = ring.size();
        if (nPts < 2) return;

        // Map ring-order index -> cell location index
        // Also build positions & per-node data
        std::vector<c_vector<double, 3>> positions(nPts);
        std::vector<double> cellTypeId(nPts, 0.0);
        std::vector<double> polarityTheta(nPts, 0.0);

        for (unsigned k = 0; k < nPts; k++)
        {
            unsigned idx = ring[k];
            c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(cellMap[idx]);
            positions[k][0] = pos[0];
            positions[k][1] = pos[1];
            positions[k][2] = (DIM == 3) ? pos[2] : 0.0;

            CellPtr pCell = cellMap[idx];
            if (pCell->GetCellProliferativeType()->template IsType<StemCellProliferativeType>())
                cellTypeId[k] = 0.0;
            else if (pCell->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
                cellTypeId[k] = 1.0;
            else
                cellTypeId[k] = 2.0;

            if (pCell->GetCellData()->HasItem("polarity_theta"))
                polarityTheta[k] = pCell->GetCellData()->GetItem("polarity_theta");
        }

        // --- 2. Write the .vtp file -----------------------------------------
        std::ostringstream fname;
        fname << "outline_" << timestep << ".vtp";
        std::string vtpFilename = fname.str();
        std::string fullPath = mOutputDirectory + vtpFilename;

        std::ofstream vtp(fullPath.c_str());
        if (!vtp.is_open())
        {
            EXCEPTION("RingOutlineWriter: cannot open " + fullPath);
        }

        unsigned nLines = nPts;  // one segment per consecutive pair, including wrap-around

        vtp << "<?xml version=\"1.0\"?>\n";
        vtp << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        vtp << "  <PolyData>\n";
        vtp << "    <Piece NumberOfPoints=\"" << nPts
            << "\" NumberOfVerts=\"0\""
            << " NumberOfLines=\"" << nLines
            << "\" NumberOfStrips=\"0\""
            << " NumberOfPolys=\"0\">\n";

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

        // --- Lines (connectivity + offsets) ----------------------------------
        vtp << "      <Lines>\n";

        // connectivity: pairs (k, k+1), wrapping last -> first
        vtp << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        vtp << "          ";
        for (unsigned k = 0; k < nPts; k++)
        {
            vtp << k << " " << ((k + 1) % nPts) << " ";
        }
        vtp << "\n        </DataArray>\n";

        // offsets: cumulative count of indices per line (each line has 2 points)
        vtp << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        vtp << "          ";
        for (unsigned k = 0; k < nLines; k++)
        {
            vtp << 2 * (k + 1) << " ";
        }
        vtp << "\n        </DataArray>\n";

        vtp << "      </Lines>\n";

        // --- PointData (per-node scalars) ------------------------------------
        vtp << "      <PointData Scalars=\"cell_type\">\n";

        vtp << "        <DataArray type=\"Float64\" Name=\"cell_type\" format=\"ascii\">\n";
        vtp << "          ";
        for (unsigned k = 0; k < nPts; k++)
        {
            vtp << cellTypeId[k] << " ";
        }
        vtp << "\n        </DataArray>\n";

        vtp << "        <DataArray type=\"Float64\" Name=\"polarity_theta\" format=\"ascii\">\n";
        vtp << "          ";
        for (unsigned k = 0; k < nPts; k++)
        {
            vtp << polarityTheta[k] << " ";
        }
        vtp << "\n        </DataArray>\n";

        vtp << "      </PointData>\n";

        // Close tags
        vtp << "    </Piece>\n";
        vtp << "  </PolyData>\n";
        vtp << "</VTKFile>\n";
        vtp.close();

        // --- 3. Append to PVD -----------------------------------------------
        if (mpPvdFile)
        {
            *mpPvdFile << "    <DataSet timestep=\"" << timestep
                       << "\" group=\"\" part=\"0\" file=\""
                       << vtpFilename << "\"/>\n";
        }
    }
};

#endif // RINGOUTLINEWRITER_HPP_

/*
 * CryptBuddingParams.hpp
 *
 * Parameter struct for crypt budding / organoid simulations.
 * Holds all tuneable knobs with model-specific defaults applied
 * by Finalise().
 *
 * Parameters can be loaded from an INI-style config file using LoadFromFile().
 * Format:
 *   # Comments start with #
 *   [Section]           # Sections are optional, for organization only
 *   parameterName = value
 *
 * Example:
 *   [Simulation]
 *   endTime = 168.0
 *   dt = 0.005
 *   enableRelaxation = true
 */
#ifndef CRYPTBUDDINGPARAMS_HPP_
#define CRYPTBUDDINGPARAMS_HPP_

#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>

struct CryptBuddingParams
{
    std::string modelType;
    double ecmStiffness;
    unsigned runNumber;
    unsigned randomSeed;

    bool enableLumenPressure;
    bool enableApicalConstriction;
    bool enableEcmGuidance;
    bool enableRelaxation;
    bool enableSloughing;
    bool enableDifferentialAdhesion;
    bool enableCellPolarity;
    bool enableEcmConfinement;
    bool enableContinuousPvd;

    double dt;
    double dtGrow;     // Phase 2 growth dt for vertex3d (separate from relaxation dt)
    double relaxationTime;
    double endTime;
    unsigned samplingMultiple;

    double organoidRadius2d;
    unsigned numCells2dNode;
    unsigned numCells2dVertex;
    double innerRadius2d;
    double outerRadius2d;
    double interactionCutoff2d;

    unsigned numCells3dNode;
    unsigned numCells3dVertex;
    double organoidRadius3d;
    double shellThickness3d;
    double interactionCutoff3d;
    double sphereRadius3dVertex;

    double bmStiffnessNode;
    double bmStiffnessVertex;
    double bmRadius2d;
    double bmRadius3d;
    double bmOffset3dVertex;
    double ecmDegradationRate;
    double ecmDiffusionCoeff;
    double ecmMaxRadius2d;
    double ecmMaxRadius3d;

    double lumenPressure;
    bool   lumenUseTargetVolume;      // Use incompressible fluid model
    double lumenVolumeGrowthRate;     // Volume growth rate (fraction per hour)
    double lumenBulkModulus;          // Fluid stiffness (pressure per fractional compression)

    double apicalConstrictionStrength;

    double springStiffness;
    double springCutoff;
    double apicalApicalAdhesion;
    double basalBasalAdhesion;
    double apicalBasalAdhesion;

    double nhMembraneSurface;
    double nhCellCellAdhesion;
    double nhBoundaryAdhesion;

    // Per-cell-type Nagai-Honda adhesion (symmetric 3×3 matrix)
    // Types: 0=stem, 1=transit, 2=differentiated
    double nhStemStemAdhesion;
    double nhStemTransitAdhesion;
    double nhStemDiffAdhesion;
    double nhTransitTransitAdhesion;
    double nhTransitDiffAdhesion;
    double nhDiffDiffAdhesion;
    // Per-type boundary adhesion
    double nhStemBoundaryAdhesion;
    double nhTransitBoundaryAdhesion;
    double nhDiffBoundaryAdhesion;

    double gammaApical;
    double gammaBasal;
    double gammaLateral;

    // Per-cell-type surface tension scaling (vertex models)
    // Stem cells are softer, Paneth/differentiated cells are stiffer
    double gammaStemScale;
    double gammaTransitScale;
    double gammaDiffScale;      // Paneth-like cells

    double quiescentFraction;
    double sloughRadiusFactor;

    double ecmDomainHalf;
    double ecmGridSpacing;
    double ecmBaseSpeed;
    std::string ecmGridType;     // "square" or "hex"

    double t1Threshold2d;
    double t2Threshold2d;

    // Curvature bending force (Drasdo 2000 - monolayer enforcement)
    bool enableCurvatureBending;
    double bendingStiffness;
    double lumenExclusionStrength;
    double minRadiusFraction;  // cells cannot go below this fraction of target radius

    // Spring neighbor strategy for node-based models
    // true  = topology-based (RingSpringForce / SurfaceSpringForce)
    // false = distance-threshold (GeneralisedLinearSpringForce)
    bool useTopologyBasedSprings;

    // Cell type proportions for uniform random distribution
    double stemFraction;
    double transitFraction;
    // remainder (1 - stemFraction - transitFraction) = differentiated/Paneth

    // Cell cycle duration parameters
    double stemCycleMin;        // Minimum total cycle for stem cells (hours)
    double stemCycleMax;        // Maximum total cycle for stem cells (hours)
    double taCycleRatio;        // TA cycle = ratio × stem cycle (0.5 or 1.0)
    // Paneth / differentiated cells never divide (handled by DifferentiatedCellProliferativeType)

    // Generational cascade (Meineke et al. 2001)
    bool enableGenerationalCascade;     // If true: Stem → TA → Differentiated cascade
    unsigned maxTransitGenerations;     // TA cells differentiate after this many divisions

    // Cell polarity parameters (monolayer enforcement via ya||a-style polarity)
    double polarityBendingStrength;     // Epithelial bending force strength
    double polarityAlignmentStrength;   // Tissue polarity alignment strength

    // Relative radius fractions — force radii derived as organoidRadius × fraction
    // These are the primary config knobs; absolute radii are auto-computed in Finalise()
    double bmRadiusFraction;            // BM target radius as fraction of organoidRadius
    double ecmMaxRadiusFraction;        // ECM degradation boundary as fraction of organoidRadius

    // Whether endTime/dt were explicitly set by user
    bool endTimeOverridden;
    bool dtOverridden;

    void SetDefaults()
    {
        modelType = "";
        ecmStiffness = 5.0;
        runNumber = 0;
        enableLumenPressure = true;
        enableApicalConstriction = true;
        enableEcmGuidance = false;
        enableRelaxation = true;
        enableSloughing = true;
        enableDifferentialAdhesion = true;
        enableEcmConfinement = true;
        enableContinuousPvd = false;
        endTimeOverridden = false;
        dtOverridden = false;
        endTime = 168.0;
        dt = 0.005;
        dtGrow = 0.002;   // Phase 2 growth dt for vertex3d (reduced from 0.006 for stability)

        relaxationTime = 10.0;

        organoidRadius2d    = 8.0;
        numCells2dNode      = 80;
        numCells2dVertex    = 40;
        innerRadius2d       = 6.0;
        outerRadius2d       = 8.0;
        interactionCutoff2d = 2.5;

        numCells3dNode       = 100;
        numCells3dVertex     = 200;
        organoidRadius3d     = 10.0;
        shellThickness3d     = 3.0;
        interactionCutoff3d  = 8.0;
        sphereRadius3dVertex = 10.0;  // derived from organoidRadius3d in Finalise()
        bmOffset3dVertex     = 1.0;

        lumenPressure          = 2.0;
        lumenUseTargetVolume   = false;   // Default to constant pressure mode
        lumenVolumeGrowthRate  = 0.05;    // 5% per hour
        lumenBulkModulus       = 10.0;    // Fluid stiffness

        apicalConstrictionStrength = 3.0;

        springStiffness       = 30.0;
        springCutoff          = 1.5;
        apicalApicalAdhesion  = 1.2;
        basalBasalAdhesion    = 1.0;
        apicalBasalAdhesion   = 0.5;

        nhMembraneSurface  = 10.0;
        nhCellCellAdhesion = 1.0;
        nhBoundaryAdhesion = 2.0;

        // Per-type adhesion defaults (all equal = uniform)
        nhStemStemAdhesion       = 1.0;
        nhStemTransitAdhesion    = 1.0;
        nhStemDiffAdhesion       = 1.0;
        nhTransitTransitAdhesion = 1.0;
        nhTransitDiffAdhesion    = 1.0;
        nhDiffDiffAdhesion       = 1.0;
        nhStemBoundaryAdhesion    = 2.0;
        nhTransitBoundaryAdhesion = 2.0;
        nhDiffBoundaryAdhesion    = 2.0;

        gammaApical  = 0.85;
        gammaBasal   = 0.85;
        gammaLateral = 0.7;

        // Per-cell-type scaling: stem cells softer, Paneth cells stiffer
        gammaStemScale    = 0.7;   // stem cells are more compliant
        gammaTransitScale = 1.0;   // baseline
        gammaDiffScale    = 1.3;   // Paneth/differentiated are stiffer

        quiescentFraction  = 0.7;
        sloughRadiusFactor = 5.0;

        ecmDomainHalf  = -1.0;  // sentinel: auto-derive in Finalise()
        ecmGridSpacing = 10.0;
        ecmBaseSpeed   = 0.3;
        ecmGridType    = "square";   // "square" or "hex"

        // Curvature bending force (Drasdo 2000 - monolayer enforcement)
        enableCurvatureBending    = true;   // Enable by default for node2d
        bendingStiffness          = 5.0;    // Bending rigidity
        lumenExclusionStrength    = 500.0;   // Strong repulsion from lumen interior
        minRadiusFraction         = 0.7;    // Cells stay outside 70% of target radius

        // Spring neighbor strategy (node-based models)
        useTopologyBasedSprings   = true;   // true = topology (ring/surface), false = distance threshold

        // Cell type proportions for uniform random distribution
        stemFraction    = 0.2;
        transitFraction = 0.5;

        // Cell cycle duration (uniform random total cycle)
        stemCycleMin  = 12.0;   // U(12, 14) h total cycle for stem cells
        stemCycleMax  = 14.0;
        taCycleRatio  = 1.0;    // TA cycle = ratio × stem cycle (set 0.5 for half)

        // Generational cascade (Meineke et al. 2001)
        enableGenerationalCascade = true;   // Enable Stem → TA → Differentiated cascade
        maxTransitGenerations = 3;          // TA cells differentiate after 3 divisions

        // Cell polarity (ya||a-style monolayer enforcement)
        enableCellPolarity          = true;   // Enable by default for node models
        polarityBendingStrength     = 0.3;    // Epithelial bending force
        polarityAlignmentStrength   = 0.1;    // Tissue polarity alignment

        // Relative radius fractions (absolute radii derived in Finalise)
        bmRadiusFraction       = 1.25;   // BM just outside cell ring
        ecmMaxRadiusFraction   = 4.0;    // ECM boundary far from organoid
    }

    void Finalise()
    {
        randomSeed = static_cast<unsigned>(ecmStiffness * 10000) + runNumber * 137;

        // Derive organoidRadius2d from numCells so spacing = rest length (1.0)
        // spacing = 2πR / N = 1.0  →  R = N / (2π)
        organoidRadius2d = numCells2dNode / (2.0 * M_PI);

        // Derive vertex-2D annulus radii so each wedge cell has area ≈ 1.0
        // Mid-radius matches node spacing: R_mid = N / (2π)
        // Ring width w = 2πR_mid / N = 1.0 (square aspect ratio)
        {
            double R_mid = numCells2dVertex / (2.0 * M_PI);
            double w = 2.0 * M_PI * R_mid / numCells2dVertex;  // = 1.0
            innerRadius2d = R_mid - 0.5 * w;
            outerRadius2d = R_mid + 0.5 * w;
        }

        // Derive organoidRadius3d from numCells so Voronoi cell area ≈ 1.0²
        // 4πR² / N = 1.0  →  R = √(N / 4π)
        organoidRadius3d = std::sqrt(numCells3dNode / (4.0 * M_PI));

        bmStiffnessNode   = ecmStiffness;
        bmStiffnessVertex = ecmStiffness * 0.5;
        ecmDegradationRate = 0.02;
        ecmDiffusionCoeff  = 0.1;    // ECM density smoothing coefficient

        // Derive all absolute radii from organoidRadius × fraction
        bmRadius2d         = organoidRadius2d * bmRadiusFraction;
        bmRadius3d         = organoidRadius3d * bmRadiusFraction;
        ecmMaxRadius2d     = organoidRadius2d * ecmMaxRadiusFraction;
        ecmMaxRadius3d     = organoidRadius3d * ecmMaxRadiusFraction;
        sphereRadius3dVertex = organoidRadius3d;

        // Auto-derive ECM domain from organoid size if not explicitly set
        if (ecmDomainHalf < 0.0)
        {
            double maxR = std::max(ecmMaxRadius2d, ecmMaxRadius3d);
            ecmDomainHalf = std::max(maxR, organoidRadius3d * 3.0) * 1.5;
            std::cout << "  ECM domain auto-derived: ecmDomainHalf = " << ecmDomainHalf
                      << " (organoidR=" << organoidRadius3d
                      << ", ecmMaxR=" << ecmMaxRadius3d << ")" << std::endl;
        }

        if (modelType == "node2d")
        {
            std::cout << "Node2D: " << numCells2dNode << " cells, R="
                      << organoidRadius2d << " (spacing="
                      << 2.0 * M_PI * organoidRadius2d / numCells2dNode << ")"
                      << "  BM=" << bmRadius2d << "  ECMmax=" << ecmMaxRadius2d << std::endl;
        }
        else if (modelType == "vertex2d")
        {
            std::cout << "Vertex2D: " << numCells2dVertex << " cells, Rin="
                      << innerRadius2d << " Rout=" << outerRadius2d
                      << " (area/cell="
                      << 0.5 * (2.0 * M_PI / numCells2dVertex)
                         * (outerRadius2d * outerRadius2d - innerRadius2d * innerRadius2d)
                      << ")"
                      << "  BM=" << bmRadius2d << "  ECMmax=" << ecmMaxRadius2d << std::endl;
        }
        else if (modelType == "node3d")
        {
            std::cout << "Node3D: " << numCells3dNode << " cells, R="
                      << organoidRadius3d << " (area/cell="
                      << 4.0 * M_PI * organoidRadius3d * organoidRadius3d / numCells3dNode
                      << ")"
                      << "  BM=" << bmRadius3d << "  ECMmax=" << ecmMaxRadius3d << std::endl;
        }
        else if (modelType == "vertex3d")
        {
            std::cout << "Vertex3D: " << numCells3dVertex << " cells, R="
                      << organoidRadius3d << " sphere=" << sphereRadius3dVertex
                      << "  BM=" << bmRadius3d << "  ECMmax=" << ecmMaxRadius3d << std::endl;
        }

        t1Threshold2d = (ecmStiffness < 2.0) ? 0.2 : 0.15;
        t2Threshold2d = 0.05;

        // Model-specific defaults for dt/endTime if not overridden
        if (modelType == "node2d")
        {
            if (!dtOverridden) dt = 0.005;
        }
        else if (modelType == "vertex2d")
        {
            if (!dtOverridden)
                dt = (ecmStiffness < 1.0) ? 0.0002
                   : (ecmStiffness < 5.0) ? 0.0005
                   :                        0.0005;
            if (!endTimeOverridden) endTime = 168.0;
        }
        else if (modelType == "node3d")
        {
            if (!dtOverridden) dt = 0.01;
            if (!endTimeOverridden) endTime = 168.0;
        }
        else if (modelType == "vertex3d")
        {
            if (!dtOverridden) dt = 0.0001;
            if (!endTimeOverridden) endTime = 100.0;
        }

        // Compute samplingMultiple so every simulation outputs exactly 50 frames
        const unsigned totalSteps = static_cast<unsigned>(std::round(endTime / dt));
        samplingMultiple = std::max(1u, totalSteps / 50);
    }

    /**
     * Load parameters from an INI-style config file.
     * File format:
     *   # Comments start with # or ;
     *   [Section]           # Sections are optional, ignored
     *   parameterName = value
     *
     * Only parameters present in the file are overwritten; others keep defaults.
     * Returns true if file was loaded successfully.
     */
    bool LoadFromFile(const std::string& filePath)
    {
        std::ifstream file(filePath);
        if (!file.is_open())
        {
            std::cerr << "Warning: Could not open config file: " << filePath << std::endl;
            return false;
        }

        std::map<std::string, std::string> configMap;
        std::string line;

        while (std::getline(file, line))
        {
            // Trim leading whitespace
            size_t start = line.find_first_not_of(" \t");
            if (start == std::string::npos) continue;
            line = line.substr(start);

            // Skip empty lines, comments, and section headers
            if (line.empty() || line[0] == '#' || line[0] == ';' || line[0] == '[')
                continue;

            // Parse key = value
            size_t eqPos = line.find('=');
            if (eqPos == std::string::npos) continue;

            std::string key = line.substr(0, eqPos);
            std::string value = line.substr(eqPos + 1);

            // Trim whitespace from key and value
            auto trim = [](std::string& s) {
                size_t start = s.find_first_not_of(" \t");
                size_t end = s.find_last_not_of(" \t");
                s = (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
            };
            trim(key);
            trim(value);

            // Remove inline comments from value
            size_t commentPos = value.find('#');
            if (commentPos != std::string::npos)
                value = value.substr(0, commentPos);
            commentPos = value.find(';');
            if (commentPos != std::string::npos)
                value = value.substr(0, commentPos);
            trim(value);

            if (!key.empty() && !value.empty())
            {
                configMap[key] = value;
            }
        }

        file.close();

        // Helper lambdas to parse values
        auto getBool = [&](const std::string& key, bool& var) {
            if (configMap.count(key)) {
                std::string v = configMap[key];
                std::transform(v.begin(), v.end(), v.begin(), ::tolower);
                var = (v == "true" || v == "1" || v == "yes" || v == "on");
            }
        };
        auto getDouble = [&](const std::string& key, double& var) {
            if (configMap.count(key)) {
                try { var = std::stod(configMap[key]); }
                catch (...) { std::cerr << "Warning: Invalid double for " << key << std::endl; }
            }
        };
        auto getUnsigned = [&](const std::string& key, unsigned& var) {
            if (configMap.count(key)) {
                try { var = static_cast<unsigned>(std::stoul(configMap[key])); }
                catch (...) { std::cerr << "Warning: Invalid unsigned for " << key << std::endl; }
            }
        };
        auto getString = [&](const std::string& key, std::string& var) {
            if (configMap.count(key)) var = configMap[key];
        };

        // Apply all parameters from config
        getString("modelType", modelType);
        getDouble("ecmStiffness", ecmStiffness);
        getUnsigned("runNumber", runNumber);
        getUnsigned("randomSeed", randomSeed);

        getBool("enableLumenPressure", enableLumenPressure);
        getBool("enableApicalConstriction", enableApicalConstriction);
        getBool("enableEcmGuidance", enableEcmGuidance);
        getBool("enableRelaxation", enableRelaxation);
        getBool("enableSloughing", enableSloughing);
        getBool("enableDifferentialAdhesion", enableDifferentialAdhesion);
        getBool("enableCurvatureBending", enableCurvatureBending);
        getBool("enableCellPolarity", enableCellPolarity);
        getBool("enableEcmConfinement", enableEcmConfinement);
        getBool("useTopologyBasedSprings", useTopologyBasedSprings);
        getBool("enableContinuousPvd", enableContinuousPvd);

        if (configMap.count("dt")) { getDouble("dt", dt); dtOverridden = true; }
        if (configMap.count("endTime")) { getDouble("endTime", endTime); endTimeOverridden = true; }
        getDouble("dtGrow", dtGrow);
        getDouble("relaxationTime", relaxationTime);

        getDouble("organoidRadius2d", organoidRadius2d);
        getUnsigned("numCells2dNode", numCells2dNode);
        getUnsigned("numCells2dVertex", numCells2dVertex);
        getDouble("innerRadius2d", innerRadius2d);
        getDouble("outerRadius2d", outerRadius2d);
        getDouble("interactionCutoff2d", interactionCutoff2d);

        getUnsigned("numCells3dNode", numCells3dNode);
        getUnsigned("numCells3dVertex", numCells3dVertex);
        getDouble("organoidRadius3d", organoidRadius3d);
        getDouble("shellThickness3d", shellThickness3d);
        getDouble("interactionCutoff3d", interactionCutoff3d);
        getDouble("sphereRadius3dVertex", sphereRadius3dVertex);
        getDouble("bmOffset3dVertex", bmOffset3dVertex);

        getDouble("bmStiffnessNode", bmStiffnessNode);
        getDouble("bmStiffnessVertex", bmStiffnessVertex);
        getDouble("ecmDegradationRate", ecmDegradationRate);
        getDouble("ecmDiffusionCoeff", ecmDiffusionCoeff);

        // Radius fractions (absolute radii derived in Finalise)
        getDouble("bmRadiusFraction", bmRadiusFraction);
        getDouble("ecmMaxRadiusFraction", ecmMaxRadiusFraction);

        getDouble("lumenPressure", lumenPressure);
        getBool("lumenUseTargetVolume", lumenUseTargetVolume);
        getDouble("lumenVolumeGrowthRate", lumenVolumeGrowthRate);
        getDouble("lumenBulkModulus", lumenBulkModulus);

        getDouble("apicalConstrictionStrength", apicalConstrictionStrength);

        getDouble("springStiffness", springStiffness);
        getDouble("springCutoff", springCutoff);
        getDouble("apicalApicalAdhesion", apicalApicalAdhesion);
        getDouble("basalBasalAdhesion", basalBasalAdhesion);
        getDouble("apicalBasalAdhesion", apicalBasalAdhesion);

        getDouble("nhMembraneSurface", nhMembraneSurface);
        getDouble("nhCellCellAdhesion", nhCellCellAdhesion);
        getDouble("nhBoundaryAdhesion", nhBoundaryAdhesion);

        // Per-type adhesion (override defaults if present)
        getDouble("nhStemStemAdhesion", nhStemStemAdhesion);
        getDouble("nhStemTransitAdhesion", nhStemTransitAdhesion);
        getDouble("nhStemDiffAdhesion", nhStemDiffAdhesion);
        getDouble("nhTransitTransitAdhesion", nhTransitTransitAdhesion);
        getDouble("nhTransitDiffAdhesion", nhTransitDiffAdhesion);
        getDouble("nhDiffDiffAdhesion", nhDiffDiffAdhesion);
        getDouble("nhStemBoundaryAdhesion", nhStemBoundaryAdhesion);
        getDouble("nhTransitBoundaryAdhesion", nhTransitBoundaryAdhesion);
        getDouble("nhDiffBoundaryAdhesion", nhDiffBoundaryAdhesion);

        getDouble("gammaApical", gammaApical);
        getDouble("gammaBasal", gammaBasal);
        getDouble("gammaLateral", gammaLateral);
        getDouble("gammaStemScale", gammaStemScale);
        getDouble("gammaTransitScale", gammaTransitScale);
        getDouble("gammaDiffScale", gammaDiffScale);

        getDouble("quiescentFraction", quiescentFraction);
        getDouble("sloughRadiusFactor", sloughRadiusFactor);

        getDouble("ecmDomainHalf", ecmDomainHalf);
        getDouble("ecmGridSpacing", ecmGridSpacing);
        getDouble("ecmBaseSpeed", ecmBaseSpeed);
        getString("ecmGridType", ecmGridType);

        getDouble("t1Threshold2d", t1Threshold2d);
        getDouble("t2Threshold2d", t2Threshold2d);

        getDouble("bendingStiffness", bendingStiffness);
        getDouble("lumenExclusionStrength", lumenExclusionStrength);
        getDouble("minRadiusFraction", minRadiusFraction);

        getDouble("stemFraction", stemFraction);
        getDouble("transitFraction", transitFraction);

        getDouble("stemCycleMin", stemCycleMin);
        getDouble("stemCycleMax", stemCycleMax);
        getDouble("taCycleRatio", taCycleRatio);

        getBool("enableGenerationalCascade", enableGenerationalCascade);
        getUnsigned("maxTransitGenerations", maxTransitGenerations);

        getDouble("polarityBendingStrength", polarityBendingStrength);
        getDouble("polarityAlignmentStrength", polarityAlignmentStrength);

        std::cout << "Loaded " << configMap.size() << " parameters from: " << filePath << std::endl;
        return true;
    }

    /**
     * Save current parameters to an INI-style config file.
     * Useful for generating a template or recording simulation settings.
     */
    bool SaveToFile(const std::string& filePath) const
    {
        std::ofstream file(filePath);
        if (!file.is_open())
        {
            std::cerr << "Error: Could not write config file: " << filePath << std::endl;
            return false;
        }

        file << "# CryptBudding Simulation Parameters\n";
        file << "# Generated config file - edit values as needed\n\n";

        file << "[Simulation]\n";
        file << "modelType = " << modelType << "\n";
        file << "endTime = " << endTime << "          # Total simulation time (hours)\n";
        file << "dt = " << dt << "              # Time step\n";
        file << "dtGrow = " << dtGrow << "          # Growth phase dt (vertex3d)\n";
        file << "relaxationTime = " << relaxationTime << "    # Relaxation phase duration\n";
        file << "runNumber = " << runNumber << "\n";
        file << "randomSeed = " << randomSeed << "\n\n";

        file << "[Features]\n";
        file << "enableLumenPressure = " << (enableLumenPressure ? "true" : "false") << "\n";
        file << "enableApicalConstriction = " << (enableApicalConstriction ? "true" : "false") << "\n";
        file << "enableEcmGuidance = " << (enableEcmGuidance ? "true" : "false") << "\n";
        file << "enableRelaxation = " << (enableRelaxation ? "true" : "false") << "\n";
        file << "enableSloughing = " << (enableSloughing ? "true" : "false") << "\n";
        file << "enableDifferentialAdhesion = " << (enableDifferentialAdhesion ? "true" : "false") << "\n";
        file << "enableCurvatureBending = " << (enableCurvatureBending ? "true" : "false") << "\n";
        file << "enableCellPolarity = " << (enableCellPolarity ? "true" : "false") << "\n";
        file << "useTopologyBasedSprings = " << (useTopologyBasedSprings ? "true" : "false") << "  # true=topology (ring/surface), false=distance threshold\n";
        file << "enableContinuousPvd = " << (enableContinuousPvd ? "true" : "false") << "  # Keep .pvd files valid during simulation\n\n";

        file << "[ECM]\n";
        file << "ecmStiffness = " << ecmStiffness << "       # Main ECM/BM stiffness parameter\n";
        file << "bmStiffnessNode = " << bmStiffnessNode << "\n";
        file << "bmStiffnessVertex = " << bmStiffnessVertex << "\n";
        file << "bmOffset3dVertex = " << bmOffset3dVertex << "\n";
        file << "ecmDegradationRate = " << ecmDegradationRate << "\n";
        file << "ecmDiffusionCoeff = " << ecmDiffusionCoeff << "  # Density smoothing coefficient\n\n";

        file << "# Radius fractions — absolute radii = organoidRadius × fraction\n";
        file << "bmRadiusFraction = " << bmRadiusFraction << "     # BM target radius (derived: 2D=" << bmRadius2d << " 3D=" << bmRadius3d << ")\n";

        file << "ecmMaxRadiusFraction = " << ecmMaxRadiusFraction << "  # ECM boundary (derived: 2D=" << ecmMaxRadius2d << " 3D=" << ecmMaxRadius3d << ")\n\n";

        file << "[Geometry2D]\n";
        file << "organoidRadius2d = " << organoidRadius2d << "\n";
        file << "numCells2dNode = " << numCells2dNode << "\n";
        file << "numCells2dVertex = " << numCells2dVertex << "\n";
        file << "innerRadius2d = " << innerRadius2d << "\n";
        file << "outerRadius2d = " << outerRadius2d << "\n";
        file << "interactionCutoff2d = " << interactionCutoff2d << "\n";
        file << "t1Threshold2d = " << t1Threshold2d << "\n";
        file << "t2Threshold2d = " << t2Threshold2d << "\n\n";

        file << "[Geometry3D]\n";
        file << "numCells3dNode = " << numCells3dNode << "\n";
        file << "numCells3dVertex = " << numCells3dVertex << "\n";
        file << "organoidRadius3d = " << organoidRadius3d << "\n";
        file << "shellThickness3d = " << shellThickness3d << "\n";
        file << "interactionCutoff3d = " << interactionCutoff3d << "\n";
        file << "sphereRadius3dVertex = " << sphereRadius3dVertex << "\n\n";

        file << "[Forces]\n";
        file << "# Spring/adhesion forces\n";
        file << "springStiffness = " << springStiffness << "\n";
        file << "springCutoff = " << springCutoff << "\n";
        file << "apicalApicalAdhesion = " << apicalApicalAdhesion << "\n";
        file << "basalBasalAdhesion = " << basalBasalAdhesion << "\n";
        file << "apicalBasalAdhesion = " << apicalBasalAdhesion << "\n\n";

        file << "# Lumen pressure\n";
        file << "lumenPressure = " << lumenPressure << "\n";
        file << "lumenUseTargetVolume = " << (lumenUseTargetVolume ? "true" : "false") << "   # Incompressible fluid model\n";
        file << "lumenVolumeGrowthRate = " << lumenVolumeGrowthRate << "   # fraction per hour\n";
        file << "lumenBulkModulus = " << lumenBulkModulus << "   # fluid stiffness\n\n";

        file << "# Apical constriction\n";
        file << "apicalConstrictionStrength = " << apicalConstrictionStrength << "\n\n";

        file << "# Curvature bending (Drasdo 2000 - monolayer enforcement)\n";
        file << "bendingStiffness = " << bendingStiffness << "\n";
        file << "lumenExclusionStrength = " << lumenExclusionStrength << "\n";
        file << "minRadiusFraction = " << minRadiusFraction << "   # Cells stay outside this fraction of target radius\n\n";

        file << "# Cell polarity (ya||a-style monolayer enforcement)\n";
        file << "polarityBendingStrength = " << polarityBendingStrength << "    # Epithelial bending force\n";
        file << "polarityAlignmentStrength = " << polarityAlignmentStrength << "  # Tissue polarity alignment\n\n";

        file << "[VertexModel]\n";
        file << "# Nagai-Honda parameters\n";
        file << "nhMembraneSurface = " << nhMembraneSurface << "\n";
        file << "nhCellCellAdhesion = " << nhCellCellAdhesion << "   # Uniform fallback\n";
        file << "nhBoundaryAdhesion = " << nhBoundaryAdhesion << "   # Uniform fallback\n";
        file << "# Per-type adhesion (symmetric)\n";
        file << "nhStemStemAdhesion = " << nhStemStemAdhesion << "\n";
        file << "nhStemTransitAdhesion = " << nhStemTransitAdhesion << "\n";
        file << "nhStemDiffAdhesion = " << nhStemDiffAdhesion << "\n";
        file << "nhTransitTransitAdhesion = " << nhTransitTransitAdhesion << "\n";
        file << "nhTransitDiffAdhesion = " << nhTransitDiffAdhesion << "\n";
        file << "nhDiffDiffAdhesion = " << nhDiffDiffAdhesion << "\n";
        file << "nhStemBoundaryAdhesion = " << nhStemBoundaryAdhesion << "\n";
        file << "nhTransitBoundaryAdhesion = " << nhTransitBoundaryAdhesion << "\n";
        file << "nhDiffBoundaryAdhesion = " << nhDiffBoundaryAdhesion << "\n\n";

        file << "# Surface tension\n";
        file << "gammaApical = " << gammaApical << "\n";
        file << "gammaBasal = " << gammaBasal << "\n";
        file << "gammaLateral = " << gammaLateral << "\n";
        file << "gammaStemScale = " << gammaStemScale << "     # Stem cells are softer\n";
        file << "gammaTransitScale = " << gammaTransitScale << "\n";
        file << "gammaDiffScale = " << gammaDiffScale << "      # Paneth/differentiated are stiffer\n\n";

        file << "[CellTypes]\n";
        file << "stemFraction = " << stemFraction << "        # Fraction of stem cells\n";
        file << "transitFraction = " << transitFraction << "     # Fraction of transit-amplifying cells\n";
        file << "# Remainder = differentiated/Paneth cells\n\n";

        file << "[CellCycle]\n";
        file << "quiescentFraction = " << quiescentFraction << "  # Contact inhibition threshold\n";
        file << "sloughRadiusFactor = " << sloughRadiusFactor << " # Sloughing boundary = radius * factor\n";
        file << "stemCycleMin = " << stemCycleMin << "        # Stem total cycle min (hours)\n";
        file << "stemCycleMax = " << stemCycleMax << "        # Stem total cycle max (hours)\n";
        file << "taCycleRatio = " << taCycleRatio << "        # TA cycle = ratio * stem (0.5 or 1)\n";
        file << "enableGenerationalCascade = " << (enableGenerationalCascade ? "true" : "false") << "  # Stem→TA→Diff cascade (Meineke 2001)\n";
        file << "maxTransitGenerations = " << maxTransitGenerations << "       # TA divisions before differentiation\n\n";

        file << "[ECMGuidance]\n";
        file << "ecmDomainHalf = " << ecmDomainHalf << "\n";
        file << "ecmGridSpacing = " << ecmGridSpacing << "\n";
        file << "ecmBaseSpeed = " << ecmBaseSpeed << "\n";
        file << "ecmGridType = " << ecmGridType << "\n";

        file.close();
        std::cout << "Saved parameters to: " << filePath << std::endl;
        return true;
    }
};

#endif // CRYPTBUDDINGPARAMS_HPP_

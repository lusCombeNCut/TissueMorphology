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
    double ecmMaxRadius2d;
    double ecmMaxRadius3d;

    double lumenPressure;
    double lumenEqRadius2d;
    double lumenEqRadius3d;

    double apicalConstrictionStrength;

    double springStiffness;
    double springCutoff;
    double apicalApicalAdhesion;
    double basalBasalAdhesion;
    double apicalBasalAdhesion;

    double nhMembraneSurface;
    double nhCellCellAdhesion;
    double nhBoundaryAdhesion;

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

    double t1Threshold2d;
    double t2Threshold2d;

    // Curvature bending force (Drasdo 2000 - monolayer enforcement)
    bool enableCurvatureBending;
    double bendingStiffness;
    double lumenExclusionStrength;
    double minRadiusFraction;  // cells cannot go below this fraction of target radius

    // Cell type proportions for uniform random distribution
    double stemFraction;
    double transitFraction;
    // remainder (1 - stemFraction - transitFraction) = differentiated/Paneth

    // Cell polarity parameters (monolayer enforcement via ya||a-style polarity)
    double polarityBendingStrength;     // Epithelial bending force strength
    double polarityAlignmentStrength;   // Tissue polarity alignment strength

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

        lumenPressure   = 2.0;
        lumenEqRadius2d = 0.0;     // derived in Finalise()
        lumenEqRadius3d = 0.0;     // derived in Finalise()

        apicalConstrictionStrength = 3.0;

        springStiffness       = 30.0;
        springCutoff          = 1.5;
        apicalApicalAdhesion  = 1.2;
        basalBasalAdhesion    = 1.0;
        apicalBasalAdhesion   = 0.5;

        nhMembraneSurface  = 10.0;
        nhCellCellAdhesion = 1.0;
        nhBoundaryAdhesion = 2.0;

        gammaApical  = 0.85;
        gammaBasal   = 0.85;
        gammaLateral = 0.7;

        // Per-cell-type scaling: stem cells softer, Paneth cells stiffer
        gammaStemScale    = 0.7;   // stem cells are more compliant
        gammaTransitScale = 1.0;   // baseline
        gammaDiffScale    = 1.3;   // Paneth/differentiated are stiffer

        quiescentFraction  = 0.7;
        sloughRadiusFactor = 5.0;

        ecmDomainHalf  = 80.0;
        ecmGridSpacing = 10.0;
        ecmBaseSpeed   = 0.3;

        // Curvature bending force (Drasdo 2000 - monolayer enforcement)
        enableCurvatureBending    = true;   // Enable by default for node2d
        bendingStiffness          = 5.0;    // Bending rigidity
        lumenExclusionStrength    = 500.0;   // Strong repulsion from lumen interior
        minRadiusFraction         = 0.7;    // Cells stay outside 70% of target radius

        // Cell type proportions for uniform random distribution
        stemFraction    = 0.2;
        transitFraction = 0.5;

        // Cell polarity (ya||a-style monolayer enforcement)
        enableCellPolarity          = true;   // Enable by default for node models
        polarityBendingStrength     = 0.3;    // Epithelial bending force
        polarityAlignmentStrength   = 0.1;    // Tissue polarity alignment
    }

    void Finalise()
    {
        randomSeed = static_cast<unsigned>(ecmStiffness * 10000) + runNumber * 137;

        bmStiffnessNode   = ecmStiffness;
        bmStiffnessVertex = ecmStiffness * 0.5;
        bmRadius2d         = organoidRadius2d + 2.0;
        bmRadius3d         = organoidRadius3d + 2.0;
        ecmDegradationRate = 0.02;
        ecmMaxRadius2d     = organoidRadius2d * 4.0;
        ecmMaxRadius3d     = organoidRadius3d * 4.0;

        // Standardised derived radii — single source of truth from organoidRadius
        sphereRadius3dVertex = organoidRadius3d;
        lumenEqRadius2d      = organoidRadius2d + 1.0;
        lumenEqRadius3d      = organoidRadius3d + 1.0;

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
        getDouble("bmRadius2d", bmRadius2d);
        getDouble("bmRadius3d", bmRadius3d);
        getDouble("ecmDegradationRate", ecmDegradationRate);
        getDouble("ecmMaxRadius2d", ecmMaxRadius2d);
        getDouble("ecmMaxRadius3d", ecmMaxRadius3d);

        getDouble("lumenPressure", lumenPressure);
        getDouble("lumenEqRadius2d", lumenEqRadius2d);
        getDouble("lumenEqRadius3d", lumenEqRadius3d);

        getDouble("apicalConstrictionStrength", apicalConstrictionStrength);

        getDouble("springStiffness", springStiffness);
        getDouble("springCutoff", springCutoff);
        getDouble("apicalApicalAdhesion", apicalApicalAdhesion);
        getDouble("basalBasalAdhesion", basalBasalAdhesion);
        getDouble("apicalBasalAdhesion", apicalBasalAdhesion);

        getDouble("nhMembraneSurface", nhMembraneSurface);
        getDouble("nhCellCellAdhesion", nhCellCellAdhesion);
        getDouble("nhBoundaryAdhesion", nhBoundaryAdhesion);

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

        getDouble("t1Threshold2d", t1Threshold2d);
        getDouble("t2Threshold2d", t2Threshold2d);

        getDouble("bendingStiffness", bendingStiffness);
        getDouble("lumenExclusionStrength", lumenExclusionStrength);
        getDouble("minRadiusFraction", minRadiusFraction);

        getDouble("stemFraction", stemFraction);
        getDouble("transitFraction", transitFraction);

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
        file << "enableCellPolarity = " << (enableCellPolarity ? "true" : "false") << "\n\n";

        file << "[ECM]\n";
        file << "ecmStiffness = " << ecmStiffness << "       # Main ECM/BM stiffness parameter\n";
        file << "bmStiffnessNode = " << bmStiffnessNode << "\n";
        file << "bmStiffnessVertex = " << bmStiffnessVertex << "\n";
        file << "bmRadius2d = " << bmRadius2d << "\n";
        file << "bmRadius3d = " << bmRadius3d << "\n";
        file << "bmOffset3dVertex = " << bmOffset3dVertex << "\n";
        file << "ecmDegradationRate = " << ecmDegradationRate << "\n";
        file << "ecmMaxRadius2d = " << ecmMaxRadius2d << "\n";
        file << "ecmMaxRadius3d = " << ecmMaxRadius3d << "\n\n";

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
        file << "lumenEqRadius2d = " << lumenEqRadius2d << "\n";
        file << "lumenEqRadius3d = " << lumenEqRadius3d << "\n\n";

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
        file << "nhCellCellAdhesion = " << nhCellCellAdhesion << "\n";
        file << "nhBoundaryAdhesion = " << nhBoundaryAdhesion << "\n\n";

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
        file << "sloughRadiusFactor = " << sloughRadiusFactor << " # Sloughing boundary = radius * factor\n\n";

        file << "[ECMGuidance]\n";
        file << "ecmDomainHalf = " << ecmDomainHalf << "\n";
        file << "ecmGridSpacing = " << ecmGridSpacing << "\n";
        file << "ecmBaseSpeed = " << ecmBaseSpeed << "\n";

        file.close();
        std::cout << "Saved parameters to: " << filePath << std::endl;
        return true;
    }
};

#endif // CRYPTBUDDINGPARAMS_HPP_

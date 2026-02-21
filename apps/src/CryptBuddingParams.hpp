/*
 * CryptBuddingParams.hpp
 *
 * Parameter struct for crypt budding / organoid simulations.
 * Holds all tuneable knobs with model-specific defaults applied
 * by Finalise().
 */
#ifndef CRYPTBUDDINGPARAMS_HPP_
#define CRYPTBUDDINGPARAMS_HPP_

#include <string>
#include <cmath>
#include <algorithm>

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

    // Cell type proportions for uniform random distribution
    double stemFraction;
    double transitFraction;
    // remainder (1 - stemFraction - transitFraction) = differentiated/Paneth

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

        // Cell type proportions for uniform random distribution
        stemFraction    = 0.2;
        transitFraction = 0.5;
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
};

#endif // CRYPTBUDDINGPARAMS_HPP_

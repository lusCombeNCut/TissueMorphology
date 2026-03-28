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
#include <iomanip>
#include <map>
#include <set>
#include <functional>

struct CryptBuddingParams
{
    std::string modelType;
    double ecmStiffness;
    unsigned runNumber;
    unsigned randomSeed;
    int verbosity;  // 0=quiet, 1=normal, 2=verbose

    // ── Unit conversion system ───────────────────────────────────
    // When usePhysicalUnits=true, parameters are converted from physical units
    // (kPa, Pa, N/m, etc.) to simulation units.
    //
    // CHASTE solves: v_sim = F_sim  (overdamped, η_sim = 1)
    // Physical mapping: v_phys = v_sim * (L0/T0),  F_phys = F_sim * η_phys * (L0/T0)
    // Therefore:  k_phys [N/m] = k_sim * η_phys / T0
    //             P_phys [Pa]  = P_sim * η_phys / (T0 * L0)
    //
    // η_phys is estimated from empirical traction force data:
    //   η_phys = F_ref / v_ref
    // where F_ref = typical traction force (nN) and v_ref = typical migration speed (µm/h).
    // For epithelial cells: ~10 nN at ~10 µm/h → η ≈ 3.6 N·s/m
    // This is ~40,000× higher than Stokes drag in water, reflecting focal adhesion coupling.
    //
    // Default values maintain backward compatibility (usePhysicalUnits=false).
    bool usePhysicalUnits;                      // Enable/disable physical unit conversion
    double cellDiameterMicrometers;             // Cell diameter in µm (default 10 µm)
    double referenceForcenN;                    // Typical traction force in nN (default 10 nN)
    double referenceVelocityUmPerH;             // Typical migration speed in µm/h (default 10 µm/h)
    double ecmElasticModulusPa;                 // ECM elastic modulus in Pa (optional; overrides ecmStiffness if > 0)

    bool enableLumenPressure;
    bool enableApicalConstriction;
    bool enableEcmGuidance;
    bool enableRelaxation;
    bool enableSloughing;
    bool enableDifferentialAdhesion;
    bool enableCellPolarity;
    bool enableEcmConfinement;
    bool enableLumenHole;
    bool enableContinuousPvd;

    double dt;
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

    double apicalConstrictionStrength;

    double springStiffness;
    double springCutoff;
    double springStiffnessTAScale;    // TA cell spring stiffness = springStiffness × this
    double springStiffnessDiffScale;  // Differentiated cell spring stiffness = springStiffness × this
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
    double ecmSpringRestLength;      // Rest length for cell-ECM springs (0 = adhesion to ECM position)
    double ecmInteractionCutoff;     // Cutoff distance for cell-ECM spring interactions

    // Ghost node ECM (off-lattice particle-based ECM)
    bool enableGhostNodeECM;            // Use ghost nodes instead of grid-based ECM
    double ghostDamping;                // ECM node viscous damping coefficient
    double ghostRemovalThreshold;       // Density below which ghost nodes are removed
    double ghostFibreRemodelingRate;    // Rate of fibre reorientation by cell traction
    double ghostAnisotropyStrength;     // Fibre anisotropy factor [0,1]
    double ghostCellGhostRestLength;    // Cell-ghost spring rest length
    double ghostCellGhostCutoff;        // Cell-ghost interaction cutoff
    double ghostSpringRestLength;       // ECM-ECM spring rest length
    double ghostGridSpacing;            // Ghost node grid spacing (defaults to ghostSpringRestLength)
    unsigned ghostRemovalCheckInterval; // Steps between removal checks

    // Viscoelastic constitutive model (generalised Maxwell: E(t) = E0 + E1*exp(-t/tau))
    bool enableViscoelasticECM;             // Use viscoelastic ghost node ECM instead of elastic
    double ghostRelaxedStiffness;           // Long-time (relaxed) modulus E0
    double ghostRelaxationModulus;          // Transient modulus E1
    double ghostRelaxationTime;             // Relaxation time tau (hours)

    // Dynamic ghost node boundary expansion
    bool enableGhostBoundaryExpansion;      // Dynamically add ghost nodes when cells approach domain boundary

    double t1Threshold2d;
    double t2Threshold2d;

    // Curvature bending force (Drasdo 2000 - monolayer enforcement)
    bool enableCurvatureBending;
    double bendingStiffness;

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

    // Stochastic 4-type cell cycle (Montes-Olivas et al. 2023)
    bool enableStochasticFourType;       // If true: use stochastic SC/TA/PC/EC transitions
    double probStemToStem;               // P(SC daughter = SC)
    double probStemToPaneth;             // P(SC daughter = PC)
    // P(SC daughter = TA) = 1 - probStemToStem - probStemToPaneth
    double probTaToTaEarly;              // P(TA daughter = TA) for t < transitionTime
    double probTaToTaLate;               // P(TA daughter = TA) for t >= transitionTime
    double transitionTime;               // Time (hours) at which TA probabilities switch (day 5 = 120h)
    double panethFraction;               // Initial fraction of Paneth cells (4-type model only)

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
    bool t1ThresholdOverridden;
    bool t2ThresholdOverridden;

    // Keys explicitly loaded from config file (populated by ApplyConfigMap)
    std::set<std::string> mLoadedKeys;

    /**
     * Compute and apply unit conversion factors if usePhysicalUnits is enabled.
     *
     * Converts parameters from physical units to simulation units using:
     *   - L0 = cellDiameterMicrometers * 1e-6 (metres), T0 = 3600 s
     *   - η_phys = referenceForcenN * 1e-9 / (referenceVelocityUmPerH * 1e-6 / T0)
     *   - k_factor = η_phys / T0,  P_factor = η_phys / (T0 * L0)
     *
     * If usePhysicalUnits is false, no conversion is applied but a table is still printed
     * showing the sim values alongside their physical equivalents (sim × factor).
     *
     * Called from ApplyConfigMap() after all parameters are loaded from config files.
     * This ensures conversion is applied only to parameters that were explicitly set.
     */
    void ComputeConversionFactors()
    {
        const double L0 = cellDiameterMicrometers * 1e-6;  // cell diameter in metres
        const double T0 = 3600.0;                           // 1 hour in seconds

        // Effective tissue damping from traction force data:
        //   η_phys = F_ref / v_ref
        const double f_ref_N  = referenceForcenN * 1e-9;
        const double v_ref_ms = (referenceVelocityUmPerH * 1e-6) / T0;
        const double etaPhys  = f_ref_N / v_ref_ms;  // N·s/m

        // Conversion factors (physical = sim × factor):
        //   k_phys [N/m]  = k_sim  × η_phys / T0
        //   P_phys [Pa]   = P_sim  × η_phys / (T0 × L0)
        //   κ_phys [N·m]  = κ_sim  × L0²   (geometric scaling)
        const double kF = etaPhys / T0;        // stiffness: N/m per sim unit
        const double pF = etaPhys / (T0 * L0); // pressure:  Pa per sim unit
        const double bF = L0 * L0;             // bending:   N·m per sim unit

        if (usePhysicalUnits)
        {
            // If elastic modulus is specified, derive ecmStiffness in physical N/m first:
            //   k_phys [N/m] = E [Pa] × L0 [m]
            // Then the common /kF below converts it to sim units.
            if (ecmElasticModulusPa > 0.0)
                ecmStiffness = ecmElasticModulusPa * L0;

            // Convert physical → sim (divide by factor)
            ecmStiffness             /= kF;
            springStiffness          /= kF;
            bendingStiffness         /= bF;
            ghostRelaxedStiffness    /= kF;
            ghostRelaxationModulus   /= kF;
            gammaApical              /= kF;
            gammaBasal               /= kF;
            gammaLateral             /= kF;
            lumenPressure            /= pF;
            // apicalConstrictionStrength is NOT converted: it is a geometric
            // stiffness coefficient (F = k * area_difference) where area is in
            // [CD²], so the units are η/(CD·h) — not pressure.
            polarityBendingStrength  /= pF;
            nhMembraneSurface        /= kF;
            nhCellCellAdhesion       /= kF;
            nhBoundaryAdhesion       /= kF;
            nhStemStemAdhesion       /= kF;
            nhStemTransitAdhesion    /= kF;
            nhStemDiffAdhesion       /= kF;
            nhTransitTransitAdhesion /= kF;
            nhTransitDiffAdhesion    /= kF;
            nhDiffDiffAdhesion       /= kF;
            nhStemBoundaryAdhesion   /= kF;
            nhTransitBoundaryAdhesion /= kF;
            nhDiffBoundaryAdhesion   /= kF;
        }

    }

    /**
     * Print a table of key force/stiffness parameters with their physical equivalents.
     * When usePhysicalUnits=false: sim → physical (useful for debugging)
     * When usePhysicalUnits=true:  already-converted sim values with the physical
     *                              inputs they came from.
     *
     * Call from CryptBuddingApp.cpp after ParseArguments() so verbosity is known.
     * Suggested: call when verbosity >= 2 (via -v) or always when usePhysicalUnits=true.
     */
    void PrintUnitTable() const
    {
        const double L0 = cellDiameterMicrometers * 1e-6;
        const double T0 = 3600.0;
        const double f_ref_N  = referenceForcenN * 1e-9;
        const double v_ref_ms = (referenceVelocityUmPerH * 1e-6) / T0;
        const double etaPhys  = f_ref_N / v_ref_ms;
        const double kF = etaPhys / T0;
        const double pF = etaPhys / (T0 * L0);
        const double bF = L0 * L0;

        auto row = [&](const std::string& name,
                       double sim, const std::string& physUnit,
                       double phys)
        {
            std::cout << "  " << std::left << std::setw(34) << name
                      << std::right << std::setw(14) << std::scientific << std::setprecision(4) << sim
                      << "  (sim)  →  "
                      << std::right << std::setw(12) << std::scientific << std::setprecision(4) << phys
                      << "  " << physUnit << "\n";
        };

        std::cout << "\n── Parameter Units Table ───────────────────────────────────────────────────\n";
        std::cout << "  Mode: " << (usePhysicalUnits ? "physical input converted to sim" : "sim input → physical equivalents") << "\n";
        std::cout << "  Cell diameter: " << cellDiameterMicrometers << " um"
                  << "   η_phys: " << etaPhys << " N·s/m"
                  << "   L0: " << L0 << " m   T0: " << T0 << " s\n";
        std::cout << "  Stiffness factor: " << kF << " N/m per sim unit"
                  << "   Pressure factor: " << pF << " Pa per sim unit\n";
        std::cout << "  " << std::string(76, '-') << "\n";
        std::cout << "  " << std::left << std::setw(34) << "Parameter"
                  << std::right << std::setw(14) << "Sim value"
                  << "           "
                  << std::right << std::setw(12) << "Physical"
                  << "  Units\n";
        std::cout << "  " << std::string(76, '-') << "\n";

        if (ecmElasticModulusPa > 0.0)
            std::cout << "  " << std::left << std::setw(34) << "ecmElasticModulusPa (input)"
                      << "  " << std::scientific << std::setprecision(4) << ecmElasticModulusPa << "  Pa\n";

        row("ecmStiffness",              ecmStiffness,              "N/m",  ecmStiffness * kF);
        row("springStiffness",           springStiffness,           "N/m",  springStiffness * kF);
        row("bendingStiffness",          bendingStiffness,          "N·m",  bendingStiffness * bF);
        row("ghostRelaxedStiffness",     ghostRelaxedStiffness,     "N/m",  ghostRelaxedStiffness * kF);
        row("ghostRelaxationModulus",    ghostRelaxationModulus,    "N/m",  ghostRelaxationModulus * kF);
        row("gammaApical",               gammaApical,               "N/m",  gammaApical * kF);
        row("gammaBasal",                gammaBasal,                "N/m",  gammaBasal * kF);
        row("gammaLateral",              gammaLateral,              "N/m",  gammaLateral * kF);
        row("lumenPressure",             lumenPressure,             "Pa",   lumenPressure * pF);
        row("apicalConstrictionStrength",apicalConstrictionStrength,"(sim)", apicalConstrictionStrength); // geometric coeff, not converted
        row("polarityBendingStrength",   polarityBendingStrength,   "Pa",   polarityBendingStrength * pF);
        row("nhMembraneSurface",         nhMembraneSurface,         "N/m",  nhMembraneSurface * kF);
        row("nhCellCellAdhesion",        nhCellCellAdhesion,        "N/m",  nhCellCellAdhesion * kF);
        row("nhBoundaryAdhesion",        nhBoundaryAdhesion,        "N/m",  nhBoundaryAdhesion * kF);
        row("nhStemStemAdhesion",        nhStemStemAdhesion,        "N/m",  nhStemStemAdhesion * kF);
        row("nhStemTransitAdhesion",     nhStemTransitAdhesion,     "N/m",  nhStemTransitAdhesion * kF);
        row("nhStemDiffAdhesion",        nhStemDiffAdhesion,        "N/m",  nhStemDiffAdhesion * kF);
        row("nhTransitTransitAdhesion",  nhTransitTransitAdhesion,  "N/m",  nhTransitTransitAdhesion * kF);
        row("nhTransitDiffAdhesion",     nhTransitDiffAdhesion,     "N/m",  nhTransitDiffAdhesion * kF);
        row("nhDiffDiffAdhesion",        nhDiffDiffAdhesion,        "N/m",  nhDiffDiffAdhesion * kF);
        row("nhStemBoundaryAdhesion",    nhStemBoundaryAdhesion,    "N/m",  nhStemBoundaryAdhesion * kF);
        row("nhTransitBoundaryAdhesion", nhTransitBoundaryAdhesion, "N/m",  nhTransitBoundaryAdhesion * kF);
        row("nhDiffBoundaryAdhesion",    nhDiffBoundaryAdhesion,    "N/m",  nhDiffBoundaryAdhesion * kF);

        std::cout << "  " << std::string(76, '-') << "\n\n";
    }

    /**
     * Initialise internal/sentinel fields only.
     *
     * IMPORTANT: Physics parameter defaults (spring stiffness, adhesion values,
     * surface tensions, ECM parameters, etc.) are defined in
     * apps/config/default_params.json — NOT here. ParseArguments() auto-loads
     * that file before applying any user config. The physics values below are
     * fallbacks only, used when default_params.json cannot be found (e.g. in
     * unit tests that construct CryptBuddingParams directly).
     *
     * Sentinel values (negative = "auto-derive in Finalise()") must remain here
     * because the JSON parser skips keys with value -1.0 that it cannot
     * distinguish from valid user intent.
     */
    void SetDefaults()
    {
        modelType = "";
        ecmStiffness = 5.0;
        runNumber = 0;
        verbosity = 1;

        // Unit conversion defaults (disabled by default for backward compatibility)
        usePhysicalUnits = false;
        cellDiameterMicrometers    = 10.0;   // Typical mammalian epithelial cell: 10 µm
        referenceForcenN           = 10.0;   // Typical epithelial traction force: ~10 nN
        referenceVelocityUmPerH    = 10.0;   // Typical migration speed: ~10 µm/h (= 1 CD/h)
        ecmElasticModulusPa        = -1.0;   // Sentinel: disabled (use ecmStiffness directly)

        enableLumenPressure = true;
        enableApicalConstriction = true;
        enableEcmGuidance = false;
        enableRelaxation = true;
        enableSloughing = false;
        enableDifferentialAdhesion = true;
        enableEcmConfinement = true;
        enableLumenHole = true;
        enableContinuousPvd = false;
        endTimeOverridden = false;
        dtOverridden = false;
        t1ThresholdOverridden = false;
        t2ThresholdOverridden = false;
        endTime = 168.0;
        dt = 0.005;


        relaxationTime = 10.0;

        organoidRadius2d    = 8.0;
        numCells2dNode      = 80;
        numCells2dVertex    = 40;
        innerRadius2d       = 6.0;
        outerRadius2d       = 8.0;
        t1Threshold2d       = 0.01;   // Chaste default; override via JSON "t1Threshold2d"
        t2Threshold2d       = 0.001;  // Chaste default; override via JSON "t2Threshold2d"
        interactionCutoff2d = 2.5;

        numCells3dNode       = 100;
        numCells3dVertex     = 200;
        organoidRadius3d     = 10.0;
        shellThickness3d     = 3.0;
        interactionCutoff3d  = 8.0;
        sphereRadius3dVertex = 10.0;  // derived from organoidRadius3d in Finalise()
        bmOffset3dVertex     = 1.0;

        // ECM degradation parameters
        ecmDegradationRate   = 0.02;   // Default ECM density loss per cell per unit time
        ecmDiffusionCoeff    = 0.1;    // ECM density smoothing coefficient

        lumenPressure          = 2.0;

        apicalConstrictionStrength = 0.2;

        springStiffness       = 30.0;
        springCutoff          = 1.5;
        springStiffnessTAScale   = 1.0;   // TA spring stiffness multiplier (relative to stem)
        springStiffnessDiffScale = 1.0;   // Differentiated/Paneth spring stiffness multiplier
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
        ecmGridType    = "hex";   // "square" or "hex"
        ecmSpringRestLength      = 0.0;    // Rest length for cell-ECM springs
        ecmInteractionCutoff     = 1.5;    // Cutoff distance for cell-ECM interactions

        // Ghost node ECM defaults
        enableGhostNodeECM          = true;
        ghostDamping                = 5.0;
        ghostRemovalThreshold       = 0.05;
        ghostFibreRemodelingRate    = 0.1;
        ghostAnisotropyStrength     = 0.5;
        ghostCellGhostRestLength    = 1.0;
        ghostCellGhostCutoff        = 1.5;
        ghostSpringRestLength       = -1.0;  // sentinel: auto-derive from grid spacing
        ghostGridSpacing            = -1.0;  // sentinel: defaults to ghostSpringRestLength (or ecmGridSpacing if rest length also unset)
        ghostRemovalCheckInterval   = 10;

        // Viscoelastic ECM (generalised Maxwell model)
        enableViscoelasticECM       = true;
        ghostRelaxedStiffness       = 5.0;    // E0: long-time modulus
        ghostRelaxationModulus      = 1.0;    // E1: transient modulus
        ghostRelaxationTime         = 1.0;    // tau: relaxation time (hours)

        // Dynamic ghost node boundary expansion
        enableGhostBoundaryExpansion = false;

        // Curvature bending force (Drasdo 2000 - monolayer enforcement)
        enableCurvatureBending    = true;   // Enable by default for node2d
        bendingStiffness          = 5.0;    // Bending rigidity

        // Spring neighbor strategy (node-based models)
        useTopologyBasedSprings   = true;   // true = topology (ring/surface), false = distance threshold

        // Cell type proportions for uniform random distribution
        stemFraction    = 0.2;
        transitFraction = 0.5;

        // Cell cycle duration (uniform random total cycle)
        stemCycleMin  = 12.0;   // U(12, 14) h total cycle for stem cells
        stemCycleMax  = 14.0;
        taCycleRatio  = 0.5;    // TA cycle = ratio × stem cycle (set 0.5 for half)

        // Generational cascade (Meineke et al. 2001)
        enableGenerationalCascade = false;  // Disabled by default (uses stochastic 4-type model)
        maxTransitGenerations = 3;          // TA cells differentiate after 3 divisions

        // Stochastic 4-type cell cycle (Montes-Olivas et al. 2023)
        enableStochasticFourType = true;    // Default: stochastic SC/TA/Paneth/EC fate model
        probStemToStem    = 0.89;           // P(SC daughter = SC)
        probStemToPaneth  = 0.09;           // P(SC daughter = PC)
        probTaToTaEarly   = 0.9;            // P(TA daughter = TA) days 1-4
        probTaToTaLate    = 0.7;            // P(TA daughter = TA) days 5-7
        transitionTime    = 120.0;          // day 5 = 120 hours
        panethFraction    = 0.09;           // Initial PC fraction (4-type model)

        // Cell polarity (ya||a-style monolayer enforcement)
        enableCellPolarity          = false;   // Enable by default for node models
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

        // For vertex2d, organoidRadius2d must match the actual mesh size,
        // not the node model size, so ECM domain and clear radius are correct.
        if (modelType == "vertex2d")
        {
            organoidRadius2d = outerRadius2d;
        }

        // Derive organoidRadius3d from numCells so Voronoi cell area ≈ 1.0²
        // 4πR² / N = 1.0  →  R = √(N / 4π)
        organoidRadius3d = std::sqrt(numCells3dNode / (4.0 * M_PI));

        bmStiffnessNode   = ecmStiffness;
        bmStiffnessVertex = ecmStiffness * 0.5;

        // Note: ecmDegradationRate and ecmDiffusionCoeff are loaded from config file
        // Only set defaults if not already set (SetDefaults handles this)

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

        // Resolve ghost node grid spacing: ghostGridSpacing > ghostSpringRestLength > ecmGridSpacing
        if (ghostGridSpacing < 0.0)
        {
            ghostGridSpacing = (ghostSpringRestLength > 0.0) ? ghostSpringRestLength : ecmGridSpacing;
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

        // Compute samplingMultiple so every simulation outputs exactly 50 frames
        const unsigned totalSteps = static_cast<unsigned>(std::round(endTime / dt));
        samplingMultiple = std::max(1u, totalSteps / 50);
    }

    /**
     * Cross-validate parameters after Finalise(). Throws EXCEPTION on
     * invalid combinations. Call after Finalise(), before entering a runner.
     */
    void Validate() const
    {
        // ── Model type must be recognised ────────────────────────────
        if (modelType != "node2d" && modelType != "vertex2d" &&
            modelType != "node3d" && modelType != "vertex3d")
        {
            EXCEPTION("Unknown modelType: \"" + modelType + "\". "
                      "Expected node2d, vertex2d, node3d, or vertex3d.");
        }

        // ── Positive time stepping ───────────────────────────────────
        if (dt <= 0.0)
            EXCEPTION("dt must be positive (got " + std::to_string(dt) + ")");
        if (endTime <= 0.0)
            EXCEPTION("endTime must be positive (got " + std::to_string(endTime) + ")");
        if (enableRelaxation && relaxationTime <= 0.0)
            EXCEPTION("enableRelaxation is true but relaxationTime <= 0");

        // ── Cell fraction sanity ─────────────────────────────────────
        double totalFraction = stemFraction + transitFraction;
        if (enableStochasticFourType)
            totalFraction += panethFraction;
        if (totalFraction > 1.0 + 1e-9)
            EXCEPTION("Cell type fractions sum to " + std::to_string(totalFraction)
                      + " (must be <= 1.0)");

        // ── Mutually exclusive cell cycle models ─────────────────────
        if (enableStochasticFourType && enableGenerationalCascade)
            EXCEPTION("enableStochasticFourType and enableGenerationalCascade "
                      "are mutually exclusive — set at most one to true.");

        // ── Cell cycle ranges ────────────────────────────────────────
        if (stemCycleMin > stemCycleMax)
            EXCEPTION("stemCycleMin (" + std::to_string(stemCycleMin)
                      + ") > stemCycleMax (" + std::to_string(stemCycleMax) + ")");

        // ── Ghost node ECM constraints ───────────────────────────────
        if (enableViscoelasticECM && !enableGhostNodeECM)
            EXCEPTION("enableViscoelasticECM requires enableGhostNodeECM=true");

        // ── Geometry sanity ──────────────────────────────────────────
        if ((modelType == "node2d" || modelType == "vertex2d") && numCells2dNode < 3 && modelType == "node2d")
            EXCEPTION("numCells2dNode must be >= 3 for a ring (got "
                      + std::to_string(numCells2dNode) + ")");
        if (modelType == "vertex2d" && numCells2dVertex < 3)
            EXCEPTION("numCells2dVertex must be >= 3 (got "
                      + std::to_string(numCells2dVertex) + ")");
        if (modelType == "node3d" && numCells3dNode < 4)
            EXCEPTION("numCells3dNode must be >= 4 for a sphere (got "
                      + std::to_string(numCells3dNode) + ")");

        // ── Stiffness/force parameters should be non-negative ────────
        if (ecmStiffness < 0.0)
            EXCEPTION("ecmStiffness must be >= 0 (got " + std::to_string(ecmStiffness) + ")");
        if (springStiffness < 0.0)
            EXCEPTION("springStiffness must be >= 0 (got " + std::to_string(springStiffness) + ")");

        std::cout << "Parameter validation passed." << std::endl;
    }

    /**
     * Print warnings for boolean feature flags that were NOT set in the config
     * file and therefore kept their (potentially unexpected) compiled-in default.
     * Only flags whose default is TRUE are flagged — these silently add forces.
     */
    void PrintDefaultsWarning() const
    {
        // Feature flags whose compiled-in default is TRUE (i.e. force/feature is ON
        // unless the user explicitly disables it in the config file).
        static const std::vector<std::pair<std::string, bool>> kFeatureDefaults = {
            {"enableLumenPressure",       true},
            {"enableApicalConstriction",  true},
            {"enableRelaxation",          true},
            {"enableDifferentialAdhesion",true},
            {"enableEcmConfinement",      true},
            {"enableLumenHole",           true},
            {"enableGhostNodeECM",        true},
            {"enableViscoelasticECM",     true},
            {"enableCurvatureBending",    true},
        };

        std::vector<std::string> warnings;
        for (const auto& kv : kFeatureDefaults)
        {
            if (mLoadedKeys.find(kv.first) == mLoadedKeys.end())
            {
                warnings.push_back(kv.first);
            }
        }

        if (!warnings.empty())
        {
            std::cerr << "\n  WARNING: The following feature flags were NOT set in your config\n"
                      << "  file and are using their compiled-in defaults (shown below).\n"
                      << "  Add them explicitly to suppress this warning.\n\n";
            for (const auto& name : warnings)
            {
                // Look up current (default) value
                bool val = false;
                if (name == "enableLumenPressure")        val = enableLumenPressure;
                else if (name == "enableApicalConstriction") val = enableApicalConstriction;
                else if (name == "enableRelaxation")      val = enableRelaxation;
                else if (name == "enableDifferentialAdhesion") val = enableDifferentialAdhesion;
                else if (name == "enableEcmConfinement")  val = enableEcmConfinement;
                else if (name == "enableLumenHole")       val = enableLumenHole;
                else if (name == "enableGhostNodeECM")    val = enableGhostNodeECM;
                else if (name == "enableViscoelasticECM") val = enableViscoelasticECM;
                else if (name == "enableCurvatureBending") val = enableCurvatureBending;
                std::cerr << "    " << std::left << std::setw(32) << name
                          << " [default: " << (val ? "true" : "false") << "]\n";
            }
            std::cerr << std::endl;
        }
    }

private:
    /**
     * Apply a flat key→value map to all parameters.
     * Used by both LoadFromFile (INI) and LoadFromJson.
     */
    void ApplyConfigMap(const std::map<std::string, std::string>& configMap)
    {
        auto getBool = [&](const std::string& key, bool& var) {
            if (configMap.count(key)) {
                std::string v = configMap.at(key);
                std::transform(v.begin(), v.end(), v.begin(), ::tolower);
                var = (v == "true" || v == "1" || v == "yes" || v == "on");
                mLoadedKeys.insert(key);
            }
        };
        auto getDouble = [&](const std::string& key, double& var) {
            if (configMap.count(key)) {
                try { var = std::stod(configMap.at(key)); mLoadedKeys.insert(key); }
                catch (...) { std::cerr << "Warning: Invalid double for " << key << std::endl; }
            }
        };
        auto getUnsigned = [&](const std::string& key, unsigned& var) {
            if (configMap.count(key)) {
                try { var = static_cast<unsigned>(std::stoul(configMap.at(key))); mLoadedKeys.insert(key); }
                catch (...) { std::cerr << "Warning: Invalid unsigned for " << key << std::endl; }
            }
        };
        auto getString = [&](const std::string& key, std::string& var) {
            if (configMap.count(key)) { var = configMap.at(key); mLoadedKeys.insert(key); }
        };

        // ── Simulation ───────────────────────────────────────────────
        getString("modelType", modelType);
        getDouble("ecmStiffness", ecmStiffness);
        getUnsigned("runNumber", runNumber);
        getUnsigned("randomSeed", randomSeed);

        // ── Unit conversion (load first, before stiffness/pressure parameters) ──
        getBool("usePhysicalUnits", usePhysicalUnits);
        getDouble("cellDiameterMicrometers", cellDiameterMicrometers);
        getDouble("referenceForcenN", referenceForcenN);
        getDouble("referenceVelocityUmPerH", referenceVelocityUmPerH);
        getDouble("ecmElasticModulusPa", ecmElasticModulusPa);

        // ── Feature toggles ──────────────────────────────────────────
        getBool("enableLumenPressure", enableLumenPressure);
        getBool("enableApicalConstriction", enableApicalConstriction);
        getBool("enableEcmGuidance", enableEcmGuidance);
        getBool("enableRelaxation", enableRelaxation);
        getBool("enableSloughing", enableSloughing);
        getBool("enableDifferentialAdhesion", enableDifferentialAdhesion);
        getBool("enableCurvatureBending", enableCurvatureBending);
        getBool("enableCellPolarity", enableCellPolarity);
        getBool("enableEcmConfinement", enableEcmConfinement);
        getBool("enableLumenHole", enableLumenHole);
        getBool("useTopologyBasedSprings", useTopologyBasedSprings);
        getBool("enableContinuousPvd", enableContinuousPvd);
        getBool("enableGhostNodeECM", enableGhostNodeECM);
        getBool("enableViscoelasticECM", enableViscoelasticECM);
        getBool("enableGhostBoundaryExpansion", enableGhostBoundaryExpansion);

        // ── Time stepping ────────────────────────────────────────────
        if (configMap.count("dt"))      { getDouble("dt", dt);           dtOverridden = true; }
        if (configMap.count("endTime")) { getDouble("endTime", endTime); endTimeOverridden = true; }
        getDouble("relaxationTime", relaxationTime);

        // ── Geometry ─────────────────────────────────────────────────
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

        getDouble("bmRadiusFraction", bmRadiusFraction);
        getDouble("ecmMaxRadiusFraction", ecmMaxRadiusFraction);

        // ── Forces ───────────────────────────────────────────────────
        getDouble("lumenPressure", lumenPressure);
        getDouble("apicalConstrictionStrength", apicalConstrictionStrength);

        getDouble("springStiffness", springStiffness);
        getDouble("springCutoff", springCutoff);
        getDouble("springStiffnessTAScale", springStiffnessTAScale);
        getDouble("springStiffnessDiffScale", springStiffnessDiffScale);
        getDouble("apicalApicalAdhesion", apicalApicalAdhesion);
        getDouble("basalBasalAdhesion", basalBasalAdhesion);
        getDouble("apicalBasalAdhesion", apicalBasalAdhesion);

        getDouble("nhMembraneSurface", nhMembraneSurface);
        getDouble("nhCellCellAdhesion", nhCellCellAdhesion);
        getDouble("nhBoundaryAdhesion", nhBoundaryAdhesion);
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

        getDouble("bendingStiffness", bendingStiffness);
        getDouble("polarityBendingStrength", polarityBendingStrength);
        getDouble("polarityAlignmentStrength", polarityAlignmentStrength);

        // ── ECM confinement ──────────────────────────────────────────
        getDouble("ecmDomainHalf", ecmDomainHalf);
        getDouble("ecmGridSpacing", ecmGridSpacing);
        getDouble("ecmBaseSpeed", ecmBaseSpeed);
        getString("ecmGridType", ecmGridType);
        getDouble("ecmSpringRestLength", ecmSpringRestLength);
        getDouble("ecmInteractionCutoff", ecmInteractionCutoff);

        // ── Ghost node ECM ───────────────────────────────────────────
        getDouble("ghostDamping", ghostDamping);
        getDouble("ghostRemovalThreshold", ghostRemovalThreshold);
        getDouble("ghostFibreRemodelingRate", ghostFibreRemodelingRate);
        getDouble("ghostAnisotropyStrength", ghostAnisotropyStrength);
        getDouble("ghostCellGhostRestLength", ghostCellGhostRestLength);
        getDouble("ghostCellGhostCutoff", ghostCellGhostCutoff);
        getDouble("ghostSpringRestLength", ghostSpringRestLength);
        getDouble("ghostGridSpacing", ghostGridSpacing);
        getUnsigned("ghostRemovalCheckInterval", ghostRemovalCheckInterval);

        // ── Viscoelastic ECM ─────────────────────────────────────────
        getDouble("ghostRelaxedStiffness", ghostRelaxedStiffness);
        getDouble("ghostRelaxationModulus", ghostRelaxationModulus);
        getDouble("ghostRelaxationTime", ghostRelaxationTime);

        // ── Vertex mesh thresholds ───────────────────────────────────
        if (configMap.count("t1Threshold2d")) { getDouble("t1Threshold2d", t1Threshold2d); t1ThresholdOverridden = true; }
        if (configMap.count("t2Threshold2d")) { getDouble("t2Threshold2d", t2Threshold2d); t2ThresholdOverridden = true; }

        // ── Cell types & cycle ───────────────────────────────────────
        getDouble("stemFraction", stemFraction);
        getDouble("transitFraction", transitFraction);
        getDouble("quiescentFraction", quiescentFraction);
        getDouble("sloughRadiusFactor", sloughRadiusFactor);
        getDouble("stemCycleMin", stemCycleMin);
        getDouble("stemCycleMax", stemCycleMax);
        getDouble("taCycleRatio", taCycleRatio);

        getBool("enableGenerationalCascade", enableGenerationalCascade);
        getUnsigned("maxTransitGenerations", maxTransitGenerations);

        getBool("enableStochasticFourType", enableStochasticFourType);
        getDouble("probStemToStem", probStemToStem);
        getDouble("probStemToPaneth", probStemToPaneth);
        getDouble("probTaToTaEarly", probTaToTaEarly);
        getDouble("probTaToTaLate", probTaToTaLate);
        getDouble("transitionTime", transitionTime);
        getDouble("panethFraction", panethFraction);

        // Apply unit conversion if enabled (must be called after all parameters loaded)
        ComputeConversionFactors();
    }

public:
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

        ApplyConfigMap(configMap);
        std::cout << "Loaded " << configMap.size() << " parameters from: " << filePath << std::endl;
        return true;
    }

    /**
     * Load parameters from a JSON config file.
     * The JSON uses a nested hierarchy (simulation/features/geometry/forces/...)
     * for clarity, but all leaf key names match the flat parameter names used
     * throughout the codebase. Object nesting is ignored during parsing —
     * only the innermost key→value pairs are applied.
     *
     * Only parameters present in the file are overwritten; others keep defaults.
     * Returns true if file was loaded successfully.
     */
    bool LoadFromJson(const std::string& filePath)
    {
        std::ifstream ifs(filePath);
        if (!ifs.is_open())
        {
            std::cerr << "Warning: Could not open JSON config file: " << filePath << std::endl;
            return false;
        }
        const std::string content((std::istreambuf_iterator<char>(ifs)),
                                   std::istreambuf_iterator<char>());
        ifs.close();

        // ── Minimal recursive-descent JSON parser ─────────────────────
        // Flattens all leaf (scalar) key→value pairs from any nesting depth
        // into a single map. Object-level keys (e.g. "forces", "SpringForce")
        // are skipped; only leaf values reach the map.
        std::map<std::string, std::string> configMap;
        size_t pos = 0;

        auto skipWs = [&]() {
            while (pos < content.size() &&
                   (content[pos]==' ' || content[pos]=='\t' ||
                    content[pos]=='\n' || content[pos]=='\r'))
                ++pos;
        };

        std::function<std::string()> parseString = [&]() -> std::string {
            ++pos; // skip opening "
            std::string result;
            while (pos < content.size() && content[pos] != '"')
            {
                if (content[pos] == '\\' && pos + 1 < content.size())
                    { ++pos; result += content[pos++]; }
                else
                    result += content[pos++];
            }
            if (pos < content.size()) ++pos; // skip closing "
            return result;
        };

        std::function<void()> parseObject = [&]() {
            skipWs();
            if (pos < content.size() && content[pos] == '{') ++pos;
            while (true)
            {
                skipWs();
                if (pos >= content.size() || content[pos] == '}')
                    { if (pos < content.size()) ++pos; break; }
                if (content[pos] == ',') { ++pos; continue; }
                if (content[pos] != '"') { ++pos; continue; } // malformed

                const std::string key = parseString();
                skipWs();
                if (pos < content.size() && content[pos] == ':') ++pos;
                skipWs();

                if (pos < content.size() && content[pos] == '{')
                {
                    parseObject(); // recurse — leaf keys bubble up into configMap
                }
                else
                {
                    std::string value;
                    if (pos < content.size() && content[pos] == '"')
                    {
                        value = parseString();
                    }
                    else
                    {
                        const size_t start = pos;
                        while (pos < content.size() &&
                               content[pos] != ',' && content[pos] != '}' &&
                               content[pos] != '\n' && content[pos] != '\r')
                            ++pos;
                        value = content.substr(start, pos - start);
                        const size_t end = value.find_last_not_of(" \t");
                        value = (end == std::string::npos) ? "" : value.substr(0, end + 1);
                    }
                    if (!key.empty() && !value.empty() && value != "null")
                        configMap[key] = value;
                }
            }
        };

        skipWs();
        parseObject();
        // ─────────────────────────────────────────────────────────────

        ApplyConfigMap(configMap);
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
        file << "relaxationTime = " << relaxationTime << "    # Relaxation phase duration\n";
        file << "runNumber = " << runNumber << "\n";
        file << "randomSeed = " << randomSeed << "\n\n";

        file << "[UnitConversion]\n";
        file << "usePhysicalUnits = " << (usePhysicalUnits ? "true" : "false") << "  # Enable/disable physical unit conversion\n";
        file << "cellDiameterMicrometers = " << cellDiameterMicrometers << "  # Cell diameter in micrometers\n";
        file << "referenceForcenN = " << referenceForcenN << "  # Reference traction force in nN\n";
        file << "referenceVelocityUmPerH = " << referenceVelocityUmPerH << "  # Reference migration speed in um/h\n\n";

        file << "[Features]\n";
        file << "enableLumenPressure = " << (enableLumenPressure ? "true" : "false") << "\n";
        file << "enableApicalConstriction = " << (enableApicalConstriction ? "true" : "false") << "\n";
        file << "enableEcmGuidance = " << (enableEcmGuidance ? "true" : "false") << "\n";
        file << "enableRelaxation = " << (enableRelaxation ? "true" : "false") << "\n";
        file << "enableSloughing = " << (enableSloughing ? "true" : "false") << "\n";
        file << "enableDifferentialAdhesion = " << (enableDifferentialAdhesion ? "true" : "false") << "\n";
        file << "enableCurvatureBending = " << (enableCurvatureBending ? "true" : "false") << "\n";
        file << "enableCellPolarity = " << (enableCellPolarity ? "true" : "false") << "\n";
        file << "enableLumenHole = " << (enableLumenHole ? "true" : "false") << "  # Clear ECM density inside organoid (creates central lumen)\n";
        file << "useTopologyBasedSprings = " << (useTopologyBasedSprings ? "true" : "false") << "  # true=topology (ring/surface), false=distance threshold\n";
        file << "enableContinuousPvd = " << (enableContinuousPvd ? "true" : "false") << "  # Keep .pvd files valid during simulation\n\n";

        file << "[ECM]\n";
        file << "ecmStiffness = " << ecmStiffness << "  # ECM spring stiffness (on-lattice: cell-ECM; off-lattice: ECM-ECM node)\n";
        file << "ecmSpringRestLength = " << ecmSpringRestLength << "  # Rest length for cell-ECM springs\n";
        file << "ecmInteractionCutoff = " << ecmInteractionCutoff << "  # Cutoff distance for cell-ECM interactions\n";
        file << "bmStiffnessNode = " << bmStiffnessNode << "  # 3D only: BasementMembraneForce stiffness (derived from ecmStiffness)\n";
        file << "bmStiffnessVertex = " << bmStiffnessVertex << "  # 3D only: BasementMembraneForce stiffness (= ecmStiffness × 0.5)\n";
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
        file << "springStiffnessTAScale = " << springStiffnessTAScale << "\n";
        file << "springStiffnessDiffScale = " << springStiffnessDiffScale << "\n";
        file << "apicalApicalAdhesion = " << apicalApicalAdhesion << "\n";
        file << "basalBasalAdhesion = " << basalBasalAdhesion << "\n";
        file << "apicalBasalAdhesion = " << apicalBasalAdhesion << "\n\n";

        file << "# Lumen pressure\n";
        file << "lumenPressure = " << lumenPressure << "\n\n";

        file << "# Apical constriction\n";
        file << "apicalConstrictionStrength = " << apicalConstrictionStrength << "\n\n";

        file << "# Curvature bending (Drasdo 2000 - monolayer enforcement)\n";
        file << "bendingStiffness = " << bendingStiffness << "\n\n";

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
        file << "maxTransitGenerations = " << maxTransitGenerations << "       # TA divisions before differentiation\n";
        file << "enableStochasticFourType = " << (enableStochasticFourType ? "true" : "false") << "  # Stochastic SC/TA/PC/EC transitions (Montes-Olivas 2023)\n";
        file << "probStemToStem = " << probStemToStem << "       # P(SC daughter = SC)\n";
        file << "probStemToPaneth = " << probStemToPaneth << "     # P(SC daughter = PC)\n";
        file << "probTaToTaEarly = " << probTaToTaEarly << "      # P(TA daughter = TA) days 1-4\n";
        file << "probTaToTaLate = " << probTaToTaLate << "       # P(TA daughter = TA) days 5-7\n";
        file << "transitionTime = " << transitionTime << "       # Switch time (hours)\n";
        file << "panethFraction = " << panethFraction << "       # Initial PC fraction\n\n";

        file << "[ECMGuidance]\n";
        file << "ecmDomainHalf = " << ecmDomainHalf << "\n";
        file << "ecmGridSpacing = " << ecmGridSpacing << "\n";
        file << "ecmBaseSpeed = " << ecmBaseSpeed << "\n";
        file << "ecmGridType = " << ecmGridType << "\n";
        file << "ecmSpringRestLength = " << ecmSpringRestLength << "   # Rest length for cell-ECM springs\n";
        file << "ecmInteractionCutoff = " << ecmInteractionCutoff << "  # Cutoff for cell-ECM interactions\n";

        file << "[GhostNodeECM]\n";
        file << "enableGhostNodeECM = " << (enableGhostNodeECM ? "true" : "false") << "  # Use off-lattice ghost node ECM instead of grid-based\n";
        file << "ghostDamping = " << ghostDamping << "                # ECM node damping coefficient\n";
        file << "ghostRemovalThreshold = " << ghostRemovalThreshold << "       # Density below which ghost nodes are removed\n";
        file << "ghostFibreRemodelingRate = " << ghostFibreRemodelingRate << "    # Fibre reorientation rate\n";
        file << "ghostAnisotropyStrength = " << ghostAnisotropyStrength << "     # Fibre anisotropy factor [0,1]\n";
        file << "ghostCellGhostRestLength = " << ghostCellGhostRestLength << "    # Cell-ghost spring rest length\n";
        file << "ghostCellGhostCutoff = " << ghostCellGhostCutoff << "        # Cell-ghost interaction cutoff\n";
        file << "ghostSpringRestLength = " << ghostSpringRestLength << "       # ECM-ECM spring rest length (-1 = auto)\n";
        file << "ghostGridSpacing = " << ghostGridSpacing << "            # Ghost node grid spacing (-1 = auto from rest length)\n";
        file << "ghostRemovalCheckInterval = " << ghostRemovalCheckInterval << "   # Steps between removal checks\n";
        file << "# Viscoelastic constitutive model (generalised Maxwell)\n";
        file << "enableViscoelasticECM = " << (enableViscoelasticECM ? "true" : "false") << "  # Use viscoelastic ghost node ECM\n";
        file << "ghostRelaxedStiffness = " << ghostRelaxedStiffness << "       # E0: relaxed modulus\n";
        file << "ghostRelaxationModulus = " << ghostRelaxationModulus << "      # E1: transient modulus\n";
        file << "ghostRelaxationTime = " << ghostRelaxationTime << "          # tau: relaxation time (hours)\n";
        file << "enableGhostBoundaryExpansion = " << (enableGhostBoundaryExpansion ? "true" : "false") << "  # Dynamically expand ghost domain when cells approach boundary\n";

        file.close();
        std::cout << "Saved parameters to: " << filePath << std::endl;
        return true;
    }

    /**
     * Save current parameters to a hierarchical JSON config file.
     * Parameters are grouped by the force/subsystem they belong to,
     * mirroring the structure of the hand-authored JSON param files.
     */
    bool SaveToJson(const std::string& filePath) const
    {
        std::ofstream f(filePath);
        if (!f.is_open())
        {
            std::cerr << "Error: Could not write JSON config file: " << filePath << std::endl;
            return false;
        }

        auto b = [](bool v) -> const char* { return v ? "true" : "false"; };

        f << "{\n";

        // ── Simulation ───────────────────────────────────────────────
        f << "  \"simulation\": {\n";
        f << "    \"modelType\": \"" << modelType << "\",\n";
        f << "    \"endTime\": " << endTime << ",\n";
        f << "    \"dt\": " << dt << ",\n";
        f << "    \"relaxationTime\": " << relaxationTime << ",\n";
        f << "    \"runNumber\": " << runNumber << ",\n";
        f << "    \"randomSeed\": " << randomSeed << "\n";
        f << "  },\n";

        // ── Unit conversion ──────────────────────────────────────────
        f << "  \"unitConversion\": {\n";
        f << "    \"usePhysicalUnits\": " << b(usePhysicalUnits) << ",\n";
        f << "    \"cellDiameterMicrometers\": " << cellDiameterMicrometers << ",\n";
        f << "    \"referenceForcenN\": " << referenceForcenN << ",\n";
        f << "    \"referenceVelocityUmPerH\": " << referenceVelocityUmPerH << "\n";
        f << "  },\n";

        // ── Feature toggles ──────────────────────────────────────────
        f << "  \"features\": {\n";
        f << "    \"enableLumenPressure\": " << b(enableLumenPressure) << ",\n";
        f << "    \"enableApicalConstriction\": " << b(enableApicalConstriction) << ",\n";
        f << "    \"enableEcmGuidance\": " << b(enableEcmGuidance) << ",\n";
        f << "    \"enableRelaxation\": " << b(enableRelaxation) << ",\n";
        f << "    \"enableSloughing\": " << b(enableSloughing) << ",\n";
        f << "    \"enableDifferentialAdhesion\": " << b(enableDifferentialAdhesion) << ",\n";
        f << "    \"enableCurvatureBending\": " << b(enableCurvatureBending) << ",\n";
        f << "    \"enableCellPolarity\": " << b(enableCellPolarity) << ",\n";
        f << "    \"enableEcmConfinement\": " << b(enableEcmConfinement) << ",\n";
        f << "    \"enableLumenHole\": " << b(enableLumenHole) << ",\n";
        f << "    \"enableGhostNodeECM\": " << b(enableGhostNodeECM) << ",\n";
        f << "    \"enableViscoelasticECM\": " << b(enableViscoelasticECM) << ",\n";
        f << "    \"enableGhostBoundaryExpansion\": " << b(enableGhostBoundaryExpansion) << ",\n";
        f << "    \"enableContinuousPvd\": " << b(enableContinuousPvd) << ",\n";
        f << "    \"useTopologyBasedSprings\": " << b(useTopologyBasedSprings) << "\n";
        f << "  },\n";

        // ── Geometry ─────────────────────────────────────────────────
        f << "  \"geometry\": {\n";
        f << "    \"numCells2dNode\": " << numCells2dNode << ",\n";
        f << "    \"numCells2dVertex\": " << numCells2dVertex << ",\n";
        f << "    \"interactionCutoff2d\": " << interactionCutoff2d << ",\n";
        f << "    \"t1Threshold2d\": " << t1Threshold2d << ",\n";
        f << "    \"t2Threshold2d\": " << t2Threshold2d << ",\n";
        f << "    \"numCells3dNode\": " << numCells3dNode << ",\n";
        f << "    \"numCells3dVertex\": " << numCells3dVertex << ",\n";
        f << "    \"shellThickness3d\": " << shellThickness3d << ",\n";
        f << "    \"interactionCutoff3d\": " << interactionCutoff3d << ",\n";
        f << "    \"bmOffset3dVertex\": " << bmOffset3dVertex << "\n";
        f << "  },\n";

        // ── Forces ───────────────────────────────────────────────────
        f << "  \"forces\": {\n";

        f << "    \"SpringForce\": {\n";
        f << "      \"springStiffness\": " << springStiffness << ",\n";
        f << "      \"springCutoff\": " << springCutoff << ",\n";
        f << "      \"springStiffnessTAScale\": " << springStiffnessTAScale << ",\n";
        f << "      \"springStiffnessDiffScale\": " << springStiffnessDiffScale << "\n";
        f << "    },\n";

        f << "    \"LumenPressureForce\": {\n";
        f << "      \"lumenPressure\": " << lumenPressure << "\n";
        f << "    },\n";

        f << "    \"ApicalConstrictionForce\": {\n";
        f << "      \"apicalConstrictionStrength\": " << apicalConstrictionStrength << "\n";
        f << "    },\n";

        f << "    \"RingSmoothingForce\": {\n";
        f << "      \"bendingStiffness\": " << bendingStiffness << ",\n";
        f << "      \"gammaStemScale\": " << gammaStemScale << ",\n";
        f << "      \"gammaTransitScale\": " << gammaTransitScale << ",\n";
        f << "      \"gammaDiffScale\": " << gammaDiffScale << "\n";
        f << "    },\n";

        f << "    \"CellPolarityForce\": {\n";
        f << "      \"polarityBendingStrength\": " << polarityBendingStrength << ",\n";
        f << "      \"polarityAlignmentStrength\": " << polarityAlignmentStrength << "\n";
        f << "    },\n";

        f << "    \"ECMConfinementForce\": {\n";
        if (ecmElasticModulusPa > 0.0)
            f << "      \"ecmElasticModulusPa\": " << ecmElasticModulusPa << ",\n";
        f << "      \"ecmStiffness\": " << ecmStiffness << ",\n";
        f << "      \"ecmSpringRestLength\": " << ecmSpringRestLength << ",\n";
        f << "      \"ecmInteractionCutoff\": " << ecmInteractionCutoff << ",\n";
        f << "      \"bmRadiusFraction\": " << bmRadiusFraction << ",\n";
        f << "      \"ecmMaxRadiusFraction\": " << ecmMaxRadiusFraction << ",\n";
        f << "      \"ecmDomainHalf\": " << ecmDomainHalf << ",\n";
        f << "      \"ecmGridSpacing\": " << ecmGridSpacing << ",\n";
        f << "      \"ecmGridType\": \"" << ecmGridType << "\",\n";
        f << "      \"ecmDegradationRate\": " << ecmDegradationRate << ",\n";
        f << "      \"ecmDiffusionCoeff\": " << ecmDiffusionCoeff << "\n";
        f << "    },\n";

        f << "    \"ECMGuidanceForce\": {\n";
        f << "      \"ecmBaseSpeed\": " << ecmBaseSpeed << "\n";
        f << "    },\n";

        f << "    \"GhostNodeECM\": {\n";
        f << "      \"ghostDamping\": " << ghostDamping << ",\n";
        f << "      \"ghostRemovalThreshold\": " << ghostRemovalThreshold << ",\n";
        f << "      \"ghostFibreRemodelingRate\": " << ghostFibreRemodelingRate << ",\n";
        f << "      \"ghostAnisotropyStrength\": " << ghostAnisotropyStrength << ",\n";
        f << "      \"ghostCellGhostRestLength\": " << ghostCellGhostRestLength << ",\n";
        f << "      \"ghostCellGhostCutoff\": " << ghostCellGhostCutoff << ",\n";
        f << "      \"ghostSpringRestLength\": " << ghostSpringRestLength << ",\n";
        f << "      \"ghostGridSpacing\": " << ghostGridSpacing << ",\n";
        f << "      \"ghostRemovalCheckInterval\": " << ghostRemovalCheckInterval << ",\n";
        f << "      \"ViscoelasticECM\": {\n";
        f << "        \"ghostRelaxedStiffness\": " << ghostRelaxedStiffness << ",\n";
        f << "        \"ghostRelaxationModulus\": " << ghostRelaxationModulus << ",\n";
        f << "        \"ghostRelaxationTime\": " << ghostRelaxationTime << "\n";
        f << "      }\n";
        f << "    },\n";

        f << "    \"NagaiHondaForce\": {\n";
        f << "      \"nhMembraneSurface\": " << nhMembraneSurface << ",\n";
        f << "      \"nhCellCellAdhesion\": " << nhCellCellAdhesion << ",\n";
        f << "      \"nhBoundaryAdhesion\": " << nhBoundaryAdhesion << ",\n";
        f << "      \"nhStemStemAdhesion\": " << nhStemStemAdhesion << ",\n";
        f << "      \"nhStemTransitAdhesion\": " << nhStemTransitAdhesion << ",\n";
        f << "      \"nhStemDiffAdhesion\": " << nhStemDiffAdhesion << ",\n";
        f << "      \"nhTransitTransitAdhesion\": " << nhTransitTransitAdhesion << ",\n";
        f << "      \"nhTransitDiffAdhesion\": " << nhTransitDiffAdhesion << ",\n";
        f << "      \"nhDiffDiffAdhesion\": " << nhDiffDiffAdhesion << ",\n";
        f << "      \"nhStemBoundaryAdhesion\": " << nhStemBoundaryAdhesion << ",\n";
        f << "      \"nhTransitBoundaryAdhesion\": " << nhTransitBoundaryAdhesion << ",\n";
        f << "      \"nhDiffBoundaryAdhesion\": " << nhDiffBoundaryAdhesion << ",\n";
        f << "      \"gammaApical\": " << gammaApical << ",\n";
        f << "      \"gammaBasal\": " << gammaBasal << ",\n";
        f << "      \"gammaLateral\": " << gammaLateral << "\n";
        f << "    },\n";

        f << "    \"DifferentialAdhesionForce\": {\n";
        f << "      \"apicalApicalAdhesion\": " << apicalApicalAdhesion << ",\n";
        f << "      \"basalBasalAdhesion\": " << basalBasalAdhesion << ",\n";
        f << "      \"apicalBasalAdhesion\": " << apicalBasalAdhesion << "\n";
        f << "    }\n";

        f << "  },\n";

        // ── Cell types ───────────────────────────────────────────────
        f << "  \"cellTypes\": {\n";
        f << "    \"stemFraction\": " << stemFraction << ",\n";
        f << "    \"transitFraction\": " << transitFraction << "\n";
        f << "  },\n";

        // ── Cell cycle ───────────────────────────────────────────────
        f << "  \"cellCycle\": {\n";
        f << "    \"quiescentFraction\": " << quiescentFraction << ",\n";
        f << "    \"sloughRadiusFactor\": " << sloughRadiusFactor << ",\n";
        f << "    \"stemCycleMin\": " << stemCycleMin << ",\n";
        f << "    \"stemCycleMax\": " << stemCycleMax << ",\n";
        f << "    \"taCycleRatio\": " << taCycleRatio << ",\n";
        f << "    \"StochasticFourTypeModel\": {\n";
        f << "      \"enableStochasticFourType\": " << b(enableStochasticFourType) << ",\n";
        f << "      \"enableGenerationalCascade\": " << b(enableGenerationalCascade) << ",\n";
        f << "      \"probStemToStem\": " << probStemToStem << ",\n";
        f << "      \"probStemToPaneth\": " << probStemToPaneth << ",\n";
        f << "      \"probTaToTaEarly\": " << probTaToTaEarly << ",\n";
        f << "      \"probTaToTaLate\": " << probTaToTaLate << ",\n";
        f << "      \"transitionTime\": " << transitionTime << ",\n";
        f << "      \"panethFraction\": " << panethFraction << "\n";
        f << "    }\n";
        f << "  }\n";

        f << "}\n";
        f.close();
        std::cout << "Saved parameters to: " << filePath << std::endl;
        return true;
    }
};

#endif // CRYPTBUDDINGPARAMS_HPP_

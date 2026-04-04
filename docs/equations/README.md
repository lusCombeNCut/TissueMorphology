# Mathematical Reference

Constitutive equations, derivations, and numerical analysis for the TissueMorphology simulation framework.

| Document | Contents |
|----------|----------|
| [ConstitutiveEquations_Forces.md](ConstitutiveEquations_Forces.md) | All 11 force models: spring, Nagai-Honda, lumen pressure, ECM confinement, ghost node, viscoelastic, apical constriction, basement membrane, differential adhesion, contact guidance |
| [ConstitutiveEquations_CellCycles.md](ConstitutiveEquations_CellCycles.md) | Cell cycle models: uniform contact inhibition, generational cascade, stochastic four-type (Montes-Olivas 2023) |
| [LumenPressureForceDerivation.md](LumenPressureForceDerivation.md) | Work-energy lumen pressure derivation for all four model types (node/vertex, 2D/3D) |
| [ECM_Forces_Comparison.md](ECM_Forces_Comparison.md) | Side-by-side comparison of grid ECM vs. dynamic ECM vs. ghost node ECM |
| [Appendix_ViscoelasticECM_NumericalConsistency.md](Appendix_ViscoelasticECM_NumericalConsistency.md) | Proof of consistency and unconditional stability of the SLS discretisation |
| [derivation.md](derivation.md) | Integrating factor solution for viscoelastic rest-length ODE |

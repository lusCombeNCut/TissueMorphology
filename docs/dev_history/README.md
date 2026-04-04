# Development History

Historical documents recording the biological motivation and design evolution of the TissueMorphology project. These files are **archived** — they may reference deprecated classes or parameter names. Use the [equations/](../equations/) folder for current mathematical reference.

## Files

| Document | Status | Summary |
|----------|--------|---------|
| [THREE_FACTOR_MODEL.md](THREE_FACTOR_MODEL.md) | Archived | Original three-factor morphogenesis hypothesis (apical constriction, differential adhesion, BM stiffness). Superseded by current 10+ force architecture. |
| [CELL_PROLIFERATION_IMPLEMENTATION.md](CELL_PROLIFERATION_IMPLEMENTATION.md) | Archived | Early contact inhibition implementation notes. Superseded by `StochasticFourTypeCellCycleModel` — see [equations/ConstitutiveEquations_CellCycles.md](../equations/ConstitutiveEquations_CellCycles.md). |
| [INVASIVE_FRONT.md](INVASIVE_FRONT.md) | Archived | Invasive front experiments with Metzcar et al. ECM fibre model. Incomplete — see thesis Section 2.3.1 for the validated replication. |

## Key design decisions

### Why off-lattice ECM over on-lattice (Section 2.3)
The initial Metzcar et al. [16] on-lattice ECM was validated against Painter et al. cell invasion (thesis Fig 2.4), but showed grid-bias artefacts and could not model viscoelastic flow-back. This motivated the ghost node off-lattice model (Section 2.3.3).

### Why Simple Linear Solid over Kelvin-Voigt or higher-order GM
Fertala et al. (2025) showed a single Prony term suffices at the cell length scale (~10 µm), with minimal error vs. experiment. Higher-order GM adds computational cost without clear gain at this scale. See thesis Section 2.3.3.

### Why four model types
Node-based and vertex-based models have different strengths (node: simpler, faster; vertex: well-defined cell boundaries, T1/T2 transitions). Running both in 2D and 3D ensures results are not artefacts of a single methodology. See thesis Section 2.1.

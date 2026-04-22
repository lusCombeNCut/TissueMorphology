# Development History

Archived design rationale. For current mathematical reference see [equations/](../equations/).

## Key design decisions

### Why off-lattice ECM over on-lattice (Section 2.3)
The initial Metzcar et al. [16] on-lattice ECM was validated against Painter et al. cell invasion (thesis Fig 2.4), but showed grid-bias artefacts and could not model viscoelastic flow-back. This motivated the ghost node off-lattice model (Section 2.3.3).

### Why Simple Linear Solid over Kelvin-Voigt or higher-order GM
Fertala et al. (2025) showed a single Prony term suffices at the cell length scale (~10 µm), with minimal error vs. experiment. Higher-order GM adds computational cost without clear gain at this scale. See thesis Section 2.3.3.

### Why four model types
Node-based and vertex-based models have different strengths (node: simpler, faster; vertex: well-defined cell boundaries, T1/T2 transitions). Running both in 2D and 3D ensures results are not artefacts of a single methodology. See thesis Section 2.1.

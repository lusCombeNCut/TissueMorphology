# Parameter Justification and Literature Sources

This document records the experimental or theoretical basis for simulation parameters,
and **suggests what values should be used** where the current values lack justification.
Where no direct experimental value exists, this is stated explicitly.

Sources are verified by reading the actual source code (`NagaiHondaForce.cpp`,
`SurfaceTensionForce.cpp`, `RunVertex2d.hpp`, `RunVertex3d.hpp`) rather than
trusting documentation or the JSON files alone.

---

## 1. Nagai–Honda Force (vertex2d: `FastNagaiHondaForce`)

The energy functional per cell α:

$$U_\alpha = \lambda_d(A_\alpha - A_\alpha^0)^2 + \lambda_m(L_\alpha - L_\alpha^0)^2 + \sum_{\langle ij\rangle} \gamma_{ij}\,\ell_{ij}$$

### λ_d — Area (deformation) elasticity

**Hardcoded value:** `100.0` in [RunVertex2d.hpp:135](../apps/src/RunVertex2d.hpp)
```
p_nh->SetDeformationEnergyParameter(100.0);
```
The source comment reads: *"Deformation energy should be ~100 (Nagai-Honda default) for area
stability. cellGhostSpringStiffness is used for ECM field stiffness, not the vertex deformation energy."*

**Source:** `NagaiHondaForce.cpp` line 41 confirms: *"This is 1.0 in the Nagai & Honda paper."*
The Chaste convention scales all NH parameters up by ×100 relative to the original
Nagai & Honda (2001) paper (*Philosophical Magazine B*, 81(7), 699–719.
doi:[10.1080/13642810108205772](https://doi.org/10.1080/13642810108205772)).
**λ_d = 100 is not set by cellGhostSpringStiffness** — the ConstitutiveEquations_Forces.md is wrong
on this point. **Do not change λ_d unless you want to leave the standard Chaste convention.**

---

### λ_m — Perimeter/membrane elasticity (`nhMembraneSurface`)

**Current JSON value:** `10.0`
**Chaste source default:** `10.0` (matches; `NagaiHondaForce.cpp` line 42: *"This is 0.1 in the Nagai & Honda paper"* — i.e., 100× paper value, consistent with λ_d scaling)

**Literature reference:** Farhadifar et al. (2007), *Current Biology* 17(24), 2095–2104.
doi:[10.1016/j.cub.2007.11.049](https://doi.org/10.1016/j.cub.2007.11.049)

Farhadifar fitted vertex model parameters to Drosophila wing disc images combined with
laser ablation (direct measurement of junction recoil). With their normalisation K_A = 1:

| Regime | Γ (perimeter) | Λ (line tension) |
|--------|--------------|-----------------|
| Soft   | 0.04         | 0.12            |
| Stiff  | 0.02         | 0.0             |

**Scaling to CryptBudding (K_A = 100):**
- Soft: Γ → **4.0**, Λ → **12.0**
- Stiff: Γ → **2.0**, Λ → **0.0**

**Assessment of current value:** `nhMembraneSurface = 10.0` is **2.5× above** the Farhadifar
soft-regime value of 4.0, placing the tissue in an unusually high-perimeter-elasticity regime.
High λ_m causes cells to strongly resist perimeter changes (very regular hexagons, less
responsive to force imbalances).

**Suggested value:** `nhMembraneSurface = 4.0`
(Farhadifar soft regime, appropriate for a proliferating epithelium that can rearrange.)
If the current simulations produce stable crypt geometry, this change will make cells
slightly less regular but more biologically realistic.

---

### γ_cc — Cell–cell line tension (`nhCellCellAdhesion`)

**Current JSON value:** `10.0`
**Chaste source default:** `0.5` (`NagaiHondaForce.cpp` line 43: *"corresponds to σ = 1.0 in
the NH paper; in the paper σ = 0.01 was set"*)

This means the CryptBudding value `10.0` is **20× the Chaste source default** of 0.5, and
roughly **2000× the original NH paper value**.

**Literature reference:** Farhadifar et al. (2007). Scaled to K_A = 100, the soft-regime
value is Λ = **12.0**.

**Assessment:** The current `γ_cc = 10.0` is close to the Farhadifar soft-regime value of 12.0.
This is coincidentally reasonable despite the motivation being unclear from the JSON alone.

**Suggested value:** `nhCellCellAdhesion = 12.0`
(Matches Farhadifar soft regime with K_A = 100.)

---

### γ_boundary — Cell–boundary line tension (`nhBoundaryAdhesion`)

**Current JSON value:** `12.0`
**Chaste source default:** `1.0`

**Assessment:** Set slightly above γ_cc (12 vs 10) following Chaste convention (Fletcher
et al., 2014, *Biophysical Journal* 106, 2291–2304. doi:[10.1016/j.bpj.2013.11.4498](https://doi.org/10.1016/j.bpj.2013.11.4498)).
The boundary is the free tissue edge where cells are less adhesive to each other.

**Suggested value:** `nhBoundaryAdhesion = 14.0`
(≈1.2 × γ_cc, maintaining the same ratio if γ_cc changes to 12.)

---

### Differential adhesion — per-cell-type γ matrix

**Current JSON values:**
- Stem–stem: 5.0 (lower → more adhesive → stems cluster at crypt base)
- Stem–transit, stem–diff: 5.0
- Transit–transit, transit–diff, diff–diff: 10.0

**No directly measured Nagai-Honda γ values exist for individual intestinal cell types.**
The values are motivated by converging experimental and computational evidence:

1. **Differential Adhesion Hypothesis:** Steinberg (1963), *Science* 141, 401–408.
   doi:[10.1126/science.141.3579.401](https://doi.org/10.1126/science.141.3579.401)

2. **EphrinB/EphB signalling gradient:** Batlle et al. (2002), *Cell* 111, 251–263.
   doi:[10.1016/S0092-8674(02)01015-2](https://doi.org/10.1016/S0092-8674(02)01015-2) —
   Lgr5⁺ stem cells express high EphB2/EphB3; the receptor gradient drives cell sorting.

3. **Cellular Potts Model of intestinal crypt** (closest to measured differential adhesion):
   van Leeuwen et al. (2009), *Cell Proliferation* 42, 617–636.
   doi:[10.1111/j.1365-2184.2009.00627.x](https://doi.org/10.1111/j.1365-2184.2009.00627.x)
   uses a J-matrix fitted to EphB expression levels:
   - J(stem–stem) = **2** (highest adhesion, lowest J)
   - J(early TA–early TA) ≈ 4–5 (moderate)
   - J(late TA–late TA) ≈ 12–15 (lowest adhesion)
   This gives a stem:differentiated adhesion ratio of **2:12 ≈ 0.17** — stems are ~6× more
   adhesive than late TA/differentiated cells.

4. **3D intestinal organoid vertex model** (most relevant experiment):
   Pérez-González et al. (2021), *Nature Cell Biology* 23, 745–757.
   doi:[10.1038/s41556-021-00699-6](https://doi.org/10.1038/s41556-021-00699-6) —
   Used TFM + laser nanosurgery on mouse intestinal organoids to measure regional tension
   differences. Found the crypt base (stem/Paneth) has **high apical contractility** driving
   inward folding, while the TA zone has **high basal contractility** driving outward expansion.
   No per-cell-type adhesion γ values reported separately.

**Assessment of current values:** The ratio nhStemStemAdhesion/nhTransitTransitAdhesion = 0.5
is more conservative than the CPM literature suggests (ratio ≈ 0.17 from EphB expression).
Stems may need to be significantly more adhesive to reproduce robust crypt-base localisation.

**Suggested values** (if γ_cc = 12.0 baseline):
- Stem–stem: `2.0–4.0` (ratio 0.17–0.33 relative to baseline, consistent with CPM literature)
- Transit–transit, Diff–Diff: `12.0` (baseline)
- Stem–transit, stem–diff: geometric mean ≈ `5.0–7.0`
- Boundary adhesion per type: scale proportionally from nhBoundaryAdhesion

**Note:** Without laser ablation measurements for each individual cell type in intestinal
crypts (as done for Drosophila by Farhadifar), these remain informed estimates.

---

## 2. Surface Tension Force (vertex3d: OrganoidChaste `SurfaceTensionForce`)

Each cell has separate surface tensions for apical (γ_a), basal (γ_b), and lateral (γ_l) faces.
The SurfaceTensionForce source default (line 48–50 of `SurfaceTensionForce.cpp`) is 0.0 for
all tensions — all values must be set explicitly.

### gammaApical, gammaBasal, gammaLateral

**Current CryptBudding values:** 0.85, 0.85, 0.70 (apical = basal; symmetric)

**OrganoidChaste published test values** (`TestSphereSimulation.hpp` line 77):
```
apical = 1.0, basal = 0.7, lateral = 0.7
```
These are from the codebase of Drozdowski & Schwarz (2024), *Physical Review Research* 6,
L022045. doi:[10.1103/PhysRevResearch.6.L022045](https://doi.org/10.1103/PhysRevResearch.6.L022045)
(the paper introducing and validating OrganoidChaste's 3D vertex model).

**Experimental calibration from intestinal organoids:**
Pérez-González et al. (2021), *Nature Cell Biology* 23, 745–757.
doi:[10.1038/s41556-021-00699-6](https://doi.org/10.1038/s41556-021-00699-6) — used laser
nanosurgery and 3D traction force microscopy on mouse intestinal organoids to show:
- **Crypt base** (stem cell zone): apical contractility dominant → inward bending
- **TA zone:** basal contractility dominant → outward expansion and elongation

This is the closest experimental calibration of apical vs. basal tensions for intestinal crypt
cells, and explicitly contradicts using equal γ_apical = γ_basal across all cell types.

**Drosophila comparison (for magnitude):** Sui et al. (2018), *Nature Communications* 9, 4620.
doi:[10.1038/s41467-018-06497-3](https://doi.org/10.1038/s41467-018-06497-3) — laser ablation
of Drosophila wing disc gives **basal edge tension = 4× apical edge tension** from measured
recoil velocities. This ratio (1:4 apical:basal) is measured directly, not fitted.

**Assessment of CryptBudding values:** `gammaApical = gammaBasal = 0.85` is inconsistent with
the experimental evidence that apical ≠ basal in epithelial cells and in intestinal crypts in
particular. The symmetric assumption may suppress crypt folding dynamics.

**Comparison:**
| Parameter | CryptBudding | D&S (2024) | Experimental basis |
|-----------|-------------|------------|-------------------|
| gammaApical | 0.85 | 1.0 | Apical actomyosin ring dominates |
| gammaBasal | 0.85 | 0.7 | Pérez-González 2021; Sui 2018 (basal > apical in TA zone) |
| gammaLateral | 0.70 | 0.7 | Both consistent |

**Note on the TA zone:** If the TA zone basal tension > apical, then per Pérez-González 2021
the transit-amplifying cells should have `gammaBasal > gammaApical`. This would require
**per-cell-type apical/basal values** rather than a single global ratio with a scalar scale
factor. The current `gammaStemScale` etc. applies the same multiplicative scale to all faces,
which cannot reproduce the crypt-base (apical dominant) vs TA-zone (basal dominant) asymmetry.

**Suggested uniform-baseline values:** Use D&S (2024) as the validated baseline:
- `gammaApical = 1.0`, `gammaBasal = 0.7`, `gammaLateral = 0.7`

For a more biologically accurate model, consider separate apical/basal scale factors for
stem vs TA cells, informed by Pérez-González 2021.

### Cell-type tension scales

**Current values:** gammaStemScale=0.7, gammaTransitScale=1.0, gammaDiffScale=1.3

**No direct experimental source** maps measured cell stiffness to vertex model γ scales for
intestinal cell types. The qualitative ordering (stem < transit < differentiated) is motivated by:
- AFM studies showing Lgr5⁺ stem cells are mechanically softer than mature enterocytes
  (reviewed in Discher et al., 2005, *Science* 310, 1139–1143.
  doi:[10.1126/science.1116995](https://doi.org/10.1126/science.1116995))
- Paneth cells (represented by gammaDiffScale) are stiffer secretory cells

The ±30% modulation (0.7/1.0/1.3) is **not experimentally calibrated**. It is a conventional
choice that produces visible but not extreme cell-type segregation.

**Suggested values:** These could remain as-is (0.7/1.0/1.3) but should be noted as
uncalibrated. A sensitivity analysis over this range would be appropriate for the thesis.

---

## 3. Viscoelastic ECM — SLS Parameters

The ECM is modelled as a Standard Linear Solid:

$$E(t) = E_0 + E_1\,e^{-t/\tau}$$

### E₀ — Relaxed stiffness (`ghostE0 = 5.0`)

**Source:** Fertala et al. (2025), bioRxiv doi:[10.1101/2025.07.02.662292](https://doi.org/10.1101/2025.07.02.662292)
*(preprint; not yet peer-reviewed)*

AFM stress-relaxation experiments on starPEG/heparin IPN hydrogel microgels (diameter
67.8 ± 2.4 µm, comparable to cell-diameter scale) fitted with a Neo-Hookean + single-term
Prony FE model (Figure 2E):

**Measured microgel values:** E₀ ≈ 0.5–2.0 kPa, median **~800 Pa**

Microgel-scale values are used in preference to bulk gel values because bulk relaxation is
dominated by poroelastic fluid migration (absent at cell-diameter scale), inflating apparent
E₀ and τ. The Chaste value `5.0` (dimensionless force/CD) cannot be directly converted to
Pa without knowing η — only the **ratios** E₁/E₀ and the relaxation time τ are
physically meaningful comparisons.

### E₁ — Relaxation modulus (`ghostE1`)

**Current JSON value:** `2.0` → E₁/E₀ = **0.40**
**C++ source default** (`ViscoelasticGhostNodeEcmField.hpp` line 165): `1.0` → E₁/E₀ = **0.20**

**Measured from Fertala et al. (2025)** (Figure 2F and parameter sweep): E₁/E₀ = **0.05–0.20**

**Assessment:** The C++ source default ratio (0.20) sits at the upper measured bound — borderline
justified. The JSON override of 2.0 (ratio 0.40) **exceeds the experimentally measured range**
by a factor of 2. This was presumably chosen to make viscoelastic effects more observable
during the simulation, but does not represent the measured hydrogel.

**Suggested value:** `ghostE1 = 1.0` (E₁/E₀ = 0.20 — upper bound of measured range)

If you want to model softer viscoelasticity closer to the median, use `0.5` (E₁/E₀ = 0.10).

### τ — Relaxation time (`ghostRelaxationTime = 1.0` h)

**Measured from Fertala et al. (2025)** (Figure 2G): τ = **0.5–1.0 s** for microgels.

In simulation time (hours): 1 s ≈ 2.8 × 10⁻⁴ h. The current `τ = 1.0 h` is therefore
**~3600× longer** than the AFM-measured molecular relaxation time.

**This mismatch is deliberate but should be acknowledged:** The AFM measurement captures
rapid physical network rearrangements (non-covalent electrostatic bonds in the IPN hydrogel
relaxing in seconds). The simulation τ is intended to represent ECM-scale remodelling
(MMP-driven degradation, fibre rearrangement) which occurs on hour timescales in vivo.
No direct measurement of this longer-timescale remodelling τ exists for intestinal
organoid hydrogels.

**Suggested value:** `ghostRelaxationTime = 1.0` h can be kept but this rationale should be
stated clearly. A parameter sweep over τ = 0.5–5.0 h would quantify sensitivity.

### ghostDamping (`5.0`) and ghostFibreRemodelingRate (`0.1`), ghostAnisotropyStrength (`0.5`)

**No experimental source.** These are equation-of-motion and phenomenological parameters
with no measured equivalents. They control numerical equilibration rate and fibre alignment
speed, not the constitutive response. Current values are reasonable starting points for
stability; sensitivity analysis is appropriate.

---

## 4. Spring Force (`springStiffness = 15.0`)

**Source:** Meineke et al. (2001), *Cell Proliferation* 34, 253–266.
doi:[10.1046/j.0960-7722.2001.00216.x](https://doi.org/10.1046/j.0960-7722.2001.00216.x)

The Meineke spring model does not provide a directly measured µ; it is fitted to reproduce
observed crypt geometry. Sandersius & Newman (2008), *Phys. Biol.* 5, 015002
doi:[10.1088/1478-3975/5/1/015002](https://doi.org/10.1088/1478-3975/5/1/015002) provide
the most systematic analysis of off-lattice spring parameters.

**Assessment:** µ = 15 force/CD is within the range (10–30) used in published intestinal
crypt simulations but lacks a single direct experimental calibration.

---

## Summary Table

| Parameter | Current JSON | Suggested | Source | Status |
|-----------|-------------|-----------|--------|--------|
| λ_d (area elasticity) | 100 (hardcoded) | 100 | Nagai & Honda (2001) Chaste convention | ✓ Keep |
| nhMembraneSurface (λ_m) | 10.0 | **4.0** | Farhadifar et al. (2007) soft regime, K_A=100 | ⚠ Too high by 2.5× |
| nhCellCellAdhesion (γ_cc) | 10.0 | **12.0** | Farhadifar et al. (2007) soft regime | ✓ Close |
| nhBoundaryAdhesion (γ_b) | 12.0 | **14.0** | Fletcher et al. (2014); scale with γ_cc | ✓ Close |
| Stem–stem adhesion ratio | 0.5 × γ_cc | **0.17–0.33 × γ_cc** | van Leeuwen CPM (2009); J=2 vs J=12-15 | ✗ Too weak; stems not adhesive enough |
| gammaApical | 0.85 | **1.0** | Drozdowski & Schwarz (2024); apical actomyosin | ⚠ Review |
| gammaBasal | 0.85 | **0.7** | Pérez-González (2021); D&S (2024); apical ≠ basal | ✗ Should not equal apical |
| gammaLateral | 0.70 | 0.7 | Drozdowski & Schwarz (2024) | ✓ OK |
| gammaStemScale (uniform) | 0.7 | — | Cannot reproduce crypt-base vs TA-zone tension asymmetry | ✗ Architecture limitation |
| ghostE0 (E₀) | 5.0 | 5.0 | Fertala et al. (2025); ratio only | ✓ OK |
| ghostE1 (E₁) | 2.0 | **1.0** | Fertala et al. (2025); current exceeds range | ✗ Reduce |
| ghostRelaxationTime (τ) | 1.0 h | 1.0 h | Fertala 0.5–1 s (molecule); ~hour for ECM remodelling | ⚠ Acknowledge |
| ghostDamping | 5.0 | 5.0 | No source; numerical | ✗ No source |
| ghostFibreRemodelingRate | 0.1 | 0.1 | No source; phenomenological | ✗ No source |
| springStiffness | 15.0 | 15.0 | Meineke (2001); Sandersius (2008) | ⚠ No single calibration |

**Legend:** ✓ = justified; ⚠ = partial/qualitative justification; ✗ = no experimental source.

---

## Key References

1. Nagai & Honda (2001). A dynamic cell model for the formation of epithelial tissues.
   *Philosophical Magazine B* 81(7), 699–719. doi:[10.1080/13642810108205772](https://doi.org/10.1080/13642810108205772)

2. Farhadifar, Röper, Aigouy, Eaton & Jülicher (2007). The influence of cell mechanics,
   cell-cell interactions, and proliferation on epithelial packing. *Current Biology* 17(24),
   2095–2104. doi:[10.1016/j.cub.2007.11.049](https://doi.org/10.1016/j.cub.2007.11.049)

3. Drozdowski & Schwarz (2024). [organoid vertex model paper]. *Physical Review Research* 6,
   L022045. doi:[10.1103/PhysRevResearch.6.L022045](https://doi.org/10.1103/PhysRevResearch.6.L022045)

4. Fertala, Uhlmann, Grigoryev, Seth, Friedrichs, Werner & Balzani (2025).
   Scale-Specific Viscoelastic Characterization of Hydrogels: Integrated AFM and Finite
   Element Modeling. *bioRxiv* doi:[10.1101/2025.07.02.662292](https://doi.org/10.1101/2025.07.02.662292)
   *(preprint)*

5. Steinberg (1963). Reconstruction of tissues by dissociated cells. *Science* 141, 401–408.
   doi:[10.1126/science.141.3579.401](https://doi.org/10.1126/science.141.3579.401)

6. Batlle et al. (2002). β-Catenin and TCF mediate cell positioning in the intestinal
   epithelium by controlling the expression of EphB/ephrinB. *Cell* 111, 251–263.
   doi:[10.1016/S0092-8674(02)01015-2](https://doi.org/10.1016/S0092-8674(02)01015-2)

7. van Leeuwen et al. (2009). An integrative computational model for intestinal tissue renewal.
   *Cell Proliferation* 42, 617–636. doi:[10.1111/j.1365-2184.2009.00627.x](https://doi.org/10.1111/j.1365-2184.2009.00627.x)

8. Pérez-González et al. (2021). Mechanical compartmentalization of the intestinal organoid
   enables crypt folding and collective cell migration. *Nature Cell Biology* 23, 745–757.
   doi:[10.1038/s41556-021-00699-6](https://doi.org/10.1038/s41556-021-00699-6)
   *(Key: TFM + laser nanosurgery on mouse intestinal organoids; regional apical/basal tensions)*

9. Yang et al. (2021). Cell fate coordinates mechano-osmotic forces in intestinal crypt
   formation. *Nature Cell Biology* 23, 733–744. doi:[10.1038/s41556-021-00700-2](https://doi.org/10.1038/s41556-021-00700-2)
   *(3D closed-monolayer vertex model; differential crypt/villus tension σ = 0.02)*

10. Sui et al. (2018). Differential lateral and basal tension drive folding of Drosophila
    wing discs through two distinct mechanisms. *Nature Communications* 9, 4620.
    doi:[10.1038/s41467-018-06497-3](https://doi.org/10.1038/s41467-018-06497-3)
    *(Laser ablation: basal edge tension = 4× apical edge tension, directly measured)*

11. Ratheesh et al. (2012). Centralspindlin and α-catenin regulate Rho signalling at the
    epithelial zonula adherens. *Nature Cell Biology* 14, 67–76. doi:[10.1038/ncb2390](https://doi.org/10.1038/ncb2390)

9. Fletcher, Osborne, Maini & Gavaghan (2014). Implementing vertex dynamics models of cell
   monolayers. *Biophysical Journal* 106, 2291–2304. doi:[10.1016/j.bpj.2013.11.4498](https://doi.org/10.1016/j.bpj.2013.11.4498)

10. Meineke, Potten & Loeffler (2001). Cell migration and organization in the intestinal
    crypt using a lattice-free model. *Cell Proliferation* 34, 253–266.
    doi:[10.1046/j.0960-7722.2001.00216.x](https://doi.org/10.1046/j.0960-7722.2001.00216.x)

11. Sandersius & Newman (2008). Modeling cell rheology with the Subcellular Element Model.
    *Physical Biology* 5, 015002. doi:[10.1088/1478-3975/5/1/015002](https://doi.org/10.1088/1478-3975/5/1/015002)

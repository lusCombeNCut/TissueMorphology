# Viscoelastic Ghost Node ECM — Methodology

## Overview

This document describes the viscoelastic constitutive model used in `ViscoelasticGhostNodeEcmField` and `ViscoelasticGhostNodeEcmForce`, which replaces the purely elastic spring model in the original `GhostNodeEcmField`. The formulation is based on the **generalized Maxwell model** (single-term Prony series) from Fertala et al. (2025), adapted for a discrete particle-based ECM representation at the cell-diameter length scale.

> **Reference**: Fertala N, Uhlmann K, Grigoryev E, et al. *Scale-Specific Viscoelastic Characterization of Hydrogels: Integrated AFM and Finite Element Modeling.* bioRxiv (2025). doi:10.1101/2025.07.02.662292

---

## 1. Problem with the Previous Implementation

The original `GhostNodeEcmField` used:

$$\eta_{\text{drag}} \, \dot{\mathbf{x}} = \mathbf{F}_{\text{spring}}(d_{ij} - s_0)$$

where $s_0$ is a **fixed, global** rest length. The drag coefficient $\eta_{\text{drag}}$ governs how quickly the system equilibrates but is **not** a constitutive material parameter — it comes from the overdamped equation of motion (low-Reynolds-number regime), not from the material's stress–strain relationship.

This system is **purely elastic**: at any constant deformation, the stress is constant (no relaxation). Remove the deformation, and the system returns to exactly its original state. The damping only controls the *rate* of return, not the *target*.

A viscoelastic material, by contrast, exhibits:
- **Stress relaxation**: stress decays at constant strain
- **Creep**: strain increases at constant stress
- **Rate-dependent stiffness**: faster loading → stiffer response

None of these are present with fixed rest lengths and equation-of-motion drag alone.

---

## 2. Continuum Constitutive Model

The generalized Maxwell model (or Standard Linear Solid with a single relaxation mode) defines a time-dependent Young's modulus:

$$E(t) = E_0 + E_1 \, e^{-t/\tau}$$

| Symbol | Meaning | Role |
|--------|---------|------|
| $E_0$ | Relaxed (equilibrium) modulus | Long-time elastic stiffness after all viscous dissipation |
| $E_1$ | Relaxation modulus | Viscous contribution that decays over time |
| $E_u = E_0 + E_1$ | Instantaneous (unrelaxed) modulus | Stiffness at $t = 0$ (immediate response) |
| $\tau = \eta_{\text{VE}} / E_1$ | Relaxation time | Timescale for stress relaxation |
| $\eta_{\text{VE}}$ | Dashpot viscosity | **Constitutive** viscosity (not EoM drag) |

### Rheological analogue

This is a spring ($E_0$) in parallel with a Maxwell element (spring $E_1$ + dashpot $\eta_{\text{VE}}$ in series):

```
         ┌─── E_0 (spring) ───┐
    ─────┤                     ├───── 
         └── E_1 ──┤├── η_VE ──┘
              spring   dashpot
              (series: Maxwell arm)
```

**Behaviour**:
- At $t = 0$: both springs resist → total stiffness $E_u = E_0 + E_1$
- At $t \gg \tau$: dashpot has fully yielded → only $E_0$ remains
- Stress relaxation at constant strain: $\sigma(t) = E_0 \varepsilon + E_1 \varepsilon \, e^{-t/\tau}$

---

## 3. Discrete Spring Analogue (Per-Pair Rest Length Evolution)

For the particle-based ghost node system, we translate the continuum model into a discrete spring-and-dashpot system per connected pair $(i, j)$.

### Instantaneous spring force

$$\mathbf{F}_{ij} = k_u \cdot (d_{ij} - s_{ij}(t)) \cdot A_{ij} \cdot \hat{\mathbf{d}}_{ij}$$

where:
- $k_u = E_0 + E_1$ is the instantaneous spring stiffness
- $d_{ij} = \|\mathbf{x}_j - \mathbf{x}_i\|$ is the current pair distance
- $s_{ij}(t)$ is the **per-pair evolving rest length**
- $A_{ij}$ is the anisotropic modulation factor
- $\hat{\mathbf{d}}_{ij}$ is the unit vector from $i$ to $j$

### Rest length evolution (Maxwell dashpot)

$$\frac{ds_{ij}}{dt} = \frac{d_{ij} - s_{ij}}{\tau}$$

This is the key equation that produces viscoelasticity. The rest length $s_{ij}$ relaxes toward the current distance $d_{ij}$ on timescale $\tau$.

**Forward Euler discretisation:**

$$s_{ij}^{n+1} = s_{ij}^n + \frac{(d_{ij}^n - s_{ij}^n) \, \Delta t}{\tau}$$

### Why this works

Consider a pair suddenly stretched from equilibrium ($d = s_0$) to a new distance $d = s_0 + \delta$:

1. **At $t = 0$**: $s_{ij} = s_0$, force $= k_u \cdot \delta$ (full instantaneous stiffness)
2. **At $t \sim \tau$**: $s_{ij}$ evolves toward $s_0 + \delta$, force decays
3. **At $t \gg \tau$**: $s_{ij} \to s_0 + \delta$, force $\to 0$

But the equilibrium arm ($E_0$) is implicit: forcing $s_{ij}$ to evolve completely to $d_{ij}$ gives a pure Maxwell element (fluid-like at long times). To recover SLS behaviour (solid-like, with equilibrium stiffness $E_0$), the effective long-time force is non-zero because:

$$\text{Effective long-time force} = E_0 \cdot (d_{ij} - s_{ij,\text{initial}})$$

In practice, the rest length clamp bounds ($0.1 \times s_{\text{initial}}$ to $5 \times s_{\text{initial}}$) and the competing forces from all neighbours prevent complete relaxation of any single pair, preserving geometric integrity.

For a true SLS where the equilibrium modulus is guaranteed, one can split the force:

$$\mathbf{F}_{ij} = \underbrace{E_0 \cdot (d_{ij} - s_0)}_{\text{equilibrium spring}} + \underbrace{E_1 \cdot (d_{ij} - s_{ij}(t))}_{\text{Maxwell arm}}$$

This variant can be enabled if the pure Maxwell element proves too fluid-like at long times by setting `mRelaxedStiffness > 0` and treating it as a separate equilibrium spring. The current implementation uses the combined approach with clamped rest lengths for simplicity.

### Anisotropic modulation

The force magnitude is modulated by fibre direction:

$$A_{ij} = 1 - a \cdot \bar{\alpha}_{ij} \cdot (1 - |\cos\theta|)$$

where:
- $a$ = anisotropy strength parameter
- $\bar{\alpha}_{ij}$ = average anisotropy of nodes $i$ and $j$
- $\theta$ = angle between the displacement vector and the average fibre direction

### Position update (equation of motion)

The position update remains overdamped:

$$\mathbf{x}^{n+1} = \mathbf{x}^n + \frac{\mathbf{F}_{\text{total}}}{\eta_{\text{drag}}} \cdot \Delta t$$

**Critically, $\eta_{\text{drag}}$ is distinct from the constitutive viscosity $\eta_{\text{VE}} = E_1 \cdot \tau$.**

| Parameter | Role | Controls |
|-----------|------|----------|
| $\eta_{\text{drag}}$ | Equation-of-motion drag | How quickly nodes move to equilibrium |
| $\eta_{\text{VE}} = E_1 \cdot \tau$ | Constitutive dashpot viscosity | Whether the material relaxes stress |

---

## 4. Cell–Ghost Node Forces

Cell–ghost spring forces use the same Chaste spring regime as the original implementation. These are **not** viscoelastic — the cell–ECM interaction is an adhesion/contact force, not a constitutive material property:

**Compression** ($d < s_0$):
$$F = k_{cg} \cdot s_0 \cdot \ln\left(1 + \frac{d - s_0}{s_0}\right)$$

**Extension** ($d > s_0$):
$$F = k_{cg} \cdot \rho_j \cdot (d - s_0) \cdot e^{-5(d - s_0) / s_0}$$

where $\rho_j \in [0,1]$ is the local ECM density, giving density-weighted adhesion.

---

## 5. Parameter Values

### Viscoelastic constitutive parameters (from Fertala et al. 2025)

The paper measured IPN hydrogel microgels (diameter ~68 µm) using AFM stress relaxation with a 10 µm spherical probe. These measurements are at the cell-diameter length scale and are appropriate for our ghost node ECM.

| Parameter | Symbol | Value range (paper) | Default (Chaste) | Notes |
|-----------|--------|---------------------|-------------------|-------|
| Relaxed modulus | $E_0$ | 500–2000 Pa | 5.0 | Chaste force units; sets long-time stiffness |
| Relaxation modulus | $E_1$ | 5–20% of $E_0$ | 1.0 | Viscous contribution; $E_1/E_0 = 0.05\text{–}0.2$ |
| Relaxation time | $\tau$ | 0.5–2.0 s | 1.0 | Timescale for stress relaxation |
| Instantaneous modulus | $E_u = E_0 + E_1$ | 525–2400 Pa | 6.0 | Full stiffness at $t=0$ |
| Constitutive viscosity | $\eta_{\text{VE}} = E_1 \cdot \tau$ | — | 1.0 | Derived, not set independently |

### Mapping paper values to Chaste force units

Chaste off-lattice simulations use a non-dimensionalised force scale where typical Chaste spring stiffnesses are $\mathcal{O}(1\text{–}30)$, not in Pa. The paper's **ratios** are the key transferable quantities:

- **$E_1 / E_0 \approx 0.1\text{–}0.2$**: the viscous contribution is 10–20% of the equilibrium stiffness
- **$\tau \approx 0.5\text{–}2.0$ s**: directly usable if simulation time is in seconds (as in Chaste)
- **$\nu = 0.45$**: nearly incompressible (relevant for continuum FE, implicit in our 3D connectivity)

### Equation-of-motion parameters (not constitutive)

| Parameter | Symbol | Default | Notes |
|-----------|--------|---------|-------|
| EoM drag | $\eta_{\text{drag}}$ | 1.0 | Controls equilibration speed, not material response |
| Cell–ghost stiffness | $k_{cg}$ | 10.0 | Cell–ECM adhesion/contact spring |
| Cell–ghost rest length | $s_0^{cg}$ | 0.0 | Cell–ECM equilibrium distance |
| Cell–ghost cutoff | $d_{\text{cut}}^{cg}$ | 1.5 | Interaction radius for cell–ECM |

### Biological parameters

| Parameter | Symbol | Default | Notes |
|-----------|--------|---------|-------|
| Degradation rate | $k_{\text{deg}}$ | 0.02 | MMP secretion rate |
| Removal threshold | $\rho_{\text{min}}$ | 0.01 | Density below which nodes are removed |
| Fibre remodeling rate | $k_{\text{remodel}}$ | 0.1 | Rate of fibre alignment with traction |
| Anisotropy strength | $a$ | 0.5 | Fibre-direction modulation of stiffness |

---

## 6. Suitability of Paper Parameters for Cell-Diameter Scale

The Fertala et al. parameters are well-suited because:

1. **Length scale match**: Their microgels (~68 µm diameter) are probed with a 10 µm sphere, comparable to cell–ECM contacts at the ~10–20 µm cell diameter scale.

2. **Microgels vs bulk**: The paper demonstrates that microgels show **rapid, localised relaxation** with minimal poroelastic contribution — matching our discrete ghost node model which does not include fluid transport. Bulk gel parameters would be inappropriate as they include poroelastic effects.

3. **Single relaxation mode sufficient**: At the cell scale, the paper shows a single Prony term captures the relaxation profile with <5% error. Multi-term Prony series or fractional models are unnecessary.

4. **$E_1/E_0$ ratio**: The 10–20% ratio means the ECM is predominantly elastic with moderate stress relaxation — consistent with biophysical observations of cell-scale ECM behaviour where the matrix provides structural support but permits slow remodeling.

5. **$\tau$ in the seconds range**: Relaxation times of 0.5–2.0 s are biologically meaningful — fast enough to respond to cell migration (minutes) and division (hours), but slow enough that cells experience significant resistance during rapid events.

### Limitations

- The paper uses a **Neo-Hookean** hyperelastic model for the elastic component. Our discrete springs are linear, which is adequate for small-to-moderate strains but does not capture strain-stiffening at large deformations.
- **Poroelastic effects** (fluid migration through the ECM pore network) are not modelled. The paper shows these are negligible at the microgel scale but significant for bulk gels.
- The paper's constitutive parameters assume a **homogeneous, isotropic** material. Our anisotropic fibre modulation extends beyond their framework.

---

## 7. Output and Visualisation

The VTP output includes a new field:

- **`rest_length_strain`**: per-node average of $(s_{ij} - s_{\text{initial}}) / s_{\text{initial}}$ across all connected pairs. Positive values indicate rest lengths that have grown (material has relaxed / yielded under sustained deformation). This allows visualisation of where the ECM has undergone viscous remodeling.

---

## 8. Summary of Equations

| Equation | Expression |
|----------|------------|
| Continuum modulus | $E(t) = E_0 + E_1 \, e^{-t/\tau}$ |
| Spring force | $\mathbf{F}_{ij} = k_u (d_{ij} - s_{ij}) \, A_{ij} \, \hat{\mathbf{d}}_{ij}$ |
| Rest length ODE | $\dot{s}_{ij} = (d_{ij} - s_{ij}) / \tau$ |
| Rest length update | $s_{ij}^{n+1} = s_{ij}^n + (d_{ij}^n - s_{ij}^n) \Delta t / \tau$ |
| Anisotropy factor | $A_{ij} = 1 - a \bar{\alpha} (1 - |\cos\theta|)$ |
| Position update (EoM) | $\mathbf{x}^{n+1} = \mathbf{x}^n + \mathbf{F} \Delta t / \eta_{\text{drag}}$ |

---

## Files

| File | Description |
|------|-------------|
| `ViscoelasticGhostNodeEcmField.hpp` | ECM field with per-pair rest lengths and Maxwell rest length evolution |
| `ViscoelasticGhostNodeEcmForce.hpp` | Force class orchestrating cell–ECM forces, degradation, remodeling, and viscoelastic update |

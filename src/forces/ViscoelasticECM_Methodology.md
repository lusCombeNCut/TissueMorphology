# Viscoelastic Ghost Node ECM — Methodology

## Overview

This document describes the constitutive model and implementation of `ViscoelasticGhostNodeEcmField` and `ViscoelasticGhostNodeEcmForce`. These classes replace the purely elastic spring model of the original `GhostNodeEcmField` with a Standard Linear Solid (SLS/Zener model) constitutive law, targeting the single-relaxation-mode viscoelastic framework of Fertala et al. (2025).

> **Reference**: Fertala N, Uhlmann K, Grigoryev E, et al. *Scale-Specific Viscoelastic Characterization of Hydrogels: Integrated AFM and Finite Element Modeling.* bioRxiv (2025). doi:10.1101/2025.07.02.662292

---

## 1. Continuum Constitutive Model

The Fertala et al. (2025) FE model combines a Neo-Hookean backbone (elastic) with a single-term Prony series (viscous) in ABAQUS. In one dimension this gives the SLS relaxation modulus:

$$E(t) = E_0 + E_1 \, e^{-t/\tau}$$

| Symbol | Meaning | Role |
|--------|---------|------|
| $E_0$ | Relaxed (equilibrium) modulus | Long-time elastic stiffness; non-zero at $t \to \infty$ |
| $E_1$ | Relaxation modulus | Viscous contribution decaying on timescale $\tau$ |
| $E_u = E_0 + E_1$ | Instantaneous (unrelaxed) modulus | Full stiffness at $t = 0$ |
| $\tau = \eta_{\text{VE}} / E_1$ | Relaxation time | Timescale for stress relaxation |
| $\eta_{\text{VE}}$ | Constitutive dashpot viscosity | Distinct from equation-of-motion drag |

The rheological analogue is a spring $E_0$ in **parallel** with a Maxwell element (spring $E_1$ + dashpot $\eta_{\text{VE}}$ in series):

```
     ┌─── E_0 (spring) ───┐
─────┤                     ├─────
     └── E_1 ──┤├── η_VE ──┘
          spring   dashpot
```

At $t \to \infty$, the dashpot has fully yielded but the parallel spring $E_0$ remains, so the long-time stress is non-zero. This is the key difference from a pure Maxwell fluid (no $E_0$ arm), where stress decays completely to zero.

---

## 2. Discrete SLS Force Decomposition

Each connected ghost-node pair $(i,j)$ carries two rest lengths:

- $s_0^{ij}$ — **permanent rest length**: set to the initial pair distance during construction, never updated; stored in `initial_rest_lengths`. Represents the covalent crosslink network of the IPN hydrogel.
- $s_{ij}(t)$ — **evolving dashpot length**: relaxes toward the current pair distance on timescale $\tau$; stored in `rest_lengths`. Represents reversible physical network rearrangements.

The force on pair $(i,j)$ is:

$$\mathbf{F}_{ij} = \left[\, E_0 \bigl(d_{ij} - s_0^{ij}\bigr) + E_1 \bigl(d_{ij} - s_{ij}(t)\bigr) \,\right] A_{ij} \, \hat{\mathbf{r}}_{ij}$$

where $d_{ij} = \|\mathbf{x}_j - \mathbf{x}_i\|$ is the current separation and $A_{ij}$ is the anisotropic modulation factor (see Section 4).

**Verification**:
- At $t = 0$: $s_{ij} = s_0^{ij}$, giving $F = E_u \cdot (d_{ij} - s_0^{ij}) \cdot A_{ij}$ ✓
- As $t \to \infty$: $s_{ij} \to d_{ij}$, giving $F = E_0 \cdot (d_{ij} - s_0^{ij}) \cdot A_{ij} \neq 0$ ✓ (solid-like)

---

## 3. Rest Length Evolution

The evolving dashpot length obeys:

$$\frac{ds_{ij}}{dt} = \frac{d_{ij} - s_{ij}}{\tau}$$

This ODE is integrated using the **exact exponential integrator**:

$$s_{ij}^{n+1} = d_{ij}^n + \bigl(s_{ij}^n - d_{ij}^n\bigr) \, e^{-\Delta t / \tau}$$

This solves $\dot{u} = -u/\tau$ (where $u = s_{ij} - d_{ij}$) analytically over each timestep, and is unconditionally stable for all $\Delta t > 0$ and $\tau > 0$.

---

## 4. Anisotropic Modulation

The force magnitude is modulated by fibre direction:

$$A_{ij} = 1 - a \cdot \bar{\alpha}_{ij} \cdot (1 - |\cos\theta|)$$

| Symbol | Meaning |
|--------|---------|
| $a$ | Anisotropy strength parameter |
| $\bar{\alpha}_{ij}$ | Mean fibre anisotropy of nodes $i$ and $j$ |
| $\theta$ | Angle between $\hat{\mathbf{r}}_{ij}$ and mean fibre direction |

Forces along the fibre direction ($\theta = 0$) are unmodified ($A = 1$); forces perpendicular to the fibre are reduced.

---

## 5. Ghost Node Position Update

Ghost node positions are updated using overdamped (inertia-free) dynamics:

$$\eta_{\text{drag}} \frac{d\mathbf{x}_j^{\text{g}}}{dt}
= \sum_{k \in \mathcal{G}(j)} \mathbf{F}_{jk}^{\text{ve}}
+ \sum_{i \in \mathcal{C}(j)} \mathbf{F}_{ji}^{\text{cg,react}}$$

where $\eta_{\text{drag}}$ is the viscous drag coefficient and $\mathcal{C}(j)$ is the set of cells exerting reaction forces on ghost node $j$. Discretised as:

$$\mathbf{x}_j^{n+1} = \mathbf{x}_j^n + \frac{\mathbf{F}_j^{\text{total}}}{\eta_{\text{drag}}} \, \Delta t$$

**$\eta_{\text{drag}}$ is an equation-of-motion parameter, not a constitutive parameter.** It controls how quickly nodes reach mechanical equilibrium, not whether the material relaxes stress. The constitutive viscosity is $\eta_{\text{VE}} = E_1 \tau$.

---

## 6. Cell–Ghost Node Forces

Cell–ECM interaction forces are not viscoelastic — they represent adhesion and contact. For a cell at position $\mathbf{x}_c$ and ghost node at distance $d$ with density $\rho$:

**Compression** ($d < s_0$):
$$F = k_{\text{cg}} \cdot s_0 \cdot \ln\!\left(1 + \frac{d - s_0}{s_0}\right)$$

**Extension** ($d > s_0$):
$$F = k_{\text{cg}} \cdot \rho \cdot (d - s_0) \cdot e^{-5(d - s_0)/s_0}$$

where $\rho \in [0,1]$ is the local ECM density, giving density-weighted adhesion. Newton's 3rd law reaction is applied to the ghost node.

---

## 7. ECM Degradation and Fibre Remodeling

### Degradation (MMP secretion)

Each cell degrades nearby ghost nodes at rate $k_{\text{deg}}$:

$$\dot{\rho}_j = -k_{\text{deg}} \cdot w_{ij} \cdot \rho_j$$

where the weight $w_{ij} = \max(0, 1 - d_{ij} / d_{\text{cut}})$ is distance-weighted within the cutoff radius. Ghost nodes with $\rho < \rho_{\min}$ are deactivated and removed from the connectivity graph.

### Fibre remodeling (cell traction)

Fibre orientation evolves toward the cell traction direction at rate $k_{\text{remodel}}$:

$$\dot{\mathbf{f}}_j = k_{\text{remodel}} \cdot w_{ij} \cdot (\mathbf{t}_i - \mathbf{f}_j)$$

where $\mathbf{t}_i = (\mathbf{x}_i - \bar{\mathbf{x}}) / \|\mathbf{x}_i - \bar{\mathbf{x}}\|$ is the outward radial traction direction from the tissue centroid $\bar{\mathbf{x}}$.

---

## 8. Parameter Values

### Viscoelastic constitutive parameters (from Fertala et al. 2025)

Measured on IPN hydrogel microgels (diameter $67.8 \pm 2.4\,\mu$m) using AFM stress relaxation with a 10 µm colloidal probe sphere, inverse-fitted to a Neo-Hookean + single-term Prony ABAQUS model:

| Parameter | Symbol | Measured (microgels) | Default (Chaste) | Notes |
|-----------|--------|----------------------|------------------|-------|
| Relaxed modulus | $E_0$ | 600–900 Pa (median ~800 Pa) | 5.0 | Long-time stiffness |
| Relaxation modulus | $E_1$ | 60–150 Pa ($E_1/E_0 \approx 0.05$–$0.20$) | 1.0 | Viscous arm |
| Relaxation time | $\tau$ | 0.5–1.0 s (median ~0.6 s) | 1.0 | Stress relaxation timescale |
| Instantaneous modulus | $E_u = E_0 + E_1$ | 660–1050 Pa | 6.0 | Stiffness at $t=0$ |

Microgel-scale parameters are used rather than bulk gel values because bulk relaxation includes poroelastic fluid migration (absent at the cell-diameter scale) that inflates apparent $E_1$ and $\tau$.

The paper's Prony series uses the notation: $g_R(t) = 1 - g_1(1 - e^{-t/\tau_1})$ with $g_1 = E_1/E_u$. The correspondence is $E_0 = E_u(1 - g_1)$ and $E_1 = E_u \cdot g_1$.

### Equation-of-motion parameters

| Parameter | Symbol | Default | Notes |
|-----------|--------|---------|-------|
| EoM drag | $\eta_{\text{drag}}$ | 1.0 | Equilibration rate, not material property |
| Cell–ghost stiffness | $k_{\text{cg}}$ | 10.0 | Cell–ECM adhesion/contact spring |
| Cell–ghost rest length | $s_0^{\text{cg}}$ | 0.0 | Cell–ECM equilibrium distance |
| Cell–ghost cutoff | $d_{\text{cut}}^{\text{cg}}$ | 1.5 | Interaction radius |

### Biological parameters

| Parameter | Symbol | Default | Notes |
|-----------|--------|---------|-------|
| Degradation rate | $k_{\text{deg}}$ | 0.02 | MMP secretion rate |
| Removal threshold | $\rho_{\min}$ | 0.01 | Density below which nodes are removed |
| Fibre remodeling rate | $k_{\text{remodel}}$ | 0.1 | Rate of fibre alignment with traction |
| Anisotropy strength | $a$ | 0.5 | Fibre-direction modulation of stiffness |

---

## 9. VTP Output Fields

| Field | Description |
|-------|-------------|
| `density` | ECM density $\rho \in [0,1]$ |
| `fibre_direction` | Unit fibre direction vector |
| `anisotropy` | Degree of fibre alignment $\in [0,1]$ |
| `rest_length_strain` | Per-node mean $(s_{ij} - s_0^{ij})/s_0^{ij}$: positive where the dashpot arm has relaxed/yielded |
| `id` | Ghost node identifier |

---

## 10. Summary of Equations

| Equation | Expression |
|----------|------------|
| Continuum modulus | $E(t) = E_0 + E_1 \, e^{-t/\tau}$ |
| Ghost–ghost force | $\mathbf{F}_{ij} = [E_0(d_{ij} - s_0^{ij}) + E_1(d_{ij} - s_{ij})] \, A_{ij} \, \hat{\mathbf{r}}_{ij}$ |
| Rest length ODE | $\dot{s}_{ij} = (d_{ij} - s_{ij}) / \tau$ |
| Rest length update | $s_{ij}^{n+1} = d_{ij}^n + (s_{ij}^n - d_{ij}^n) \, e^{-\Delta t/\tau}$ |
| Anisotropy factor | $A_{ij} = 1 - a \bar{\alpha}_{ij} (1 - |\cos\theta|)$ |
| Position update (EoM) | $\mathbf{x}^{n+1} = \mathbf{x}^n + \mathbf{F}_{\text{total}} \Delta t / \eta_{\text{drag}}$ |

---

## Files

| File | Description |
|------|-------------|
| `ViscoelasticGhostNodeEcmField.hpp` | ECM field: ghost node data, SLS force computation, rest length evolution, degradation, remodeling |
| `ViscoelasticGhostNodeEcmForce.hpp` | Force class: cell–ghost forces, ECM dynamics orchestration, periodic node removal |

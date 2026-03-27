# Viscoelastic Ghost Node ECM — Development History

This document records the design decisions, errors found, and corrections made during development of the viscoelastic ECM implementation. It is a record of the *process*, not the current state. For a description of the current code, see `ViscoelasticECM_Methodology.md`.

---

## Stage 1: On-Lattice Elastic Model (GhostNodeEcmField)

The original implementation used an on-lattice ghost node grid with **fixed rest-length linear springs** for ghost–ghost interactions:

$$\mathbf{F}_{\text{spring}} = k \cdot (d_{ij} - s_0) \cdot \hat{\mathbf{r}}_{ij}$$

where $s_0$ is a single global rest length, fixed for all time. The drag coefficient $\eta_{\text{drag}}$ in the overdamped position update provides rate-of-return control but **no viscoelastic constitutive behaviour**:

$$\eta_{\text{drag}} \, \dot{\mathbf{x}} = \mathbf{F}_{\text{spring}}(d_{ij} - s_0)$$

**Limitations identified:**

- **Purely elastic**: at constant deformation the stress is constant — no relaxation, no creep, no rate-dependent stiffness.
- **No permanent remodeling**: if the organoid expanded and degraded surrounding ECM, the vacated region had no restorative force; it remained a permanent void.
- Biologically, hydrogels do exhibit stress relaxation at the cell-diameter scale (Fertala et al. 2025; Lai 2025).

This motivated the introduction of a viscoelastic constitutive law.

---

## Stage 2: First Viscoelastic Attempt — Single Evolving Rest Length (Maxwell Fluid)

### Design intent

To introduce stress relaxation, the rest length was made **evolving**: each ghost–ghost pair $(i,j)$ stores $s_{ij}(t)$, which relaxes toward the current distance on timescale $\tau$:

$$\dot{s}_{ij} = \frac{d_{ij} - s_{ij}}{\tau}$$

The spring force was:

$$\mathbf{F}_{ij} = k_u \cdot (d_{ij} - s_{ij}) \cdot A_{ij} \cdot \hat{\mathbf{r}}_{ij}$$

where $k_u = E_0 + E_1$ (the "instantaneous stiffness"). The intent was to use the constitutive model of Fertala et al. (2025):

$$E(t) = E_0 + E_1 \, e^{-t/\tau}$$

with $E_0$ providing a non-zero long-time stiffness.

### Bug: this is a Maxwell fluid, not an SLS

**Proof**: Let $x = d_{ij} - s_{ij}$. Under constant deformation ($\dot{d} = 0$):

$$\dot{x} = \dot{d} - \dot{s}_{ij} = 0 - \frac{x}{\tau} = -\frac{x}{\tau}$$

Solution: $x(t) = x_0 \, e^{-t/\tau}$. Therefore:

$$F(t) = k_u \cdot x_0 \, e^{-t/\tau} \xrightarrow{t \to \infty} 0$$

The force decays **completely to zero** at long times, regardless of the value of $E_0$. This is a **pure Maxwell fluid** — the material flows indefinitely under sustained loading. The parameter $E_0$ has no role in the long-time behaviour because both the elastic and viscous arms share the same evolving rest length.

This contradicts the Fertala et al. (2025) model, in which $E_0$ is the "fully relaxed Young's modulus" measured from the plateau of the stress relaxation curve — explicitly non-zero.

### Masking via rest-length clamps

The implementation included clamps on the rest length:

```cpp
double min_rest = mInitialSpacing * 0.1;
double max_rest = mInitialSpacing * 5.0;
new_s = std::max(min_rest, std::min(max_rest, new_s));
```

These prevented the rest length from fully reaching $d_{ij}$, which partially masked the Maxwell fluid behaviour by preventing complete relaxation. However, this is a numerical workaround, not a physical model: the clamp range is arbitrary and the residual stiffness depends on grid spacing rather than the constitutive parameter $E_0$.

### Additional bug: conditionally stable Forward Euler integrator

The rest length ODE was integrated with Forward Euler:

$$s_{ij}^{n+1} = s_{ij}^n + \frac{(d_{ij}^n - s_{ij}^n) \, \Delta t}{\tau}$$

This has amplification factor $(1 - \Delta t/\tau)$. Stability requires:

$$\left|1 - \frac{\Delta t}{\tau}\right| < 1 \implies \Delta t < 2\tau$$

For stiff gels (small $\tau$) or large timesteps, the rest length oscillates with growing amplitude. The rest-length clamps partially masked this by limiting the divergence, but the underlying integrator remained conditionally unstable.

---

## Stage 3: SLS Split-Force Fix + Exact Integrator (Current Implementation)

### Fix 1: SLS split-force decomposition

The root cause of the Maxwell fluid bug is that a single evolving rest length cannot represent an SLS: the permanent elastic arm $E_0$ and the relaxing Maxwell arm $E_1$ must evolve **independently**. The fix stores two rest lengths per pair:

- $s_0^{ij}$ — **permanent** rest length: set to the initial pair distance at construction, never changed. Stored in `initial_rest_lengths`. Represents the covalent crosslink network.
- $s_{ij}(t)$ — **evolving dashpot length**: relaxes toward the current distance on timescale $\tau$. Stored in `rest_lengths`. Represents reversible physical network rearrangements.

The force is decomposed into the two SLS arms:

$$\mathbf{F}_{ij} = \left[E_0(d_{ij} - s_0^{ij}) + E_1(d_{ij} - s_{ij}(t))\right] A_{ij} \, \hat{\mathbf{r}}_{ij}$$

**Verification**:
- At $t = 0$: $s_{ij} = s_0^{ij}$, giving $F = (E_0 + E_1)(d - s_0^{ij}) = E_u \cdot \varepsilon$ ✓
- As $t \to \infty$: $s_{ij} \to d_{ij}$, giving $F = E_0(d_{ij} - s_0^{ij}) \neq 0$ ✓ (solid-like)

The rest-length clamps were removed — the $E_0$ arm naturally prevents tissue collapse and provides the correct equilibrium resistance without any ad hoc numerical limits.

### Fix 2: Exact exponential integrator

The Forward Euler integrator was replaced with the **exact exponential integrator**. Writing $u = s_{ij} - d_{ij}$, the ODE $\dot{u} = -u/\tau$ has exact solution $u(t) = u_0 e^{-t/\tau}$, giving:

$$s_{ij}^{n+1} = d_{ij}^n + (s_{ij}^n - d_{ij}^n) \, e^{-\Delta t / \tau}$$

This is unconditionally stable for all $\Delta t > 0$ and $\tau > 0$. The computational cost is identical to Forward Euler (one `std::exp` call per pair per timestep).

### Code changes made

| Location | Change |
|----------|--------|
| `ViscoelasticGhostNode` struct | Added `initial_rest_lengths` field |
| `BuildNeighbourConnectivity` | Initialises `initial_rest_lengths` alongside `rest_lengths` |
| `UpdateGhostNodePositions` force | Replaced `k_u*(d-s)` with `E_0*(d-s0) + E_1*(d-s)` |
| `UpdateGhostNodePositions` integrator | Replaced Forward Euler with exact exponential integrator |
| `UpdateGhostNodePositions` | Removed rest-length clamps |
| `RemoveDepletedNodes` | Added `initial_rest_lengths` erase/clear alongside `rest_lengths` |
| `WriteOutput` | Updated `rest_length_strain` to use per-pair `initial_rest_lengths[nb]` instead of `mInitialSpacing` |
| File/class header comments | Updated to describe SLS correctly |
| `ViscoelasticGhostNodeEcmForce.hpp` header | Updated discrete analogue comment |

---

## Connection to Fertala et al. (2025)

The Fertala et al. ABAQUS model is:
- **Neo-Hookean backbone** → provides the permanent elastic stiffness $E_0$ (their "fully relaxed Young's modulus")
- **Single-term Prony series** → provides the decaying term $E_1 e^{-t/\tau}$

This is precisely an SLS. Our discrete split-force formulation maps this directly:
- $E_0$ arm with permanent rest length $s_0^{ij}$ → Neo-Hookean backbone
- $E_1$ arm with evolving rest length $s_{ij}(t)$ → Prony series dashpot

The Stage 2 single-evolving-rest-length formulation had no equivalent of the Neo-Hookean backbone: it was mathematically equivalent to a Maxwell element with stiffness $k_u = E_0 + E_1$, not an SLS.

---

## Summary of Development Stages

| Stage | Model | Long-time force | Integrator | Status |
|-------|-------|-----------------|------------|--------|
| 1: on-lattice | Fixed rest-length spring | Always non-zero (elastic, no relaxation) | Forward Euler (position) | Superseded |
| 2: single evolving rest length | Maxwell fluid | → 0 (full stress decay) | Forward Euler (rest length, conditionally stable) | Superseded |
| 3: split-force SLS | Standard Linear Solid | $E_0(d - s_0^{ij}) \neq 0$ | Exact exponential (unconditionally stable) | **Current** |

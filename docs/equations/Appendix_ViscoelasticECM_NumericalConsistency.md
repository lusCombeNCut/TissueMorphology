# Appendix A — Numerical Consistency of the Viscoelastic ECM Scheme

**Appendix to:** *Computational Modelling of Tissue Morphology using a Viscoelastic Extracellular Matrix*

---

## A.1 Overview

This appendix provides a rigorous mathematical demonstration that the discrete numerical
scheme implemented in `ViscoelasticGhostNodeEcmField` and `ViscoelasticGhostNodeEcmForce`
is consistent with the Standard Linear Solid (SLS) constitutive equations from which it is
derived. Specifically, it is shown that:

1. The discrete force law follows directly from the SLS mechanical model (Section A.3).
2. The internal state variable ODE governing rest length evolution is derived from the dashpot
   constitutive equation of the Maxwell arm (Section A.4).
3. The exact exponential integrator is the analytic solution to this ODE and not a numerical
   approximation (Section A.5).
4. A step-strain consistency test confirms that the scheme recovers the correct relaxation
   modulus $E(t) = E_0 + E_1 e^{-t/\tau}$ (Section A.6).
5. The scheme is unconditionally stable for all admissible timestep sizes (Section A.7).

---

## A.2 Rheological Model

The ghost–ghost ECM springs are modelled as a **Standard Linear Solid** (Zener model),
consisting of an equilibrium spring $E_0$ placed in parallel with a **Maxwell arm** — a
spring $E_1$ in series with a dashpot of viscosity $\eta$ (Figure A.1).

```
        ── E_0 ──────────────────────────
       │                                 │
 ──────┤                                 ├──────
       │                                 │
        ── E_1 ──── [η] ─────────────────
```
*Figure A.1. Standard Linear Solid (Zener model). An equilibrium spring $E_0$ acts in
parallel with a Maxwell arm comprising spring $E_1$ in series with dashpot $\eta$.*

For a connected ghost node pair $(i,j)$, define:

| Symbol | Meaning | Code variable |
|--------|---------|---------------|
| $d_{ij}(t)$ | Current node separation | `dist` |
| $s^0_{ij}$ | Permanent rest length (set at initialisation) | `initial_rest_lengths[nb]` |
| $s_{ij}(t)$ | Evolving dashpot rest length | `rest_lengths[nb]` |
| $E_0$ | Relaxed (equilibrium) stiffness | `mRelaxedStiffness` |
| $E_1$ | Relaxation modulus | `mRelaxationModulus` |
| $\tau$ | Relaxation time, $\tau = \eta / E_1$ | `mRelaxationTime` |

---

## A.3 Split-Force Constitutive Law

At any instant, the force carried by the equilibrium arm is:

$$F^{(0)}_{ij} = E_0 \left( d_{ij} - s^0_{ij} \right)$$

The force carried by the Maxwell arm spring (equal to the dashpot force, since they are in
series) is:

$$F^{(1)}_{ij} = E_1 \left( d_{ij} - s_{ij}(t) \right)$$

Summing both parallel arms gives the **total pair force**:

$$\boxed{F_{ij}(t) = E_0 \left( d_{ij} - s^0_{ij} \right) + E_1 \left( d_{ij} - s_{ij}(t) \right)}
\tag{A.1}$$

This is the split-force formulation implemented at
[`ViscoelasticGhostNodeEcmField.hpp:412–413`](../src/forces/ViscoelasticGhostNodeEcmField.hpp#L412-L413).

**Limiting behaviour** confirms physical correctness:

- **Instantaneous response** ($t = 0$): At initialisation $s_{ij}(0) = s^0_{ij}$, so
  $F = (E_0 + E_1)(d - s^0_{ij})$, giving unrelaxed stiffness $E_u = E_0 + E_1$.
- **Long-time response** ($t \to \infty$): As $s_{ij} \to d_{ij}$, the Maxwell arm
  contribution vanishes and $F \to E_0(d_{ij} - s^0_{ij})$, giving residual solid-like
  stiffness $E_0 \neq 0$.

---

## A.4 Derivation of the Internal State Variable ODE

The Maxwell arm places spring $E_1$ and dashpot $\eta$ in **series**; they therefore carry
the same force $F^{(1)}_{ij}$.

The dashpot constitutive law states that the dashpot force equals viscosity times the rate of
change of its internal length. Since the dashpot length is the component of the arm
deformation not stored in the spring — that is, $s_{ij}(t)$ — this gives:

$$F^{(1)}_{ij} = \eta\,\dot{s}_{ij}$$

Setting equal to the Maxwell arm spring force from (A.1):

$$\eta\,\dot{s}_{ij} = E_1 \left( d_{ij} - s_{ij} \right)$$

Dividing by $\eta$ and substituting $\tau = \eta / E_1$:

$$\boxed{\frac{ds_{ij}}{dt} = \frac{d_{ij} - s_{ij}}{\tau}}
\tag{A.2}$$

Equation (A.2) is a **first-order linear autonomous relaxation ODE**. It states that the
dashpot rest length $s_{ij}$ relaxes toward the instantaneous node separation $d_{ij}$ on the
characteristic timescale $\tau$. The internal variable $s_{ij}$ encodes the accumulated
inelastic strain of the dashpot; no explicit dashpot force is ever computed in the
implementation. This formulation follows the internal state variable approach of Taylor,
Pister & Goudreau (1970, Eq. 40) and Simo & Hughes (1998, §10.2).

---

## A.5 Exact Exponential Integrator

Over a single timestep $[t_n,\, t_n + \Delta t]$, the node separation $d_{ij}$ is treated as
**quasi-static** (frozen at $d^n_{ij}$). This approximation is valid when the timescale of
geometric change is large relative to $\tau$, which holds for the overdamped dynamics used
throughout the simulation.

Define the error variable $e = s_{ij} - d_{ij}$. Equation (A.2) becomes:

$$\dot{e} = -\frac{e}{\tau}$$

This is a pure exponential decay with closed-form solution:

$$e(t_n + \Delta t) = e(t_n)\,e^{-\Delta t / \tau}$$

Substituting back $e^n = s^n_{ij} - d^n_{ij}$:

$$\boxed{s^{n+1}_{ij} = d^n_{ij} + \left( s^n_{ij} - d^n_{ij} \right) e^{-\Delta t / \tau}}
\tag{A.3}$$

This is the **exact exponential integrator** — not a finite-difference approximation of (A.2)
but its analytic solution over the frozen-geometry interval. It is implemented at
[`ViscoelasticGhostNodeEcmField.hpp:443`](../src/forces/ViscoelasticGhostNodeEcmField.hpp#L443):

```cpp
double new_s = dist + (s_ij - dist) * std::exp(-dt / mRelaxationTime);
```

---

## A.6 Step-Strain Consistency Test

To verify that (A.1) and (A.3) together reproduce the declared relaxation modulus
$E(t) = E_0 + E_1 e^{-t/\tau}$, consider an idealised **step-strain** experiment.

**Setup.** Apply a constant strain $\varepsilon$ at $t = 0$ to a pair with unit natural length
($s^0_{ij} = 0$, $d_{ij} = \varepsilon = \text{const}$). The initial condition is
$s_{ij}(0) = 0$.

**Solution.** With $d_{ij}$ constant, (A.3) applied repeatedly over time gives the continuous
limit:

$$s_{ij}(t) = \varepsilon \left( 1 - e^{-t/\tau} \right)$$

**Force.** Substituting into (A.1):

$$F(t) = E_0\,\varepsilon + E_1 \left[ \varepsilon - \varepsilon\left(1 - e^{-t/\tau}\right) \right]
       = E_0\,\varepsilon + E_1\,\varepsilon\,e^{-t/\tau}$$

**Effective modulus:**

$$E(t) \equiv \frac{F(t)}{\varepsilon} = E_0 + E_1\,e^{-t/\tau} \qquad \checkmark
\tag{A.4}$$

Equation (A.4) is identical to the constitutive model stated in the file header of
`ViscoelasticGhostNodeEcmField.hpp` (line 10), confirming full consistency between the
continuous constitutive law and the discrete implementation.

---

## A.7 Numerical Stability

The stability of the exact integrator follows immediately from (A.3). The **amplification
factor** — the ratio $|e^{n+1}|/|e^n|$ — is:

$$G = e^{-\Delta t/\tau}$$

Since $\Delta t > 0$ and $\tau > 0$, it follows that $G \in (0, 1)$ for all admissible
parameter values. The scheme is therefore **unconditionally stable**: errors decay
monotonically regardless of timestep size.

For comparison, the Forward Euler discretisation of (A.2) gives:

$$s^{n+1}_{ij} = s^n_{ij} + \frac{\Delta t}{\tau}\left(d^n_{ij} - s^n_{ij}\right)$$

with amplification factor $G_{\text{FE}} = 1 - \Delta t/\tau$. Stability requires
$|G_{\text{FE}}| \leq 1$, i.e. $\Delta t \leq 2\tau$ — a **conditional** stability
constraint that would impose an additional restriction on the simulation timestep. The exact
integrator avoids this entirely.

---

## A.8 Summary

Table A.1 summarises the correspondence between each property of the continuous SLS model and
its discrete implementation.

*Table A.1. Correspondence between the continuous SLS constitutive model and the discrete
numerical scheme.*

| Property | Continuous model | Discrete scheme | Consistent? |
|----------|-----------------|-----------------|-------------|
| Force law | $F = E_0(d-s^0) + E_1(d-s)$ | Eq. (A.1); code lines 412–413 | ✓ |
| Internal variable ODE | $\dot{s} = (d-s)/\tau$ | Derived from dashpot law; Eq. (A.2) | ✓ |
| Time integrator | Analytic solution of ODE | Exact exponential integrator; Eq. (A.3) | ✓ |
| Instantaneous stiffness | $E_u = E_0 + E_1$ | $s(0) = s^0 \Rightarrow F = E_u(d - s^0)$ | ✓ |
| Long-time stiffness | $E_\infty = E_0$ | $s \to d \Rightarrow F = E_0(d - s^0)$ | ✓ |
| Relaxation modulus | $E(t) = E_0 + E_1 e^{-t/\tau}$ | Step-strain test; Eq. (A.4) | ✓ |
| Numerical stability | Unconditional | $G = e^{-\Delta t/\tau} \in (0,1)$ | ✓ |

---

## References

- Taylor, R. L., Pister, K. S., & Goudreau, G. L. (1970). Thermomechanical analysis of
  viscoelastic solids. *International Journal for Numerical Methods in Engineering*, 2(1),
  45–59. [Eq. 40 — original recursive exponential update for viscoelastic internal variables]

- Simo, J. C., & Hughes, T. J. R. (1998). *Computational Inelasticity*. Springer, New York.
  [§10.2 — internal state variable formulation for finite viscoelasticity]

- Liu, Y., et al. (2024). A computational framework for finite strain viscoelasticity.
  *Computer Methods in Applied Mechanics and Engineering*, 429, 117157.
  [§2.5, Remark 9 — internal state variable $\Gamma$ as tensor generalisation of $s_{ij}$]

- Fertala, A., et al. (2025). Mechanical properties of collagen type I fibrils.
  *bioRxiv* 2025.07.02.662292. [Constitutive model $E(t) = E_0 + E_1 e^{-t/\tau}$]

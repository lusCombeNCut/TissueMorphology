# Unit Conversion Guide: Translating CHASTE Simulations to Physical Units

This guide explains how to convert simulation parameters from **natural units** (cell diameters and hours) to **SI units** (pascals, newtons, etc.) for comparison with experimental data.

---

## Fundamental Unit System

### Simulation Units
- **Length:** Cell diameter (CD), where 1 CD ≈ **10 µm** (typical epithelial cell)
- **Time:** Hours (h), where 1 h = 3600 s
- **Force:** Arbitrary — absorbed into the damping coefficient η = 1

### CHASTE Equation of Motion

CHASTE solves the overdamped equation of motion with $\eta_{\text{sim}} = 1$:

$$\mathbf{v}_{\text{sim}} = \mathbf{F}_{\text{sim}}$$

The physical mapping is:

$$x_{\text{phys}} = x_{\text{sim}} \cdot L_0 \qquad (L_0 = \text{cell diameter in metres})$$

$$t_{\text{phys}} = t_{\text{sim}} \cdot T_0 \qquad (T_0 = 3600\ \text{s})$$

$$v_{\text{phys}} = v_{\text{sim}} \cdot \frac{L_0}{T_0}$$

Substituting into the physical drag equation $F_{\text{phys}} = \eta_{\text{phys}} \cdot v_{\text{phys}}$:

$$\boxed{F_{\text{phys}} = F_{\text{sim}} \cdot \eta_{\text{phys}} \cdot \frac{L_0}{T_0}}$$

This is the fundamental conversion relation.

---

## Conversion Factors

### Stiffness (Hooke's Law: $F = k \cdot \Delta x$)

$$k_{\text{phys}}\ [\text{N/m}] = \frac{F_{\text{phys}}}{\Delta x_{\text{phys}}}
= \frac{F_{\text{sim}} \cdot \eta_{\text{phys}} \cdot L_0/T_0}{\Delta x_{\text{sim}} \cdot L_0}
= k_{\text{sim}} \cdot \frac{\eta_{\text{phys}}}{T_0}$$

**Note:** $L_0$ cancels out — stiffness conversion does NOT depend on cell size directly, only on $\eta_{\text{phys}}$ and $T_0$.

### Pressure (Force per area: $P = F / A$)

For parameters with units [force/CD²]:

$$P_{\text{phys}}\ [\text{Pa}] = \frac{F_{\text{phys}}}{A_{\text{phys}}}
= \frac{F_{\text{sim}} \cdot \eta_{\text{phys}} \cdot L_0/T_0}{A_{\text{sim}} \cdot L_0^2}
= P_{\text{sim}} \cdot \frac{\eta_{\text{phys}}}{T_0 \cdot L_0}$$

### Bending Stiffness (Moment: $\kappa = \text{force} \times \text{length}$)

For parameters with units [force $\times$ CD]:

$$\kappa_{\text{phys}}\ [\text{N} \cdot \text{m}] = F_{\text{phys}} \cdot L_{\text{phys}}
= \left(F_{\text{sim}} \cdot \eta_{\text{phys}} \cdot \frac{L_0}{T_0}\right) \cdot (L_{\text{sim}} \cdot L_0)
= \kappa_{\text{sim}} \cdot \frac{\eta_{\text{phys}} \cdot L_0^2}{T_0}$$

---

## Estimating Effective Tissue Damping ($\eta_{\text{phys}}$)

**Do not use Stokes drag** ($\eta_{\text{Stokes}} = 6\pi\mu r \approx 9.4 \times 10^{-5}\ \text{N}\cdot\text{s}/\text{m}$ for a 5 µm cell in water). This assumes a free sphere in fluid, which is orders of magnitude wrong for epithelial cells adhered to ECM via focal adhesions.

**Use traction force data instead:**

Epithelial cells move at ~10 µm/h by exerting ~10 nN traction forces. Therefore:

$$\eta_{\text{phys}} = \frac{F_{\text{ref}}}{v_{\text{ref}}}
= \frac{10 \times 10^{-9}\ \text{N}}{10 \times 10^{-6}\ \text{m} / 3600\ \text{s}}
= \frac{10^{-8}\ \text{N}}{2.78 \times 10^{-9}\ \text{m}\cdot\text{s}^{-1}}
\approx 3.6\ \text{N}\cdot\text{s}/\text{m}$$

This is ~40,000× higher than Stokes drag, reflecting the dominant role of focal adhesion coupling rather than fluid viscosity.

**Default JSON values** (`referenceForcenN: 10`, `referenceVelocityUmPerH: 10`) give $\eta_{\text{phys}} \approx 3.6\ \text{N}\cdot\text{s}/\text{m}$.

---

## Worked Examples

Using defaults: $\eta_{\text{phys}} = 3.6\ \text{N}\cdot\text{s}/\text{m}$, $T_0 = 3600\ \text{s}$, $L_0 = 10^{-5}\ \text{m}$.

### ECM Stiffness (`ecmStiffness = 50.0`)

$$k_{\text{phys}} = 50.0 \times \frac{3.6}{3600} = 0.05\ \text{N/m}$$

Converting to elastic modulus (crypt geometry: $L = 50\ \mu\text{m}$, $A = 100\ \mu\text{m}^2$):

$$E = \frac{k_{\text{phys}} \cdot L}{A}
= \frac{0.05 \times 50 \times 10^{-6}}{100 \times 10^{-12}}
= 25{,}000\ \text{Pa} = 25\ \text{kPa}$$

**Literature range for intestinal mucosa: 5–15 kPa.** The default of 50.0 is slightly stiff (25 kPa ≈ fibrotic range). Reducing `ecmStiffness` to ~20–30 would target healthy intestine.

### Lumen Pressure (`lumenPressure = 5.0`)

$$P_{\text{phys}} = 5.0 \times \frac{3.6}{3600 \times 10^{-5}}
= 5.0 \times 100
= 500\ \text{Pa}$$

**Biological cavity pressures: 100–1000 Pa.** 500 Pa is right in range.

### Spring Stiffness (`springStiffness = 15.0`)

$$k_{\text{phys}} = 15.0 \times \frac{3.6}{3600} = 0.015\ \text{N/m}$$

### Bending Stiffness (`bendingStiffness = 30.0`)

$$\kappa_{\text{phys}} = 30.0 \times \frac{3.6 \times (10^{-5})^2}{3600}
= 30.0 \times \frac{3.6 \times 10^{-10}}{3600}
= 30.0 \times 10^{-13}
= 3 \times 10^{-12}\ \text{N\cdot m}$$

---

## Python Helper Class

```python
import numpy as np

class ChasteTissueUnitConverter:
    """
    Convert CHASTE simulation parameters to SI units.

    Derivation:
        F_phys = F_sim * eta_phys * (L0/T0)
        k_phys [N/m]  = k_sim  * eta_phys / T0
        P_phys [Pa]   = P_sim  * eta_phys / (T0 * L0)
        kappa_phys [N·m] = kappa_sim * eta_phys * (L0**2) / T0
    """

    def __init__(self, cell_diameter_um=10, reference_force_nN=10, reference_velocity_um_per_h=10):
        """
        Parameters
        ----------
        cell_diameter_um : float
            Cell diameter in µm (default: 10 µm)
        reference_force_nN : float
            Typical epithelial traction force in nN (default: 10 nN)
        reference_velocity_um_per_h : float
            Typical migration speed in µm/h (default: 10 µm/h = 1 CD/h)
        """
        self.L0 = cell_diameter_um * 1e-6   # metres
        self.T0 = 3600.0                     # seconds per hour

        # Effective tissue damping from traction data
        f_N  = reference_force_nN * 1e-9
        v_ms = (reference_velocity_um_per_h * 1e-6) / self.T0
        self.eta_phys = f_N / v_ms          # N·s/m

        print(f"Effective tissue damping η = {self.eta_phys:.3f} N·s/m")

    def stiffness_to_N_per_m(self, k_sim):
        """k_sim [force/CD] → N/m"""
        return k_sim * (self.eta_phys / self.T0)

    def stiffness_to_elastic_modulus_kPa(self, k_sim, domain_length_um=50):
        """k_sim [force/CD] → elastic modulus in kPa.
        Uses E = k * L / A where A = cell cross-section area."""
        k_SI = self.stiffness_to_N_per_m(k_sim)
        L = domain_length_um * 1e-6
        A = self.L0 ** 2
        return (k_SI * L / A) / 1e3

    def pressure_to_Pa(self, p_sim):
        """p_sim [force/CD²] → Pa"""
        return p_sim * (self.eta_phys / (self.T0 * self.L0))

    def bending_stiffness_to_N_m(self, kappa_sim):
        """kappa_sim [force × CD] → N·m (flexural rigidity)"""
        return kappa_sim * (self.eta_phys * self.L0 ** 2 / self.T0)

    def force_to_N(self, f_sim):
        """f_sim → N (for a body moving at 1 CD/h)"""
        return f_sim * self.eta_phys * (self.L0 / self.T0)


# Example: default parameters
c = ChasteTissueUnitConverter()

print(f"\necmStiffness=50.0:    {c.stiffness_to_N_per_m(50.0):.4f} N/m  →  {c.stiffness_to_elastic_modulus_kPa(50.0):.1f} kPa")
print(f"springStiffness=15.0: {c.stiffness_to_N_per_m(15.0):.4f} N/m")
print(f"lumenPressure=5.0:    {c.pressure_to_Pa(5.0):.1f} Pa")
print(f"bendingStiffness=30.0:{c.bending_stiffness_to_N_m(30.0):.2e} N·m")
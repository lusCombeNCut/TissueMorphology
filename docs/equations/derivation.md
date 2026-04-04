This is a standard first-order linear ODE. Here's the full derivation:

The ODE:

$$\frac{ds}{dt} = \frac{d - s}{\tau}$$

where $d$ (current distance) is treated as a constant over one timestep.

Step 1 — Rearrange into standard form

$$\frac{ds}{dt} + \frac{s}{\tau} = \frac{d}{\tau}$$

This is a first-order linear ODE of the form $s' + P(t)s = Q(t)$ with $P = 1/\tau$, $Q = d/\tau$.

Step 2 — Integrating factor

Multiply both sides by $\mu(t) = e^{t/\tau}$:

$$e^{t/\tau}\frac{ds}{dt} + \frac{e^{t/\tau}}{\tau}s = \frac{d}{\tau}e^{t/\tau}$$

The left side is the exact derivative of a product:

$$\frac{d}{dt}\left(s \cdot e^{t/\tau}\right) = \frac{d}{\tau}e^{t/\tau}$$

Step 3 — Integrate both sides

$$\int \frac{d}{dt}\left(s \cdot e^{t/\tau}\right),dt = \int \frac{d}{\tau}e^{t/\tau},dt$$

$$s \cdot e^{t/\tau} = d \cdot e^{t/\tau} + C$$

Step 4 — Solve for $s(t)$

Divide by $e^{t/\tau}$:

$$s(t) = d + C e^{-t/\tau}$$

Step 5 — Apply initial condition $s(0) = s^n$

At $t = 0$:

$$s^n = d + C e^{0} = d + C \implies C = s^n - d$$

Substituting back:

$$\boxed{s(t) = d + (s^n - d),e^{-t/\tau}}$$

Step 6 — Evaluate at $t = \Delta t$ to get the timestep update

$$s^{n+1} = d + (s^n - d),e^{-\Delta t/\tau}$$

This is exactly line 443:
# Guiding Center

By solving the trajectories of particles, we can calculate the actual guiding center orbits by following the definition. This is supported directly via `get_gc`. `TestParticle.jl` also provides several methods to trace the guiding center, serving different purposes.

## Diagnostic Methods (Reference-based)

These methods require a full particle trajectory solution to calculate drifts. They are useful for analyzing which drift components dominate or checking the validity of drift approximations.

* [`trace_gc_drifts!`](@ref): Calculates the guiding center trajectory using standard analytical drift formulas (ExB, Gradient-B, Curvature drifts). It uses the parallel and perpendicular velocities from the full particle simulation.
* [`trace_gc_exb!`](@ref): Calculates the trajectory considering only the $\mathbf{E} \times \mathbf{B}$ drift and parallel motion. Useful for identifying ExB dominance.
* [`trace_gc_flr!`](@ref): Includes Finite Larmor Radius (FLR) corrections to the ExB drift, suitable for nonuniform electric fields.

## Self-consistent Solvers (GCA)

[`trace_gc!`](@ref) solves the 1st order Guiding Center Approximation (GCA) equations. This approximation holds when the characteristic scale of field variations ``L`` (often represented by the curvature radius) is much larger than the Larmor radius ``\rho`` (``L \gg \rho``).

**State Variables**: The solver evolves the 4D state ``[\mathbf{R}, v_\parallel]``:

- ``\mathbf{R}``: The spatial position of the guiding center.
- ``v_{\parallel}``: The velocity parallel to the magnetic field.

**Parameters**:

- ``\mu = p_{\perp}^2 / (2mB)``: The magnetic moment (adiabatic invariant). It is calculated from initial conditions and treated as a constant parameter.

**Effective Fields**:

To incorporate geometric drifts (curvature) and mirror forces naturally, we define the effective magnetic and electric fields:

```math
\begin{aligned}
\mathbf{B}^* &= \mathbf{B} + \frac{m v_{\parallel}}{q} \nabla \times \mathbf{b} \\
\mathbf{E}^* &= \mathbf{E} - \frac{1}{q} \nabla ( \mu B )
\end{aligned}
```
where ``\mathbf{b} = \mathbf{B}/B`` is the unit vector along the magnetic field.

**Equations of Motion**:

The time evolution of the guiding center is governed by:

```math
\begin{aligned}
\dot{\mathbf{R}} &= \frac{1}{B^*_{\parallel}} \left[ v_{\parallel} \mathbf{B}^* + \mathbf{E}^* \times \mathbf{b} \right] \\
\dot{v}_{\parallel} &= \frac{q}{m B^*_{\parallel}} \mathbf{B}^* \cdot \mathbf{E}^*
\end{aligned}
```
where ``B^*_{\parallel} = \mathbf{b} \cdot \mathbf{B}^*``.

* Term 1 of ``\dot{\mathbf{R}}`` (``\propto \mathbf{B}^\ast``) represents motion along the field line plus curvature drift.
* Term 2 of ``\dot{\mathbf{R}}`` (``\propto \mathbf{E}^\ast \times \mathbf{b}``) represents the ``\mathbf{E} \times \mathbf{B}`` drift and the ``\nabla B`` drift.

In theoretical treatments, individual drift terms are often derived separately (see below). By solving the GCA equations, all these drifts are captured self-consistently to first order.

## Summary of guiding center drifts

General force:

```math
\mathbf{v}_f = \frac{1}{q}\frac{\mathbf{F}\times\mathbf{B}}{B^2}
```

Electric field:

```math
\mathbf{v}_E = \frac{\mathbf{E}\times\mathbf{B}}{B^2}
```

Gravitational field:

```math
\mathbf{v}_g = \frac{m}{q}\frac{\mathbf{g}\times\mathbf{B}}{B^2}
```

Nonuniform electric field:

```math
\mathbf{v}_E = \Big( 1+\frac{1}{4}r_L^2 \nabla^2 \Big)\frac{\mathbf{E}\times\mathbf{B}}{B^2}
```

Nonuniform magnetic field:

Grad-B:

```math
\begin{aligned}
\mathbf{v}_{\nabla B} &= \pm \frac{1}{2}v_\perp r_L\frac{\mathbf{B}\times\nabla B}{B^2} \\
&=\frac{mv_\perp^2}{2qB^3}\mathbf{B}\times\nabla B
\end{aligned}
```

Curvature drift:

```math
\begin{aligned}
\mathbf{v}_c &= \frac{mv_\parallel^2}{q}\frac{\mathbf{R}_c \times\mathbf{B}}{R_c^2 B^2} \\
&= \frac{mv_\parallel^2}{qB^4}\mathbf{B}\times\left[ \mathbf{B}\cdot\nabla\mathbf{B} \right] \\
&= \frac{m v_\parallel^2}{qB^3}\mathbf{B}\times\nabla B\,\text{(current-free)}
\end{aligned}
```

Curved vacuum field:

```math
\begin{aligned}
\mathbf{v}_c + \mathbf{v}_{\nabla B} &= \frac{m}{q}\Big( v_\parallel^2 + \frac{1}{2}v_\perp^2 \Big) \frac{\mathbf{R}_c \times\mathbf{B}}{R_c^2 B^2} \\
&= \frac{m}{q}\frac{\mathbf{B}\times\nabla B}{B^3}\Big( v_\parallel^2 + \frac{1}{2}v_\perp^2 \Big)
\end{aligned}
```

Polarization drift:[^2]

```math
\mathbf{v}_p = \pm \frac{1}{\omega_c B}\frac{d\mathbf{E}}{dt}
```

[^2]: In common test particle models, we assume static EM fields, so the polarization drift as well as adiabatic heating is not present. However, it is easily achieveable here in this package.

## Adiabatic Invariants

See more thorough [notes](https://henry2004y.github.io/KeyNotes/contents/single.html#sec-adiabatic-invariant) on the adiabatic invariants.

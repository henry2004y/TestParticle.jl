# Tutorial

What makes plasmas particularly difficult to analyze is the fact that the densities fall in an intermediate range. Fluids like water are so dense that the motions of individual molecules do not have to be considered. Collisions dominate, and the simple equations of ordinary fluid dynamics suffice. At the other extreme in very low-density devices, only single-particle trajectories need to be considered; collective effects are often unimportant. Plasma behaves sometimes like fluids, and sometimes like a collection of individual particles. The first step in learning how to deal with this schizophrenic personality is to understand how single particles behave in electric and magnetic fields.

Here we assume that the EM fields are prescribed and not affected by the charged particles. References can be found in classic textbooks like [Introduction to Plasma Physics and Controlled Fusion](https://link.springer.com/book/10.1007/978-3-319-22309-4) by F.F.Chen, and [Fundamentals of Plasma Physics](https://doi.org/10.1017/CBO9780511807183) by Paul Bellan. For more complete notes corresponding to the derivation online, please check out [Single-Particle Motions](https://henry2004y.github.io/KeyNotes/contents/single.html).

## Choice of numerical algorithms

By default DifferentialEquations.jl applies `Tsit5` to an ODE problem.
However, it is not always guaranteed to work. For example, the demo case of electron tracing in the magnetic bottle with strong magnetic field is tested to work only with fixed timestep algorithms like `Euler` and the Adams-Bashforth family.
Take you some time to figure out which algorithm works for your problem!

## Multiple particles tracing

There are two ways to trace multiple particles simultaneously:

1. Extracting the solution in a loop with varying initial conditions. See the example [demo_ensemble](@ref).
2. Constructing the [Ensemble Simulations](https://diffeq.sciml.ai/stable/features/ensemble/). One example can be found [here](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_ensemble.jl). However, note that currently the ensemble type replicates the parameters for each solution, which is very memory inefficient for tracing in a numeric field.

## Summary of Guiding Center Drifts

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
\mathbf{v}_{\nabla B} = \pm \frac{1}{2}v_\perp r_L\frac{\mathbf{B}\times\nabla B}{B^2}
```

Curvature drift:

```math
\mathbf{v}_R = \frac{mv_\parallel^2}{q}\frac{\mathbf{R}_c \times\mathbf{B}}{R_c^2 B^2}
```

Curved vacuum field:

```math
\mathbf{v}_R + \mathbf{v}_{\nabla B} = \frac{m}{q}\Big( v_\parallel^2 + \frac{1}{2}v_\perp^2 \Big) \frac{\mathbf{R}_c \times\mathbf{B}}{R_c^2 B^2}
```

Polarization drift:[^2]

```math
\mathbf{v}_p = \pm \frac{1}{\omega_c B}\frac{d\mathbf{E}}{dt}
```

[^2]: In common test particle models, we assume static EM fields, so the polarization drift as well as adiabatic heating is not present. However, it is easily achieveable here in this package.

## Adiabatic Invariants

See more thorough [notes](https://henry2004y.github.io/KeyNotes/contents/single.html#sec-adiabatic-invariant) on the adiabatic invariants.

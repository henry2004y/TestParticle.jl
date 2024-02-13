# Tutorial

What makes plasmas particularly difficult to analyze is the fact that the densities fall in an intermediate range. Fluids like water are so dense that the motions of individual molecules do not have to be considered. Collisions dominate, and the simple equations of ordinary fluid dynamics suffice. At the other extreme in very low-density devices, only single-particle trajectories need to be considered; collective effects are often unimportant. Plasma behaves sometimes like fluids, and sometimes like a collection of individual particles. The first step in learning how to deal with this schizophrenic personality is to understand how single particles behave in electric and magnetic fields.

Here we assume that the EM fields are prescribed and not affected by the charged particles. References can be found in classic textbooks like [Introduction to Plasma Physics and Controlled Fusion](https://link.springer.com/book/10.1007/978-3-319-22309-4) by F.F.Chen, and [Fundamentals of Plasma Physics](https://doi.org/10.1017/CBO9780511807183) by Paul Bellan. For more complete notes corresponding to the derivation online, please check out [Single-Particle Motions](https://henry2004y.github.io/KeyNotes/contents/single.html).

## Choice of numerical algorithms

By default DifferentialEquations.jl applies `Tsit5` to an ODE problem. If only OrdinaryDiffEq.jl is imported, then we need to explicitly select a solver.
However, not all solvers are guaranteed to work. For example, the demo case of electron tracing in the magnetic bottle with strong magnetic field is tested to work only with fixed timestep algorithms like `Euler` and the Adams-Bashforth family.

Currently we recommend `Vern9` as a starting point for adaptive timestepping, with additional fine tuning by changing `reltol` if needed. You can also try out the native implementation of the Boris method in TestParticle.jl, with a constraint of using a fixed time step.

Take you some time to figure out which algorithm works for your problem!

## Unit conversions

By default SI units are applied within the package. However, users can also define their own units by setting the particle mass and charge (2 constants) and providing the basic scales of magnetic field, length, time, or velocity (3 reference scales required; length, time and velocity are interchangable).

In the dimensionless natural unit system, we have the most number of ones which in turn gives the simplest form for computation. Check out the demos on unit conversions for more details.

## Tracing backwards in time

It is easy to trace a charged particle backwards in time. In a forward tracing problem, we set `tspan = (t1, t2)` where `t1 < t2`. In a backward tracing problem, we simply set `t1 > t2`, e.g. `tspan = (0.0, -1.0)`. Note that tracing backwards in time is different from inversing the velocity because of the cross product in the Lorentz force. More specifically, the drift velocities will not change sign if one inverses velocity.

## Multiple particles tracing

There are two ways to trace multiple particles simultaneously:

1. Extracting the solution in a loop with varying initial conditions. See the example [demo_ensemble](@ref).
2. Constructing the [Ensemble Simulations](https://diffeq.sciml.ai/stable/features/ensemble/). One example can be found [here](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_ensemble.jl). However, note that by default the ensemble type replicates the parameters for each solution, which is very memory inefficient for tracing in a numeric field. We need to set `safetycopy=false` to make the field as a reference in the parameter of each trajectory.

The Boris pusher follows a similar interface for multithreading.

## Guiding center drifts

By solving the trajectories of particles, we can calculate the actual guiding center orbits by following the definition. This is supported directly via `get_gc`. In theoretical treatments, all kinds of drifts have been derived with analytical formulas listed below. With additional ODEs for solving the drifts, we can separate the different effects or check the deviation of actual orbits from the low order analytical formulas.

### Summary of guiding center drifts

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

## Presentations

Please checkout:

* [TestParticle.jl: A New Tool for An Old Problem](https://henry2004y.github.io/pluto_playground/testparticle_202401.html)
* [基于开源工具链的测试粒子模型](https://henry2004y.github.io/pluto_playground/testparticle_202212.html)
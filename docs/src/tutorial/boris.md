# [Boris Pusher](@id Boris-Pusher)

The native Boris pusher is a widely used algorithm for tracing charged particles in magnetic fields. It is second order, volume preserving, and exactly conserves energy in a static, uniform magnetic field, while being highly efficient for long-time tracing. A basic usage example is shown in [Boris Method](@ref); this page covers the full family of Boris solvers available in `TestParticle.jl`.

## 1. Standard Boris Method

The Boris method advances the particle in a staggered (leapfrog) manner: the velocity is first rotated by the magnetic field, then the position is pushed using the updated velocity, and the process repeats. Because the magnetic rotation is applied as an *exact* rotation rather than an approximate update, energy is conserved to machine precision in a static, uniform field.

The method requires a fixed time step $\Delta t$ much smaller than the gyroperiod $T = 2\pi / \Omega$ to accurately capture the orbital phase. The phase error scales as $\mathcal{O}(\Delta t^2)$, so a finer step is needed when phase fidelity matters more than wall-clock time.

To reach extremely high fidelity over long time scales without shrinking $\Delta t$ drastically, `TestParticle.jl` implements several advanced variants: **Multistep Boris**, **Hyper Boris**, and **Adaptive Boris**.

## 2. Multistep Boris ($n$-step Boris)

The Multistep Boris integrator artificially subdivides the standard update cycle into $n$ smaller sub-steps. The electric and magnetic fields are evaluated once at the full time step location, but the velocity rotation is broken into $n$ smaller partial rotations.

Let the normalized rotation vector for the standard Boris step be:
```math
\mathbf{t} = \frac{q \Delta t}{2 m} \mathbf{B}
```

In the Multistep Boris scheme, the rotation vector is divided by $n$:
```math
\mathbf{t}_n = \frac{1}{n} \mathbf{t}
```
The rotation part of the Boris method is recursively applied $n$ times using $\mathbf{t}_n$. While reducing the numerical detuning and improving rotation accuracy, it does not intrinsically change the mathematical order of the phase error (it remains 2nd order) but substantially decreases the absolute error coefficient.

## 3. Hyper Boris (Gyrophase Correction)

The Hyper Boris integrator ([Zenitani & Kato, 2025](https://arxiv.org/abs/2505.02270)) achieves higher-order accuracy in tracking the precise gyro-orbital phase by introducing specialized correction factors to the electric and magnetic drift terms.

Given an order $N$ ($N \in \{4, 6\}$), the normalized vectors are modified prior to the rotation:
```math
\begin{aligned}
\mathbf{t}_n^\prime &= f_N \, \mathbf{t}_n \\
\mathbf{e}_n^\prime &= f_N \, \mathbf{e}_n + c_N (\mathbf{e}_n \cdot \mathbf{t}_n) \mathbf{t}_n
\end{aligned}
```
where $\mathbf{e}_n = \frac{q \Delta t}{2 m n} \mathbf{E}$.

The factors $f_N$ and $c_N$ are Taylor-expanded coefficients based on the magnitude $t_{mag}^2 = |\mathbf{t}_n|^2$:

**4th-Order Hyper Boris ($N=4$):**
```math
\begin{aligned}
f_4 &= 1 + \frac{1}{3} t_{mag}^2 \\
c_4 &= -\frac{1}{3}
\end{aligned}
```

**6th-Order Hyper Boris ($N=6$):**
```math
\begin{aligned}
f_6 &= 1 + \frac{1}{3} t_{mag}^2 + \frac{2}{15} t_{mag}^4 \\
c_6 &= -\frac{1}{3} - \frac{2}{15} t_{mag}^2
\end{aligned}
```

These analytically tuned correctors virtually eliminate phase disparity and energy fluctuation across standard $\Delta t$ bounds.

## 4. Adaptive Boris Method

The `AdaptiveBoris` solver adjusts the time step $\Delta t$ dynamically based on the local cyclotron frequency $\Omega_c = |q B / m|$. This is particularly useful in systems with strong magnetic field gradients, such as magnetic mirrors or planetary magnetospheres, where the required resolution varies significantly along the particle's trajectory.

The time step is determined by:
```math
\Delta t = \eta T_c = \eta \frac{2\pi}{\Omega_c}
```
where $\eta$ is a safety factor (typically between 0.01 and 0.1). This represents the ratio of the time step to the local gyroperiod.

### Maintaining Time Reversibility and Energy Conservation

The standard Boris method is a second-order, volume-preserving integrator that exactly conserves energy in a static, uniform magnetic field. These properties are closely linked to its **time-reversibility**. However, naive adaptive time-stepping usually breaks this reversibility because it disrupts the staggered "leapfrog" synchronization between position and velocity.

To preserve the leapfrog structure and maintain energy conservation, `TestParticle.jl` employs a **velocity resync** procedure whenever the time step is updated. This technique involves moving the velocity back to the node (position location) and then repositioning it based on the new time step. Given a step from $\mathbf{x}_n$ to $\mathbf{x}_{n+1}$ with $\Delta t_{old}$:

1.  **Advance position**: The position is updated to $\mathbf{x}_{n+1}$ using the velocity $\mathbf{v}_{n+1/2}$ and $\Delta t_{old}$.
2.  **Update time step**: A new $\Delta t_{new}$ is calculated based on the magnetic field at the new position $\mathbf{x}_{n+1}$.
3.  **Resynchronize velocity**: The velocity $\mathbf{v}_{n+1/2}^{old}$ (centered at $t_{n+1} - \frac{1}{2}\Delta t_{old}$) is moved to the node $t_{n+1}$ by a half-step Boris push, and then moved to a new half-step location $t_{n+1} - \frac{1}{2}\Delta t_{new}$ by a half-step backward push.

This ensures that the velocity is always correctly centered relative to the current $\Delta t$ before each step. This "re-centering" at the nodes allows the integrator to remain practically time-reversible and maintains excellent energy conservation even as the time step changes by orders of magnitude.

## 5. Using the Solvers

### Fixed-step Solvers

Fixed-step solvers are specified as the second argument to `solve` (after the problem) and require a `dt` keyword:

- `Boris()`: The standard Boris pusher.
- `MultistepBoris2(n)`: Sub-cycling division count `n`.
- `MultistepBoris4(n)`: 4th-order Hyper-Boris with sub-cycling `n`.
- `MultistepBoris6(n)`: 6th-order Hyper-Boris with sub-cycling `n`.

```julia
# Standard Boris
sol = TestParticle.solve(prob, Boris(); dt)

# Multistep Boris (n=2)
sol = TestParticle.solve(prob, MultistepBoris2(n=2); dt)

# Hyper Boris (N=4, n=2)
sol = TestParticle.solve(prob, MultistepBoris4(n=2); dt)
```

Combining both $n > 1$ and higher-order correction ($N > 2$) ensures ultra-high stability tracking over drastically varying gradient fields.

One trajectory returns an `ODESolution`, the same type any other SciML solver returns.

### Adaptive Boris

The adaptive solver adjusts the time step automatically based on the local gyroperiod.

```julia
# Adaptive Boris with safety factor 0.05 (20 steps per period)
alg = AdaptiveBoris(safety=0.05)
sol = TestParticle.solve(prob, alg)
```

### Ensembles

Multiple particles are traced through a SciML [`EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/), so the whole ensemble interface applies: `EnsembleThreads()`, `EnsembleDistributed()`, `EnsembleSplitThreads()`, `trajectories`, `seed`, `batch_size`, `pmap_batch_size`, and the `output_func` and `reduction` hooks.

```julia
eprob = EnsembleProblem(prob; prob_func, safetycopy = false)
sols = TestParticle.solve(eprob, Boris(), EnsembleThreads();
    dt, trajectories = 1000, seed = 1234)
```

`prob_func(prob, ctx)` prepares each trajectory. It receives an `EnsembleContext` whose `ctx.rng` is derived from `seed`, which keeps a run reproducible. Because the `prob_func` of a `TraceProblem` is used to build the same ensemble, the shorter form is equivalent:

```julia
prob = TraceProblem(stateinit, tspan, param; prob_func)
sols = TestParticle.solve(prob, Boris(), EnsembleThreads();
    dt, trajectories = 1000)
```

`sols.u` is then a vector of `ODESolution`s, one per trajectory.

The ensemble algorithm is a required third argument in the shorter form, because
`solve(prob, alg)` without one means a single trajectory. It cannot have a
default either: SciML calls exactly that two-argument form once per trajectory it
builds, so a default would ask every trajectory to build an ensemble of its own.

### Saving the Output

Every accepted step is saved by default. To choose the output times explicitly, pass `saveat`, either as a collection of times or as an interval:

```julia
# Save at the given times
sol = TestParticle.solve(prob, Boris(); dt, saveat = 0.0:5.0e-10:3.0e-8)

# Save every 5.0e-10 across the time span
sol = TestParticle.solve(prob, Boris(); dt, saveat = 5.0e-10)
```

`saveat` does not shorten any step, so the trajectory is identical to a run without it and only the reporting changes. Inside a step the state is interpolated **linearly**, which is the dense output these methods admit: a Boris step advances the position linearly with the half-step velocity $\mathbf{v}_{n+1/2}$, so the interpolant reproduces the stored positions exactly, whereas the velocity is only as accurate as the method itself. Requesting a value inside a step is therefore consistent with the solver's own order, not better or worse than the step values around it.

`save_start` and `save_end` (both `true` by default) add the ends of the time span to the requested times. `save_fields = true` and `save_work = true` keep appending their columns to every saved state.

!!! warning "Removed keyword: `savestepinterval`"
    `savestepinterval = k`, which saved every $k$-th step regardless of how the
    times were named, has been removed. Use `saveat = k * dt` to report the same
    times in a fixed-step run. Passing `savestepinterval` is now an ordinary
    unsupported-keyword error.

## 6. Where the solvers live

The Boris family is implemented in its own package,
[OrdinaryDiffEqBoris](https://github.com/henry2004y/TestParticle.jl/tree/master/lib/OrdinaryDiffEqBoris),
which follows the SciML convention for an algorithm package: the methods are
ordinary SciML algorithms, driven by the SciML loop, and they can be used on
their own. TestParticle.jl depends on it and adds the parts that are specific to
particle tracing: the `TraceProblem` container, the boundary callback, and the
field and work columns.

Used directly, the solvers only need a parameter container that answers three
questions, `get_q2m`, `get_EField` and `get_BField`. The default methods assume
the layout `(q2m, m, E, B, ...)` that `prepare` produces, and any other type can
support the solvers by adding methods to those three functions:

```julia
using OrdinaryDiffEqBoris, StaticArrays

param = (q2m, m, Efunc, Bfunc)   # or any type with the three accessors
prob = ODEProblem((u, p, t) -> nothing, SA[0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan, param)
sol = solve(prob, AdaptiveBoris(safety = 0.1); dt)
```

An `ODEProblem` needs a right-hand side `f(u, p, t)`, but a Boris method is a
map from one state to the next, not a differential equation, so `f` is never
called: at every step the method takes the charge-to-mass ratio and the fields
from `p`, and `(u, p, t) -> nothing` is enough. Adaptivity is selected by the
`adaptive` solve keyword, as for any other SciML solver, and the adaptive methods
follow the local gyroperiod rather than an error estimate, since a Boris step
forms none.

What TestParticle.jl still owns is the layer around the loop, in
`src/boris/boris_solve.jl`: the conversion of a `TraceProblem` into an
`ODEProblem`, the gyroperiod-based initial step of the adaptive methods, the
`isoutside` callback, and the appended field and work columns.

## 7. What a step costs, and running on a device

A Boris method carries the velocity at the half step, so a step advances the
staggered pair `(r, v_{n+1/2})` with the fields taken at the node it starts from.
To report `(r_{n+1}, v_{n+1})`, the node velocity is reconstructed with a half
Boris rotation.

In `OrdinaryDiffEqBoris`, carrying evaluated fields across step boundaries ensures
each step spends only one field evaluation. Furthermore, when `save_everystep=false`
(or when saving selectively without callbacks), node velocity is computed lazily
only for saved states, skipping redundant Boris rotations on intermediate steps.

In `OrdinaryDiffEqBoris`, fixed-step Boris solves on standard domains automatically take
a specialized fast path inside `solve!`. This eliminates generic SciML integrator loop
bookkeeping and preallocates trajectory arrays, achieving raw native performance (~5 ns/step
on analytic fields). Whenever adaptive stepping (`AdaptiveBoris`), boundary checking
(`isoutside`), intermediate save targets (`saveat`), or custom SciML callbacks are present,
`solve!` routes transparently through the full SciML integrator loop.

For pushing large ensembles on accelerators or multiple CPU threads, pass a
`KernelAbstractions` backend to `solve(prob, alg, backend, ...)`: see
[GPU Ensemble Tracing](@ref GPU-Ensemble-Tracing). The device runs fixed-step
solvers (`Boris()` and `MultistepBoris{N}`), while adaptive methods decide their step
sizes on the host through the SciML loop.

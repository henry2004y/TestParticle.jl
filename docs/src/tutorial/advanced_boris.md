# Advanced Boris Methods

The standard Boris method is a widely used algorithm for tracing charged particles because it exactly conserves energy in a static, uniform magnetic field. However, it requires a time step size $\Delta t$ much smaller than the gyroperiod $T = 2\pi / \Omega$ to accurately capture the orbital phase (the phase error scales as $\mathcal{O}(\Delta t^2)$).

To address scenarios requiring extremely high fidelity over long time scales without decreasing the simulation $\Delta t$ too drastically, `TestParticle.jl` implements two advanced variants: **Multistep Boris** and **Hyper Boris**.

## 1. Multistep Boris ($n$-step Boris)

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

## 2. Hyper Boris (Gyrophase Correction)

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

### Using the Advanced Solvers

You can access these features in `TestParticle.jl` by explicitly defining the parameters `n` and `N` in the `solve` parameters.

- `n`: Sub-cycling division count (default `1`).
- `N`: Evaluator correction order ($2, 4,$ or $6$). `2` denotes standard uncorrected Boris.

```julia
# Standard Boris (n=1, N=2)
sol = TestParticle.solve(prob; dt)

# Multistep Boris (n=2, N=2)
sol = TestParticle.solve(prob; dt, n=2)

# Hyper Boris (n=2, N=4)
sol = TestParticle.solve(prob; dt, n=2, N=4)
```

Combining both $n > 1$ and $N > 2$ ensures ultra-high stability tracking over drastically varying gradient fields, such as inside magnetic traps or collisionless shocks!

## 3. Adaptive Boris Method

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

## 4. Using the Advanced Solvers

### Fixed-step Advanced Boris

You can access multistep and hyper features in the default solver by explicitly defining the parameters `n` and `N` in the `solve` parameters.

- `n`: Sub-cycling division count (default `1`).
- `N`: Evaluator correction order ($2, 4,$ or $6$). `2` denotes standard uncorrected Boris.

```julia
# Multistep Boris (n=2, N=2)
sol = TestParticle.solve(prob; dt, n=2)

# Hyper Boris (n=2, N=4)
sol = TestParticle.solve(prob; dt, n=2, N=4)
```

### Adaptive Boris

You can use the adaptive solver by passing an `AdaptiveBoris` object as the second argument to `solve`.

```julia
# Adaptive Boris with safety factor 0.05 (20 steps per period)
alg = AdaptiveBoris(safety=0.05)
sol = TestParticle.solve(prob, alg)[1]
```

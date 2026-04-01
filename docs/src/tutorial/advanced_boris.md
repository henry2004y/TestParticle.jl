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

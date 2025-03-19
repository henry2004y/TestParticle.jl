# ---
# title: Cosmic Ray Tracing
# id: demo_cosmic_ray
# date: 2025-03-19
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.11.4
# description: Tracing cosmic ray charged particles
# ---

# This example shows how to trace cosmic rays in a background magnetic field. In practice, at least in the interstellar medium, the effect of the electric field over high-energy (multi TeV) cosmic rays can be neglected. Therefore, energy is conserved and we are interested in looking at the scattering and diffusion processes.
# In this demo, we are following the normalization procedures in [Numerical Study of Cosmic Ray Diffusion in MHD Turbulence](https://iopscience.iop.org/article/10.1088/0004-637X/728/1/60).
# The Lorentz equation for each particle of charge ``q`` and mass ``m``. The particle has a momentum ``\mathbf{p} = \gamma m \mathbf{v}`` and a normalized velocity ``\mathbf{v}_c = \mathbf{v} / c`` and propagates in an electromagnetic field ``\mathbf{E}`` (no mean electric field), ``\mathbf{B} = \delta \mathbf{B} + \mathbf{B}_0``:
# ```math
# \begin{aligned}
# \frac{\mathrm{d}\mathbf{p}}{\mathrm{d} t} &= q (\mathbf{E} + \mathbf{v}_c times \mathbf{B}) \\
# \frac{\mathrm{d}\mathbf{x}}{\mathrm{d} t} &= \mathbf{v}
# \end{aligned}
# ```
# Each particle is injected with a Lorentz factor ``\gamma_0``. Physically, one can think of ``gamma_0`` as a measure of the relativity of the particle, i.e., for small ``gamma_0`` we will recover nonrelativistic equations, and for large ``gamma_0`` --- ultra-relativistic equations. ``\gamma_0`` also defines the initial Larmor radius
# ```math
# r_{L0} = \gamma_0 m c^2 / (e B_0)
# ```
# where ``B_0`` is background magnetic field strength.
# After normalization, we have
# ```math
# \begin{aligned}
# \frac{\mathrm{d}\mathbf{v}^\prime}{\mathrm{d} t^\prime} &= \gamma^\prime \mathbf{E}^\prime + \mathbf{v}^\prime \times \mathbf{B}^\prime \\
# \frac{\mathrm{d}\mathbf{x}^\prime}{\mathrm{d} t^\prime} &= \mathbf{v}^\prime
# \end{aligned}
# ```
# where ``\mathbf{v}^\prime`` is the normalized space component of the 4-velocity ``(\gamma c, \gamma \mathbf{v}), \mathbf{v}^\prime = \mathbf{v} / \gamma_0``, and ``\mathbf{x}^\prime = \mathbf{x} / r_{L0}`` is the normalized location, ``\mathbf{E}^\prime = \mathbf{E} / (c B_0)`` is the normalized electric field, and ``\mathbf{B}^\prime = \mathbf{B} / B_0`` is the normalized magnetic field. ``t^\prime = e B_0 / (m c^2) t$ is self-time measured in cyclotron frequency units. A particle with pitch angle cosine ``mu = 0`` will make a full orbit in the ``B_0`` field in ``2 \pi`` time.
# ```math
# \gamma^\prime = \sqrt{\frac{1}{\gamma_0^2} + {\mathbf{u}^\prime}^2}
# ```
# Most of the time, ``1 / \gamma_0 \rightarrow 0`` since cosmic rays are of high energy.

# Let us first look at a simpler case without the electric field:
# ```math
# \begin{aligned}
# \frac{\mathrm{d}\mathbf{v}^\prime}{\mathrm{d} t^\prime} &= \mathbf{v}^\prime \times \mathbf{B}^\prime \\
# \frac{\mathrm{d}\mathbf{x}^\prime}{\mathrm{d} t^\prime} &= \mathbf{v}^\prime
# \end{aligned}
# ```
# This is simply the normalized Lorentz equation. By taking ``q=1, m=1``,
# ```math
# \begin{aligned}
# r_{L0} = \gamma_0 / B_0 \\
# \Omega_0 = 1 / r_{L0}
# \end{aligned}
# ```
# In the original MHD solution, the length scale is also dimensionless. There, the ratio between the initial Larmor radius ``r_{L0}`` and the simulation box scale length ``L`` becomes important, as this determines how many discrete points are involved within a gyroradius. The smaller ``r_{L0} / L`` is, the more likely a particle will experience an inhomogeneous magnetic field during gyration, and the more likely it will get scattered.

# Now let's demonstrate this with `trace_normalized!`.

import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type = "png") #hide

## Number of cells for the field along each dimension
nx, ny, nz = 4, 6, 2
## Unit conversion factors between dimensional and dimensionless units
γ0 = 10.0
B0 = 10.0
rL0 = γ0 / B0
L = rL0 * 4.0
Ω0 = 1 / rL0
## All quantities are in dimensionless units
x = range(-L/2-0.01, L/2+0.01, length=nx) # [Ω0]
y = range(-L-0.01, 0.01, length=ny) # [Ω0]
z = range(-10, 10, length=nz) # [Ω0]

B = fill(0.0, 3, nx, ny, nz) # [B0]
B[3,:,:,:] .= 1.0
E = fill(0.0, 3, nx, ny, nz) # [E₀]

param = prepare(x, y, z, E, B; species=User)

## Initial condition
stateinit = let
   x0 = [0.0, 0.0, 0.0] # initial position [l₀]
   u0 = [2.0, 0.0, 0.0] # initial velocity [v₀] -> r = 2 * rL0
   [x0..., u0...]
end
## Time span
tspan = (0.0, π) # half gyroperiod

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)
sol = solve(prob, Vern9())

### Visualization
f = Figure(fontsize = 18)
ax = Axis(f[1, 1],
   title = "Proton trajectory",
   xlabel = "X",
   ylabel = "Y",
   limits = (-2.1, 2.1, -4.1, 0.1),
   aspect = DataAspect()
)

lines!(ax, sol, idxs=(1,2))

xgrid = [i for i in x, _ in y]
ygrid = [j for _ in x, j in y]
scatter!(ax, xgrid[:], ygrid[:], color=:tomato)

f = DisplayAs.PNG(f) #hide


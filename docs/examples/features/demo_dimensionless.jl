# # Dimensionless Units and Normalization
#
# This example shows how to trace charged particles in dimensionless units.
#
# Fields from scale-free models such as MHD come in dimensionless units, so we can trace
# in those units and convert to any physical scale afterwards. Normalizing by the
# gyrofrequency also makes the grid spacing directly comparable to the gyroradius, which
# sets how many grid points a particle samples per gyration.
#
# For the cosmic-ray-specific form of this normalization — the pitch angle,
# energy-conserving tracing, and the ``r_{L0}/L`` resolution argument — see the
# [Cosmic Ray Propagation](@ref Cosmic-Ray) example.
#
# ## Normalization
#
# After normalization, ``q=1, B=1, m=1`` so that the gyroradius ``r_L = mv_\perp/qB = v_\perp``.
# All the quantities are given in dimensionless units: if the magnetic field is homogeneous and the initial perpendicular velocity ``v_{\perp 0}^\prime`` is 4, then the gyroradius is 4.
# To convert them to the original units, ``v_\perp = v_{\perp}^\prime U_0`` and ``r_L = r_L^\prime l_0 = 4 l_0``.
#
# The Lorentz equation in SI units is written as
#
# ```math
# \frac{\mathrm{d}\mathbf{v}}{\mathrm{d}t} = \frac{q}{m}\left( \mathbf{v}\times\mathbf{B} + \mathbf{E} \right)
# ```
#
# It can be normalized to
#
# ```math
# \frac{\mathrm{d}\mathbf{v}^\prime}{\mathrm{d}t^\prime} = \mathbf{v}^\prime\times\mathbf{B}^\prime + \mathbf{E}^\prime
# ```
#
# with the following transformation
#
# ```math
# \begin{aligned}
# \mathbf{v} &= \mathbf{v}^\prime U_0 \\
# t &= t^\prime t_0 = t^\prime \Omega^{-1} = t^\prime \frac{m}{qB_0} \\
# \mathbf{x} &= \mathbf{x}^\prime l_0 = \mathbf{x}^\prime U_0 t_0 \\
# \mathbf{B} &= \mathbf{B}^\prime B_0 \\
# \mathbf{E} &= \mathbf{E}^\prime E_0 = \mathbf{E}^\prime U_0 B_0
# \end{aligned}
# ```
#
# where all the coefficients with subscript 0 are expressed in SI units. All the variables with a prime are written in the dimensionless units.
#
# The normalized equation is the general one with ``q = m = 1``, so tracing in
# dimensionless units only needs `prepare(E, B; q = 1, m = 1)`.
#
# ## Dimensionless and Dimensional Tracing
#
# We first solve the Lorentz equation in SI units, and then convert the quantities to normalized units and solve it again in dimensionless units.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
using TestParticle: qᵢ, mᵢ, c
using CairoMakie
CairoMakie.activate!(type = "png") #hide

## Unit conversion factors between SI and dimensionless units
B_dim = 1.0e-8             # [T]
U_dim = c                # [m/s]
E_dim = U_dim * B_dim    # [V/m]
Ω_dim = abs(qᵢ) * B_dim / mᵢ # [1/s]
t_dim = 1 / Ω_dim        # [s]
l_dim = U_dim * t_dim    # [m]
## Electric field magnitude in SI units
Emag_dim = 1.0e-8           # [V/m]
### Solving in SI units
B_field(x) = SA[0, 0, B_dim]
E_field(x) = SA[Emag_dim, 0.0, 0.0]

## Initial conditions
x0_dim = [0.0, 0.0, 0.0] # [m]
v0_dim = [0.0, 0.01c, 0.0] # [m/s]
stateinit1 = [x0_dim..., v0_dim...]
tspan1 = (0, 2π * t_dim) # [s]

param1 = prepare(E_field, B_field, species = Proton)
prob1 = ODEProblem(trace!, stateinit1, tspan1, param1)
sol1 = solve(prob1, Vern9(); reltol = 1.0e-4, abstol = 1.0e-6)

### Solving in dimensionless units
B_normalize(x) = SA[0, 0, B_dim / B_dim]
E_normalize(x) = SA[Emag_dim / E_dim, 0.0, 0.0]
## For full EM problems, the normalization of E and B should be done separately.
param2 = prepare(E_normalize, B_normalize; m = 1, q = 1)
## Scale initial conditions by the conversion factors
x0_norm = x0_dim ./ l_dim
v0_norm = v0_dim ./ U_dim
tspan2 = (0, 2π)
stateinit2 = [x0_norm..., v0_norm...]

prob2 = ODEProblem(trace!, stateinit2, tspan2, param2)
sol2 = solve(prob2, Vern9(); reltol = 1.0e-4, abstol = 1.0e-6)

### Visualization
f = Figure(fontsize = 18)
ax = Axis(
    f[1, 1],
    xlabel = "x [km]",
    ylabel = "y [km]",
    aspect = DataAspect()
)

lines!(ax, sol1, idxs = (1, 2))
## Interpolate dimensionless solutions and map back to SI units
xp, yp = let trange = range(tspan2..., length = 40)
    sol2.(trange, idxs = 1) .* l_dim, sol2.(trange, idxs = 2) .* l_dim
end
lines!(ax, xp, yp, linestyle = :dashdot, linewidth = 5, color = Makie.wong_colors()[2])
invL = inv(1.0e3)
scale!(ax.scene, invL, invL)

f = DisplayAs.PNG(f) #hide

# We see that the results are almost identical, with only floating point numerical errors. Tracing in dimensionless units usually allows larger timesteps, which leads to faster computation.

# ## Relativistic Tracing
#
# In the relativistic case,
# * Velocity is normalized by speed of light c, $V = V^\prime c$;
# * Magnetic field is normalized by the characteristic magnetic field, $B = B^\prime B_0$;
# * Electric field is normalized by $E_0 = c\,B_0$, $E = E^\prime E_0$;
# * Location is normalized by the $L = c / \Omega_0$, where $\Omega_0 = q\,B_0 / m$, and
# * Time is normalized by $\Omega_0^{-1}$, $t = t^\prime * \Omega_0^{-1}$.
#
# Hence `c = 1` in `prepare`; `relativistic = true` makes the state `(x, γv)`.
# At ``v = 0.9`` the orbit closes after ``2\pi\gamma`` with gyroradius ``\gamma v \approx 2.06``,
# instead of after ``2\pi`` with radius 0.9 as in the non-relativistic case.

param = prepare(xu -> SA[0.0, 0.0, 0.0], xu -> SA[0.0, 0.0, 1.0]; m = 1, q = 1, c = 1)
v = 0.9
γ = 1 / sqrt(1 - v^2)
stateinit = [0.0, 0.0, 0.0, γ * v, 0.0, 0.0]
tspan = (0.0, 2π * γ) # one gyroperiod
prob = TraceProblem(stateinit, tspan, param; relativistic = true)
sol = solve(prob, Vern9())

### Visualization
f = Figure(fontsize = 18)
ax = Axis(
    f[1, 1],
    title = "Relativistic particle trajectory",
    xlabel = "X",
    ylabel = "Y",
    aspect = DataAspect()
)

lines!(ax, sol, idxs = (1, 2))

f = DisplayAs.PNG(f) #hide

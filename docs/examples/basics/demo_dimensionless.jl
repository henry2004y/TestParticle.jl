# ---
# title: Dimensionless Units
# id: demo_dimensionless
# date: 2023-12-14
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.4
# description: Tracing charged particle in dimensionless units
# ---

# This example shows how to trace charged particles in dimensionless units.
# After normalization, ``q=1, B=1, m=1`` so that the gyroradius ``r_L = mv_\perp/qB = v_\perp``. All the quantities are given in dimensionless units.
# If the magnetic field is homogeneous and the initial perpendicular velocity is 4, then the gyroradius is 4.
# To convert them to the original units, ``v_\perp = 4*U_0`` and ``r_L = 4*l_0``.
# Check [Demo: single tracing with additional diagnostics](@ref demo_savingcallback) for explaining the unit conversion.

# Tracing in dimensionless units is beneficial for many scenarios. For example, MHD simulations do not have intrinsic scales. Therefore, we can do dimensionless particle tracing in MHD fields, and then convert to any scale we would like.

# Now let's demonstrate this with `trace_normalized!`.

import DisplayAs #hide
using TestParticle
using TestParticle: qᵢ, mᵢ
using OrdinaryDiffEq

using CairoMakie
CairoMakie.activate!(type = "png") #hide

## Number of cells for the field along each dimension
nx, ny, nz = 4, 6, 8
## Unit conversion factors between SI and dimensionless units
B₀ = 10e-9            # [T]
Ω = abs(qᵢ) * B₀ / mᵢ # [1/s]
t₀ = 1 / Ω            # [s]
U₀ = 1.0              # [m/s]
l₀ = U₀ * t₀          # [m]
E₀ = U₀*B₀            # [V/m]
## All quantities are in dimensionless units
x = range(-10, 10, length=nx) # [l₀]
y = range(-10, 10, length=ny) # [l₀]
z = range(-10, 10, length=nz) # [l₀]

B = fill(0.0, 3, nx, ny, nz) # [B₀]
B[3,:,:,:] .= 1.0
E = fill(0.0, 3, nx, ny, nz) # [E₀]

param = prepare(x, y, z, E, B; species=User)

x0 = [0.0, 0.0, 0.0] # initial position [l₀]
u0 = [4.0, 0.0, 0.0] # initial velocity [v₀]
stateinit = [x0..., u0...]
tspan = (0.0, π) # half gyroperiod

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)
sol = solve(prob, Vern9())

### Visualization
f = Figure(fontsize = 18)
ax = Axis(f[1, 1],
   title = "Proton trajectory",
   xlabel = "X",
   ylabel = "Y",
   limits = (-4.1, 4.1, -8.1, 0.1),
   aspect = DataAspect()
)

lines!(ax, sol, vars=(1,2))

f = DisplayAs.PNG(f) #hide
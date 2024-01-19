# ---
# title: Tracing with Dimensionless Units and Periodic Boundary
# id: demo_dimensionless_periodic
# date: 2023-12-15
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.4
# description: Tracing charged particle in dimensionless units and periodic boundary
# ---

# This example shows how to trace charged particles in dimensionless units and EM fields with periodic boundaries.
# For details about dimensionless units, please check another [example](@ref demo_dimensionless).

# Now let's demonstrate this with `trace_normalized!`.

import DisplayAs #hide
using TestParticle
using TestParticle: qᵢ, mᵢ
using OrdinaryDiffEq
using Meshes
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png")

x = range(-10, 10, length=15)
y = range(-10, 10, length=20)
B = fill(0.0, 3, length(x), length(y)) # [B₀]

B[3,:,:] .= 10e-9
## This is the magnitude of B
B₀ = 10e-9 # [nT]

Ω = abs(qᵢ) * B₀ / mᵢ
t₀ = 1 / Ω  # [s]
U₀ = 1.0    # [m/s]
l₀ = U₀ * t₀ # [m]
E₀ = U₀*B₀ # [V/m]

## For full EM problems, the normalization of E and B should be done separately.
B ./= B₀
E(x) = SA[0.0/E₀, 0.0/E₀, 0.0/E₀]

## If bc == 1, we set a NaN value outside the domain (default);
## If bc == 2, we set periodic boundary conditions.
param = prepare(E, B; species=User, bc=2)

# Note that we set a radius of 10, so the trajectory extent from -20 to 0 in y, which is beyond the original y range.

x0 = [0.0, 0.0, 0.0] # initial position [l₀]
u0 = [10, 0.0, 0.0] # initial velocity [v₀]
stateinit = [x0..., u0...]
tspan = (0.0, 1.5π) # 3/4 gyroperiod

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)
sol = solve(prob, Vern9())

### Visualization
f = Figure(fontsize = 18)
ax = Axis(f[1, 1],
   title = "Proton trajectory",
   xlabel = "X",
   ylabel = "Y",
   limits = (-10.1, 10.1, -20.1, 0.1),
   aspect = DataAspect()
)

lines!(ax, sol, vars=(1,2))

f = DisplayAs.PNG(f) #hide
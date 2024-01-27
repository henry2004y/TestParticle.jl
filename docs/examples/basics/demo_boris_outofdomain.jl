# ---
# title: Boris tracing with domain check
# id: demo_boris_domain
# date: 2024-01-27
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.0
# description: Boris tracing with field out of domain check
# ---

# This example shows how to trace charged particles using the Boris method in dimensionless units with additionally boundary check.
# If the particle travels out of the domain specified by the field, the tracing will stop.
# Check [Demo: Dimensionless Units](@ref demo_dimensionless) for explaining the unit conversion, and [Demo: Boris Method](@ref demo_boris) for introducing the Boris method.

import DisplayAs #hide
using TestParticle
using TestParticle: qᵢ, mᵢ
using StaticArrays
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type = "png")

uniform_B(x) = SA[0.0, 0.0, 0.01]
uniform_E(x) = SA[0.0, 0.0, 0.0]

function isoutofdomain(xv)
   if isnan(xv[1])
      return true
   else
      return false
   end
end

## Number of cells for the field along each dimension
nx, ny = 4, 6
## Unit conversion factors between SI and dimensionless units
B₀ = 10e-9            # [T]
Ω = abs(qᵢ) * B₀ / mᵢ # [1/s]
t₀ = 1 / Ω            # [s]
U₀ = 1.0              # [m/s]
l₀ = U₀ * t₀          # [m]
E₀ = U₀*B₀            # [V/m]

x = range(0, 11, length=nx) # [l₀]
y = range(-21, 0, length=ny) # [l₀]

B = fill(0.0, 3, nx, ny) # [B₀]
B[3,:,:] .= 1.0

E(x) = SA[0.0, 0.0, 0.0] # [E₀]

## If bc == 1, we set a NaN value outside the domain (default);
## If bc == 2, we set periodic boundary conditions.
param = prepare(x, y, E, B; species=User, bc=1);

# Note that we set a radius of 10, so the trajectory extent from -20 to 0 in y, and -10 to 10 in x.
# After half a cycle, the particle will move into the region where is field is not defined.
# The tracing will stop with the final step being all NaNs.

x0 = [0.0, 0.0, 0.0] # initial position [l₀]
u0 = [10.0, 0.0, 0.0] # initial velocity [v₀]
stateinit = [x0..., u0...]
tspan = (0.0, 1.5π) # 3/4 gyroperiod

paramBoris = BorisMethod(param)
dt = 0.1
prob = TraceProblem(stateinit, tspan, dt, paramBoris)

sol = trace_trajectory(prob; savestepinterval=1, isoutofdomain)

f = Figure(fontsize = 18)
ax = Axis(f[1, 1],
   title = "Proton trajectory",
   xlabel = "X",
   ylabel = "Y",
   limits = (-10.1, 10.1, -20.1, 0.1),
   aspect = DataAspect()
)

@views lines!(ax, sol[1].u[1,:], sol[1].u[2,:])

f = DisplayAs.PNG(f) #hide
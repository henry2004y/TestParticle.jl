import DisplayAs #hide
using TestParticle
using TestParticle: qᵢ, mᵢ
using OrdinaryDiffEq, StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# Number of cells for the field along each dimension
nx, ny, nz = 4, 6, 8
# Unit conversion factors between SI and dimensionless units
B₀ = 10e-9            # [T]
Ω = abs(qᵢ) * B₀ / mᵢ # [1/s]
t₀ = 1 / Ω            # [s]
U₀ = 1.0              # [m/s]
l₀ = U₀ * t₀          # [m]
E₀ = U₀*B₀            # [V/m]
# All quantities are in dimensionless units
x = range(-10, 10, length=nx) # [l₀]
y = range(-10, 10, length=ny) # [l₀]
z = range(-10, 10, length=nz) # [l₀]

B = fill(0.0, 3, nx, ny, nz) # [B₀]
B[3,:,:,:] .= 1.0
E = fill(0.0, 3, nx, ny, nz) # [E₀]

param = prepare(x, y, z, E, B; species=User)

# Initial condition
stateinit = let
   x0 = [0.0, 0.0, 0.0] # initial position [l₀]
   u0 = [4.0, 0.0, 0.0] # initial velocity [v₀]
   [x0..., u0...]
end
# Time span
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

lines!(ax, sol, idxs=(1,2))

f = DisplayAs.PNG(f) #hide

param = prepare(xu -> SA[0.0, 0.0, 0.0], xu -> SA[0.0, 0.0, 1.0]; species=User)
tspan = (0.0, π) # half period
stateinit = [0.0, 0.0, 0.0, 0.01, 0.0, 0.0]
prob = ODEProblem(trace_relativistic_normalized!, stateinit, tspan, param)
sol = solve(prob, Vern9())

### Visualization
f = Figure(fontsize = 18)
ax = Axis(f[1, 1],
   title = "Relativistic particle trajectory",
   xlabel = "X",
   ylabel = "Y",
   #limits = (-0.6, 0.6, -1.1, 0.1),
   aspect = DataAspect()
)

lines!(ax, sol, idxs=(1,2))

f = DisplayAs.PNG(f) #hide

param = prepare(xu -> SA[0.0, 0.0, 0.0], xu -> SA[0.0, 0.0, 1.0]; species=User)
tspan = (0.0, π) # half period
stateinit = [0.0, 0.0, 0.0, 0.9, 0.0, 0.0]
prob = ODEProblem(trace_relativistic_normalized!, stateinit, tspan, param)
sol = solve(prob, Vern9())

### Visualization
f = Figure(fontsize = 18)
ax = Axis(f[1, 1],
   title = "Relativistic particle trajectory",
   xlabel = "X",
   ylabel = "Y",
   aspect = DataAspect()
)

lines!(ax, sol, idxs=(1,2))

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

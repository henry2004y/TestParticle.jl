import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# Number of cells for the field along each dimension
nx, ny, nz = 4, 6, 2
# Unit conversion factors between dimensional and dimensionless units
rL0 = 1.0
L = rL0 * 4.0
Ω0 = 1 / rL0
# All quantities are in dimensionless units
x = range(-L/2-0.01, L/2+0.01, length=nx) # [rL0]
y = range(-L-0.01, 0.01, length=ny) # [rL0]
z = range(-10, 10, length=nz) # [rL0]

B = fill(0.0, 3, nx, ny, nz) # [B0]
B[3,:,:,:] .= 1.0
E = fill(0.0, 3, nx, ny, nz) # [E₀]

param = prepare(x, y, z, E, B; species=User)

# Initial condition
stateinit = let
   x0 = [0.0, 0.0, 0.0] # initial position [l₀]
   u0 = [1.0, 0.0, 0.0] # initial velocity [v₀] -> r = 1 * rL0
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
   limits = (-2.1, 2.1, -4.1, 0.1),
   aspect = DataAspect()
)

lines!(ax, sol, idxs=(1,2))

xgrid = [i for i in x, _ in y]
ygrid = [j for _ in x, j in y]
scatter!(ax, xgrid[:], ygrid[:], color=:tomato)

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using Meshes

using CairoMakie
CairoMakie.activate!(type = "png")

x = range(-10, 10, length=15)
y = range(-10, 10, length=20)
z = range(-10, 10, length=25)
B = fill(0.0, 3, length(x), length(y), length(z)) # [B₀]
E = fill(0.0, 3, length(x), length(y), length(z)) # [E₀]

B₀ = 10e-9

B[3,:,:,:] .= 1.0

Δx = x[2] - x[1]
Δy = y[2] - y[1]
Δz = z[2] - z[1]

mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
   (x[1], y[1], z[1]),
   (Δx, Δy, Δz))

param = prepare(mesh, E, B, B₀; species=Proton)

Ω = param[1]
U₀ = 1.0

x0 = [0.0, 0.0, 0.0] # initial position [l₀]
u0 = [4*Ω*U₀, 0.0, 0.0] # initial velocity [v₀]
stateinit = [x0..., u0...]
tspan = (0.0, 2π/Ω) # one gyroperiod

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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

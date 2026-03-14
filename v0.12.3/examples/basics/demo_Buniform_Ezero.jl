import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

uniform_B(x) = SA[0.0, 0.0, 1e-8]
uniform_E(x) = SA[0.0, 0.0, 0.0]

# Initial condition
stateinit = let x0 = [1.0, 0, 0], v0 = [0.0, 1.0, 0.1]
   [x0..., v0...]
end
# Time span
tspan = (0, 18)

param = prepare(uniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())

# Visualization
f = Figure(fontsize=18)
ax = Axis3(f[1, 1],
   title = "Helix Trajectory",
   xlabel = "x [m]",
   ylabel = "y [m]",
   zlabel = "z [m]",
   aspect = :data,
)

plot!(ax, sol, idxs=(1, 2, 3))

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

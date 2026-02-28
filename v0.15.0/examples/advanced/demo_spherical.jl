import DisplayAs #hide
import TestParticle
using OrdinaryDiffEq, StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# Define the magnetic field in spherical coordinates
r = logrange(0.1, 10.0, length=4)
θ = range(0, π, length=8)
ϕ = range(0, 2π, length=8)

B₀ = 1e-8
B = zeros(3, length(r), length(θ), length(ϕ))

for (iθ, θ_val) in enumerate(θ)
    B[1, :, iθ, :] .= B₀ * cos(θ_val)
    B[2, :, iθ, :] .= -B₀ * sin(θ_val)
end

B_field = TestParticle.getinterp(TestParticle.SphericalNonUniformR(), B, r, θ, ϕ)
zero_E = TestParticle.ZeroField()

# Initial condition
stateinit = let x0 = [3.0, 0.0, -1.0], v0 = [0.0, 1.0, 0.1]
   [x0..., v0...]
end
# Time span
tspan = (0, 18)

param = TestParticle.prepare(zero_E, B_field, species=TestParticle.Proton)
prob = ODEProblem(TestParticle.trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())

# Visualization
f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Helix Trajectory in Spherical B-field",
   xlabel = "x [m]",
   ylabel = "y [m]",
   zlabel = "z [m]",
   aspect = :data
)

plot!(ax, sol, idxs = (1, 2, 3))

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

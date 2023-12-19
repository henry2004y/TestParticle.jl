import DisplayAs #hide
using TestParticle
using TestParticle: c, qᵢ, mᵢ
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png")

# Unit conversion factors
const B₀ = 1e-8 # [T]
const U₀ = c    # [m/s]
const E₀ = U₀*B₀ # [V/m]
const Ω = abs(qᵢ) * B₀ / mᵢ
const t₀ = 1 / Ω  # [s]
const l₀ = U₀ * t₀ # [m]

const Emag = 1e-8 # [V/m]
### Solving in SI units
B(x) = SA[0, 0, B₀]
E(x) = SA[Emag, 0.0, 0.0]

x0 = [0.0, 0.0, 0.0] # [m]
v0 = [0.0, 0.01c, 0.0] # [m/s]
stateinit1 = [x0..., v0...]
tspan1 = (0, 2π/Ω) # [s]

param1 = prepare(E, B, species=Proton)
prob1 = ODEProblem(trace!, stateinit1, tspan1, param1)
sol1 = solve(prob1, Vern9(); reltol=1e-4, abstol=1e-6)

### Solving in dimensionless units
B_normalize(x) = SA[0, 0, 1.0]
E_normalize(x) = SA[Emag/E₀, 0.0, 0.0]
# Default User type has q=1, m=1. Here we also set B=0.
param2 = prepare(E_normalize, B_normalize, 1.0; species=User)
# Scale initial quantities by the conversion factors
x0 ./= l₀
v0 ./= U₀
tspan2 = (0, 2π/Ω/t₀)
stateinit2 = [x0..., v0...]

prob2 = ODEProblem(trace_normalized!, stateinit2, tspan2, param2)
sol2 = solve(prob2, Vern9(); reltol=1e-4, abstol=1e-6)

### Visualization
f = Figure(fontsize=18)
ax = Axis3(f[1, 1],
    xlabel = "x",
    ylabel = "y",
    zlabel = "z",
    aspect = :data,
)

lines!(ax, sol1)
# Interpolate dimensionless solutions and map back to SI units
xp, yp, zp = let trange = range(tspan2..., length=40)
   sol2.(trange, idxs=1) .* l₀, sol2.(trange, idxs=2) .* l₀, sol2.(trange, idxs=3) .* l₀
end
lines!(ax, xp, yp, zp, linestyle=:dashdot, linewidth=5)

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

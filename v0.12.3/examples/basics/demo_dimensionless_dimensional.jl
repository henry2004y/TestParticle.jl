import DisplayAs #hide
using TestParticle
using TestParticle: c, qᵢ, mᵢ
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# Unit conversion factors between SI and dimensionless units
const B₀ = 1e-8             # [T]
const U₀ = c                # [m/s]
const E₀ = U₀*B₀            # [V/m]
const Ω = abs(qᵢ) * B₀ / mᵢ # [1/s]
const t₀ = 1 / Ω            # [s]
const l₀ = U₀ * t₀          # [m]
# Electric field magnitude in SI units
const Emag = 1e-8           # [V/m]
### Solving in SI units
B(x) = SA[0, 0, B₀]
E(x) = SA[Emag, 0.0, 0.0]

# Initial conditions
x0 = [0.0, 0.0, 0.0] # [m]
v0 = [0.0, 0.01c, 0.0] # [m/s]
stateinit1 = [x0..., v0...]
tspan1 = (0, 2π*t₀) # [s]

param1 = prepare(E, B, species=Proton)
prob1 = ODEProblem(trace!, stateinit1, tspan1, param1)
sol1 = solve(prob1, Vern9(); reltol=1e-4, abstol=1e-6)

### Solving in dimensionless units
B_normalize(x) = SA[0, 0, B₀/B₀]
E_normalize(x) = SA[Emag/E₀, 0.0, 0.0]
# For full EM problems, the normalization of E and B should be done separately.
param2 = prepare(E_normalize, B_normalize; species=User)
# Scale initial conditions by the conversion factors
x0 ./= l₀
v0 ./= U₀
tspan2 = (0, 2π)
stateinit2 = [x0..., v0...]

prob2 = ODEProblem(trace_normalized!, stateinit2, tspan2, param2)
sol2 = solve(prob2, Vern9(); reltol=1e-4, abstol=1e-6)

### Visualization
f = Figure(fontsize=18)
ax = Axis(f[1, 1],
    xlabel = "x [km]",
    ylabel = "y [km]",
    aspect = DataAspect(),
)

lines!(ax, sol1, idxs=(1, 2))
# Interpolate dimensionless solutions and map back to SI units
xp, yp = let trange = range(tspan2..., length=40)
   sol2.(trange, idxs=1) .* l₀, sol2.(trange, idxs=2) .* l₀
end
lines!(ax, xp, yp, linestyle=:dashdot, linewidth=5, color=Makie.wong_colors()[2])
invL = inv(1e3)
scale!(ax.scene, invL, invL)

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

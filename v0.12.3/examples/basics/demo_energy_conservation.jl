import DisplayAs #hide
using TestParticle
import TestParticle as TP
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra: ×
using CairoMakie
CairoMakie.activate!(type = "png") #hide

const B₀ = 1e-8 # [T]
const E₀ = 3e-2 # [V/m]

"f2"
function location!(dx, v, x, p::TP.TPTuple, t)
   dx .= v
end

"f1"
function lorentz!(dv, v, x, p::TP.TPTuple, t)
   q2m, E, B = p
   dv .= q2m*(E(x, t) + v × (B(x, t)))
end

### Initialize field

function uniform_B(x)
   return SA[0, 0, B₀]
end

function uniform_E(x)
   return SA[E₀, 0.0, 0.0]
end

function zero_B(x)
   return SA[0.0, 0.0, 0.0]
end

function zero_E(x)
   return SA[0.0, 0.0, 0.0]
end

"Check energy conservation."
E(dx, dy, dz) = 1 // 2 * (dx^2 + dy^2 + dz^2)

### Initialize particles

x0 = [0.0, 0, 0]
v0 = [0.0, 1e2, 0.0]
stateinit = [x0..., v0...]
tspan_proton = (0.0, 2000.0);

param_proton = prepare(zero_E, uniform_B, species=Proton)

### Solve for the trajectories

prob_p = DynamicalODEProblem(lorentz!, location!, v0, x0, tspan_proton, param_proton)

Ωᵢ = TP.qᵢ * B₀ / TP.mᵢ
Tᵢ = 2π / Ωᵢ
println("Number of gyrations: ", tspan_proton[2] / Tᵢ)

sol = solve(prob_p, ImplicitMidpoint(), dt = Tᵢ/15)

f = Figure(fontsize=18)
ax = Axis(f[1, 1],
    title = "Proton in a uniform B field and zero E field",
    xlabel = "x",
    ylabel = "y",
    aspect = 1,
)

lines!(ax, sol, idxs=(1, 2))

f = DisplayAs.PNG(f) #hide

param_proton = prepare(uniform_E, zero_B, species=Proton)

# acceleration, [m/s²]
a = param_proton[1] * E₀
# predicted final speed, [m/s]
v_final_predict = a * tspan_proton[2]
# predicted travel distance, [m/s]
d_final_predict = 0.5 * tspan_proton[2] * v_final_predict
# predicted energy gain, [eV]
E_predict = E₀ * d_final_predict

prob_p = DynamicalODEProblem(lorentz!, location!, v0, x0, tspan_proton, param_proton)

sol = solve(prob_p, Vern6())

energy = map(x -> E(x[1:3]...), sol.u) .* TP.mᵢ;

println("predicted final speed: $v_final_predict [m/s]") #hide

println("simulated final speed: $(sol.u[end][1]) [m/s]") #hide

println("predicted travel distance: $d_final_predict [m]") #hide

println("simulated travel distance: $(sol.u[end][4]) [m]") #hide

println("predicted energy gain: $E_predict [eV]") #hide

println("simulated final energy: $(energy[end] / TP.qᵢ) [eV]") #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

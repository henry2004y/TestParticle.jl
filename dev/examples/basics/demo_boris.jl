import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
using TestParticle: ZeroField, get_BField
import TestParticle as TP
using CairoMakie
CairoMakie.activate!(type = "png") #hide

function plot_trajectory(sol_boris, sol1, sol2)
   f = Figure(size = (700, 600), fontsize = 18)
   ax = Axis(f[1, 1], aspect = 1, limits = (-3, 1, -2, 2),
      xlabel = "X",
      ylabel = "Y")
   idxs = (1, 2)
   l0 = lines!(ax, sol_boris[1]; idxs, linewidth = 2, label = "Boris")
   l1 = lines!(ax, sol1; idxs,
      color = Makie.wong_colors()[2], linewidth = 2, linestyle = :dashdot, label = "Tsit5 fixed")
   l2 = linesegments!(ax, sol2; idxs,
      color = Makie.wong_colors()[3], linewidth = 2, linestyle = :dot, label = "Tsit5 adaptive")

   scale!(ax.scene, invrL, invrL)

   axislegend(position = :rt, framevisible = false)

   f
end

const Bmag = 0.01
uniform_B(x) = SA[0.0, 0.0, Bmag]
zero_E = ZeroField()

x0 = [0.0, 0.0, 0.0]
v0 = [0.0, 1e5, 0.0]
stateinit = [x0..., v0...]

param = prepare(zero_E, uniform_B, species = Electron)

# Reference parameters
const tperiod = 2π / (abs(param[1]) *
                 sqrt(sum(x -> x^2, get_BField(param)([0.0, 0.0, 0.0], 0.0))))
const rL = sqrt(v0[1]^2 + v0[2]^2 + v0[3]^2) / (abs(param[1]) * Bmag)
const invrL = 1 / rL;

tspan = (0.0, tperiod)
dt = tperiod / 4

prob = TraceProblem(stateinit, tspan, param)

sol_boris = TP.solve(prob; dt, savestepinterval = 1);

prob = ODEProblem(trace!, stateinit, tspan, param)
sol1 = solve(prob, Tsit5(); adaptive = false, dt, dense = false, saveat = dt);
sol2 = solve(prob, Tsit5());

# Visualization
f = plot_trajectory(sol_boris, sol1, sol2)
f = DisplayAs.PNG(f) #hide

dt = tperiod / 8

prob = TraceProblem(stateinit, tspan, param)

sol_boris = TP.solve(prob; dt, savestepinterval = 1);

prob = ODEProblem(trace!, stateinit, tspan, param)
sol1 = solve(prob, Tsit5(); adaptive = false, dt, dense = false, saveat = dt);

# Visualization
f = plot_trajectory(sol_boris, sol1, sol2)
f = DisplayAs.PNG(f) #hide

tspan = (0.0, 200*tperiod)
dt = tperiod / 12

prob_boris = TraceProblem(stateinit, tspan, param)
prob = ODEProblem(trace!, stateinit, tspan, param)

sol_boris = TP.solve(prob_boris; dt, savestepinterval = 10);
sol1 = solve(prob, Tsit5(); adaptive = false, dt, dense = false, saveat = dt);
sol2 = solve(prob, Tsit5());
sol3 = solve(prob, Vern7());
sol4 = solve(prob, Vern9());

# Visualization
f = plot_trajectory(sol_boris, sol1, sol2)
f = DisplayAs.PNG(f) #hide

E(vx, vy, vz) = 1 // 2 * (vx^2 + vy^2 + vz^2)
f = Figure(size = (800, 400), fontsize = 18)
ax = Axis(f[1, 1],
   xlabel = "time [period]",
   ylabel = "Normalized Kinetic Energy")

sols_to_plot = [
   (sol_boris[1], "Boris"),
   (sol1, "Tsit5 fixed"),
   (sol2, "Tsit5 adaptive"),
   (sol3, "Vern7 adaptive"),
   (sol4, "Vern9 adaptive")
]

for (sol, label) in sols_to_plot
   energy = map(x -> E(x[4:6]...), sol.u)
   lines!(ax, sol.t ./ tperiod, energy ./ energy[1], label=label)
end

f = DisplayAs.PNG(f) #hide

@time sol_boris = TP.solve(prob_boris; dt, savestepinterval = 10)[1];
@time sol1 = solve(prob, Tsit5(); adaptive = false, dt, dense = false, saveat = dt);
@time sol2 = solve(prob, Tsit5());
@time sol3 = solve(prob, Vern7());
@time sol4 = solve(prob, Vern9());

t = tspan[2] / 2
sol_boris(t)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

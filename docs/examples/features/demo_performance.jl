# # Performance Benchmark
#
# This example compares the performance of the native Boris method against various ODE solvers from `OrdinaryDiffEq.jl`.
# We use [Chairmarks.jl](https://github.com/LilithHafner/Chairmarks.jl) for benchmarking.

using Chairmarks
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
using Statistics
using Random
import DisplayAs #hide
import TestParticle as TP

CairoMakie.activate!(type = "png") #hide

# ## Simulation Setup
#
# We use a simple setup: an electron in a uniform magnetic field.

const Bmag = 0.01
uniform_B(x) = SA[0.0, 0.0, Bmag]
zero_E = TP.ZeroField()

x0 = SA[0.0, 0.0, 0.0]
v0 = SA[0.0, 1e5, 0.0]
stateinit = SA[x0..., v0...]
## (q2m, m, E, B, F)
param = prepare(zero_E, uniform_B, species = Electron)
q2m = TP.get_q2m(param)

## Reference parameters
const tperiod = 2π / (abs(q2m) * Bmag)

tspan = (0.0, 200 * tperiod)
dt = tperiod / 12

## Native Boris method requires a mutable state (Vector or MVector)
prob_boris = TraceProblem(MVector(stateinit), tspan, param)
## ODE solvers from DifferentialEquations.jl are optimized for StaticArrays (SVector)
prob_ode = ODEProblem(trace, stateinit, tspan, param)

# ## Benchmark
#
# We benchmark the following solvers:
#
# | Solver | Description |
# | :--- | :--- |
# | Boris (n=1) | Native Boris with n=1 |
# | Boris (n=2) | Native Multistep Boris with n=2 |
# | Tsit5 (fixed) | `OrdinaryDiffEq` Tsit5 with fixed step |
# | Tsit5 (adaptive) | `OrdinaryDiffEq` Tsit5 with adaptive step |
# | Vern7 (fixed) | `OrdinaryDiffEq` Vern7 with fixed step |
# | Vern7 (adaptive) | `OrdinaryDiffEq` Vern7 with adaptive step |
# | Vern9 (fixed) | `OrdinaryDiffEq` Vern9 with fixed step |
# | Vern9 (adaptive) | `OrdinaryDiffEq` Vern9 with adaptive step |
#
# To simulate realistic applications, we save the solution at fixed intervals for all solvers.

"""
Helper functions to extract the median execution time and memory allocation.
"""
function get_median_time_memory(b)
   mb = median(b)
   return mb.time, mb.bytes
end

solvers = [
   ("Boris (n=1)", () -> TP.solve(prob_boris; dt)),
   ("Boris (n=2)", () -> TP.solve(prob_boris; dt, n = 2)),
   ("Tsit5 (fixed)", () -> solve(prob_ode, Tsit5(); adaptive = false, dt, dense = false)),
   ("Tsit5 (adaptive)", () -> solve(prob_ode, Tsit5(); saveat = dt)),
   ("Vern7 (fixed)", () -> solve(prob_ode, Vern7(); adaptive = false, dt, dense = false)),
   ("Vern7 (adaptive)", () -> solve(prob_ode, Vern7(); saveat = dt)),
   ("Vern9 (fixed)", () -> solve(prob_ode, Vern9(); adaptive = false, dt, dense = false)),
   ("Vern9 (adaptive)", () -> solve(prob_ode, Vern9(); saveat = dt))
]

n_solvers = length(solvers)
results_time = Vector{Float64}(undef, n_solvers)
results_mem = Vector{Float64}(undef, n_solvers)
names = Vector{String}(undef, n_solvers)

for (i, (name, func)) in enumerate(solvers)
   println("Benchmarking $name...")
   b = @be $func() seconds=1
   mt, mm = get_median_time_memory(b)
   results_time[i] = mt
   results_mem[i] = mm
   names[i] = name
end

# Normalize results
min_time = minimum(results_time)
min_mem = minimum(results_mem)

results_time_norm = results_time ./ min_time
results_mem_norm = results_mem ./ min_mem;

# ## Visualization

f = Figure(size = (1200, 1000), fontsize = 24)

ax = Axis(f[1, 1],
   title = "Solver Efficiency (Time vs. Memory)",
   xlabel = "Relative Time (1.0 = Fastest)",
   ylabel = "Relative Memory (1.0 = Lowest)",
   xgridstyle = :dash,
   ygridstyle = :dash,
   xminorticksvisible = true,
   yminorticksvisible = true,
   xminorticks = IntervalsBetween(5),
   yminorticks = IntervalsBetween(5)
)

sc = scatter!(ax, results_time_norm, results_mem_norm,
   color = 1:n_solvers,
   colormap = :tab10,
   markersize = 30,
   strokewidth = 1,
   strokecolor = :black
)

## Add annotations with random offsets to fix overlaps
rng = Random.MersenneTwister(42)
offsets = [(10, rand(rng, -30:30)) for _ in 1:n_solvers]

text!(ax, results_time_norm, results_mem_norm,
   text = names,
   align = (:left, :center),
   offset = offsets,
   fontsize = 24
)

## Highlight the "Utopia Point" (Theoretical Best)
scatter!(ax, [1.0], [1.0],
   marker = :star5,
   markersize = 20,
   color = :red,
   label = "Ideal Limit"
)
text!(ax, 1.0, 1.0, text = "Utopia Point", align = (:right, :top),
   offset = (55, -5), color = :red, fontsize = 15)

## Add "Iso-Efficiency" curves (Optional visual aid)
## Curves where Time * Memory = Constant (Cost invariant)
x_range = range(
   minimum(results_time_norm) * 0.8, stop = maximum(results_time_norm) * 1.1, length = 100)
lines!(ax, x_range, 5.0 ./ x_range, color = (:gray, 1.0), linestyle = :dot)
text!(ax, maximum(x_range), 5.0 / maximum(x_range),
   text = "Iso-cost", fontsize = 20, color = :black)

xlims!(ax, minimum(results_time_norm) * 0.9, maximum(results_time_norm) * 1.2)
ylims!(ax, minimum(results_mem_norm) * 0.9, maximum(results_mem_norm) * 1.1)

f = DisplayAs.PNG(f) #hide

# In practice, it is pretty hard to find an optimal algorithm. The native Boris method is good if you want a fixed time step.
# When calling OrdinaryDiffEq.jl, we recommend using `Vern9()` as a starting point instead of `Tsit5()`, especially combined with adaptive timestepping. Further fine-grained control includes setting `dtmax`, `reltol`, and `abstol` in the `solve` method.

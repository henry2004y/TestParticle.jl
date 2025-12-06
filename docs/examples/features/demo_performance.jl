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
import DisplayAs
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
# We define helper functions to extract the median execution time and memory allocation.

function get_median_time_memory(b)
    mb = median(b)
    return mb.time, mb.bytes
end

# We benchmark the following solvers:
# 1. Native Boris (n=1)
# 2. Native Multistep Boris (n=2)
# 3. Tsit5 (fixed step)
# 4. Tsit5 (adaptive)
# 5. Vern7 (fixed step)
# 6. Vern7 (adaptive)
# 7. Vern9 (fixed step)
# 8. Vern9 (adaptive)

# We use a reduced sampling time/output size to focus on computation speed.
# For a fair comparison, we only save the start and end states for all solvers.

solvers = [
    ("Boris (n=1)", () -> TP.solve(prob_boris; dt, save_everystep = false)),
    ("Boris (n=2)", () -> TP.solve(prob_boris; dt, save_everystep = false, n = 2)),
    ("Tsit5 (fixed)", () -> solve(prob_ode, Tsit5(); adaptive = false, dt, dense = false, save_everystep = false)),
    ("Tsit5 (adaptive)", () -> solve(prob_ode, Tsit5(); save_everystep = false)),
    ("Vern7 (fixed)", () -> solve(prob_ode, Vern7(); adaptive = false, dt, dense = false, save_everystep = false)),
    ("Vern7 (adaptive)", () -> solve(prob_ode, Vern7(); save_everystep = false)),
    ("Vern9 (fixed)", () -> solve(prob_ode, Vern9(); adaptive = false, dt, dense = false, save_everystep = false)),
    ("Vern9 (adaptive)", () -> solve(prob_ode, Vern9(); save_everystep = false)),
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

# ## Visualization

f = Figure(size = (800, 800), fontsize = 18)

ax1 = Axis(f[1, 1],
    title = "Solver Performance Comparison: Time",
    ylabel = "Median Time (s)",
    xticklabelrotation = π/4,
    xticklabelcolor = :transparent
)
barplot!(ax1, eachindex(results_time), results_time, color = eachindex(results_time), colormap = :viridis)
ax1.xticks = (eachindex(names), names)

ax2 = Axis(f[2, 1],
    title = "Solver Performance Comparison: Memory",
    ylabel = "Median Memory (Bytes)",
    xticklabelrotation = π/4
)
barplot!(ax2, eachindex(results_mem), results_mem, color = eachindex(results_mem), colormap = :viridis)
ax2.xticks = (eachindex(names), names)

f = DisplayAs.PNG(f) #hide

# The Boris method is typically faster and consumes less memory. However, in practice, it is pretty hard to find an optimal algorithm.
# When calling OrdinaryDiffEq.jl, we recommend using `Vern9()` as a starting point instead of `Tsit5()`, especially combined with adaptive timestepping. Further fine-grained control includes setting `dtmax`, `reltol`, and `abstol` in the `solve` method.

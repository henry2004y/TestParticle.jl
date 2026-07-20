# # Multithreaded Parallelization
#
# This demo shows the single-process ensemble parallelization strategies:
# `EnsembleSerial()` (no parallelism) and `EnsembleThreads()` (trajectories split
# across the Julia threads of the current session). Both run entirely inside one
# Julia process, so no extra setup is required. The achievable `EnsembleThreads()`
# speedup is bounded by `Threads.nthreads()`, so start Julia with e.g.
# `julia -t auto` (or set `JULIA_NUM_THREADS`) to use more than one thread.
#
# We benchmark one SciML solver (`Vern6`) and one Boris solver (`Boris`) to show
# that the *same* `ensemblealg` argument drives both.

import DisplayAs #hide
using TestParticle
import TestParticle as TP
using OrdinaryDiffEq
using StaticArrays
using Statistics
using Chairmarks
using CairoMakie
CairoMakie.activate!(type = "png") #hide

nthreads = Threads.nthreads()
println("Running with $nthreads Julia thread(s).")

# ## Simulation Setup
#
# We trace an ensemble of electrons in a uniform magnetic field, sampling each
# particle's initial velocity from a Maxwellian distribution. A fixed random seed
# makes the sampling reproducible and identical across the two parallelization
# modes and both solvers.

const Bmag = 0.01
uniform_B(x) = SA[0.0, 0.0, Bmag]
const zero_E = ZeroField()

const param = prepare(zero_E, uniform_B; species = Electron)
const q2m = TP.get_q2m(param)

const tperiod = 2π / (abs(q2m) * Bmag)
const vth = 1.0e5

const tspan = (0.0, 20 * tperiod)
const dt = tperiod / 12
const stateinit = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

const trajectories = 512
const seed = 1234;

# `prob_func` receives an `EnsembleContext` `ctx`, whose `ctx.rng` is a per-trajectory
# random number generator seeded deterministically from `seed`. The same
# `prob_func` works for both the SciML `ODEProblem` and the Boris `TraceProblem`.
# We draw the initial velocity from a 3D Maxwellian (independent Gaussian
# components, std = vth) for a self-contained, dependency-free setup.

function prob_func(prob, ctx)
    v = SVector{3}(
        vth * randn(ctx.rng),
        vth * randn(ctx.rng),
        vth * randn(ctx.rng),
    )
    return remake(prob; u0 = vcat(SVector{3}(prob.u0[1:3]), v))
end;

# ## SciML Solver
#
# For the SciML path we wrap the base `ODEProblem` in an `EnsembleProblem` and pass
# the desired `ensemblealg` as the third positional argument of `solve`.

const prob_ode = ODEProblem(trace, stateinit, tspan, param)
const ensemble_ode = EnsembleProblem(prob_ode; prob_func, safetycopy = false)

run_ode(ealg) = solve(
    ensemble_ode, Vern6(), ealg;
    trajectories, saveat = dt, seed
);

# ## Boris Solver
#
# The native Boris pusher uses `TraceProblem` and accepts the identical
# `ensemblealg` argument. It additionally requires a fixed timestep `dt`.

const prob_boris = TraceProblem(stateinit, tspan, param; prob_func)

run_boris(ealg) = TP.solve(
    prob_boris, Boris(), ealg;
    dt, trajectories, savestepinterval = 1, seed
)

# ## Benchmarking
#
# We use [Chairmarks.jl](https://github.com/LilithHafner/Chairmarks.jl) to measure
# the median wall-clock time of each configuration.

median_time(b) = median(b).time

benchmarks = [
    ("Vern6", "Serial", () -> run_ode(EnsembleSerial())),
    ("Vern6", "Threads", () -> run_ode(EnsembleThreads())),
    ("Boris", "Serial", () -> run_boris(EnsembleSerial())),
    ("Boris", "Threads", () -> run_boris(EnsembleThreads())),
]

labels = String[]
modes = String[]
times = Float64[]

for (solver, mode, func) in benchmarks
    println("Benchmarking $solver ($mode)...")
    t = median_time(@be func())
    push!(labels, "$solver\n$mode")
    push!(modes, mode)
    push!(times, t)
end

speedup_ode = times[1] / times[2]
speedup_boris = times[3] / times[4]
println("Vern6 speedup: $(round(speedup_ode, digits = 2))x")
println("Boris speedup: $(round(speedup_boris, digits = 2))x")

# ## Results
#
# The bar chart shows the elapsed time for each solver in serial vs. threaded mode.
# With more than one thread, the threaded bars should be noticeably shorter.

colors = Makie.wong_colors()
mode_color = Dict("Serial" => colors[1], "Threads" => colors[2])

f = Figure(size = (1000, 600), fontsize = 24)
ax = Axis(
    f[1, 1],
    title = "Ensemble tracing time ($trajectories trajectories, $nthreads threads)",
    xlabel = "Solver / Mode",
    ylabel = "Elapsed time (s)",
    xticks = (1:length(labels), labels),
    ygridstyle = :dash
)

barplot!(
    ax, 1:length(times), times;
    color = [mode_color[m] for m in modes]
)

for (i, t) in enumerate(times)
    text!(
        ax, i, t; text = string(round(t, digits = 3), " s"),
        align = (:center, :bottom), offset = (0, 4), fontsize = 16
    )
end

elements = [
    PolyElement(polycolor = mode_color["Serial"]),
    PolyElement(polycolor = mode_color["Threads"]),
]
Legend(
    f[1, 1], elements, ["Serial", "Threads"], "Mode";
    framevisible = false, halign = :right, valign = :top,
    tellwidth = false, tellheight = false
)

f = DisplayAs.PNG(f) #hide

# When scaling beyond a single machine (or a single NUMA node) the same interface
# extends to multiple processes; see the [Distributed parallelization](@ref)
# demo for the `EnsembleDistributed` / `EnsembleSplitThreads` workflow.

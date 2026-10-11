"""
Shared configuration, problem setup and CSV I/O for the parallel scaling benchmarks.

`run_scaling.jl` measures every point with the same problem and writes one common
schema, so the `threads`, `distributed` and `multinode` modes produce directly
comparable tables. The field functions and `prob_func` live here rather than in the
driver because `TestParticle` serializes the problem to the workers, which load this
file again through `run_scaling.jl`.

Problem size and repetition are read from the environment so that a sweep can be
reduced without editing code:

- `TP_N_PARTICLES`: trajectories per timed run (default 16384).
- `TP_N_STEPS`: Boris steps per trajectory (default 250000).
- `TP_N_SAMPLES`: timed repetitions, median is reported (default 3).
- `TP_N_WARMUP`: trajectories in the small compile warmup (default 16). Every point then
  runs one full-size untimed warmup, see `measure`.
"""
module ScalingCommon

using TestParticle
using StaticArrays
using Statistics
using Printf
using DelimitedFiles

export ScalingPoint, CSV_HEADER, N_PARTICLES, N_STEPS, N_SAMPLES, throughput,
    make_problem, measure, speedups, write_points, read_points, progress

env_int(name, default) = parse(Int, get(ENV, name, string(default)))

const N_PARTICLES = env_int("TP_N_PARTICLES", 16_384)
# Sized so that the 128-thread point lands near 3 s and the 1-thread baseline near
# 280 s on a 128-core node, about 4e9 particle-steps per run. A sweep costs roughly
# twice its slowest point, so this keeps one near half an hour. Raise it back towards
# 1_000_000 for a final run at the historical problem size.
const N_STEPS = env_int("TP_N_STEPS", 250_000)
const N_SAMPLES = env_int("TP_N_SAMPLES", 3)
const N_WARMUP = env_int("TP_N_WARMUP", 16)

const DT = 1.0e-9
const TSPAN = (0.0, DT * N_STEPS)
# 100 output steps per particle.
const SAVEAT = TSPAN[2] / 100

uniform_B(x) = SA[0.0, 0.0, 1.0e-8]
uniform_E(x) = SA[0.0, 0.0, 0.0]

prob_func(prob, ctx) = remake(
    prob;
    u0 = [prob.u0[1], prob.u0[2], prob.u0[3], (ctx.sim_id / 1000.0) * 1.0e5, 0.0, 0.0]
)

"""
    make_problem() -> TraceProblem

Uniform field tracing problem shared by every benchmark mode. `prob_func` only changes
the parallel velocity, so every trajectory costs the same.
"""
function make_problem()
    param = prepare(uniform_E, uniform_B; species = Proton)
    stateinit = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
    return TraceProblem(stateinit, TSPAN, param; prob_func)
end

# Millions of particle-steps per second.
throughput(t_s) = N_PARTICLES * N_STEPS / (t_s * 1.0e6)

"""
    ScalingPoint

One measured point of a scaling curve. `count` is the swept axis and `units` names it
(`"threads"`, `"workers"` or `"nodes"`); `workers` and `threads_per_worker` record how
the point was provisioned, and `ensemble_alg` is the algorithm that was timed.
`particles` and `steps` describe the workload, so measurements taken at different
problem sizes stay distinguishable.
"""
struct ScalingPoint
    mode::String
    units::String
    count::Int
    workers::Int
    threads_per_worker::Int
    ensemble_alg::String
    particles::Int
    steps::Int
    median_s::Float64
end

const CSV_HEADER = [
    "mode", "units", "count", "workers", "threads_per_worker", "ensemble_alg",
    "particles", "steps", "median_s", "msteps_per_s",
]

# Progress lines should appear promptly even when SLURM redirects stdout to a file.
function progress(msg::AbstractString)
    println(msg)
    flush(stdout)
    return nothing
end

"""
    measure(ensemblealg, prob) -> Float64

Median wall time in seconds over `N_SAMPLES` runs. Two untimed runs come first: a small
one with `N_WARMUP` trajectories that compiles the solver on every worker, then one at the
full `N_PARTICLES` size. Without the second, the first timed sample is also the first time
the process handles the whole ensemble, and it can come out twice as slow as the rest,
which biases the median. Every sample is reported as it finishes, since one sample of the
default problem takes minutes.
"""
function measure(ensemblealg, prob)
    progress(@sprintf("  warmup: %d trajectories (compile)", N_WARMUP))
    solve(prob, Boris(), ensemblealg; trajectories = N_WARMUP, dt = DT, saveat = SAVEAT)
    warmup_s = @elapsed solve(
        prob, Boris(), ensemblealg; trajectories = N_PARTICLES, dt = DT, saveat = SAVEAT
    )
    progress(
        @sprintf("  warmup: %d trajectories, %.3f s, discarded", N_PARTICLES, warmup_s)
    )
    samples = Vector{Float64}(undef, N_SAMPLES)
    for i in eachindex(samples)
        GC.gc()
        samples[i] = @elapsed solve(
            prob, Boris(), ensemblealg;
            trajectories = N_PARTICLES, dt = DT, saveat = SAVEAT
        )
        progress(
            @sprintf(
                "  sample %d/%d: %.3f s (%.2f M particle-steps/s)",
                i, N_SAMPLES, samples[i], throughput(samples[i])
            )
        )
    end
    t = median(samples)
    progress(@sprintf("  median: %.3f s", t))
    return t
end

"""
    speedups(points) -> Vector{Float64}

Speedup of every point over the first point with the same `ensemble_alg`, so each curve
is normalised to its own smallest configuration.
"""
function speedups(points)
    baseline = Dict{String, Float64}()
    for p in points
        get!(baseline, p.ensemble_alg, p.median_s)
    end
    return [baseline[p.ensemble_alg] / p.median_s for p in points]
end

"""
    write_points(path, points)

Write the header and one row per point as plain CSV, which keeps the file readable for
`plot_scaling.jl`.
"""
function write_points(path::AbstractString, points)
    open(path, "w") do io
        println(io, join(CSV_HEADER, ','))
        for p in points
            row = (
                p.mode, p.units, string(p.count), string(p.workers),
                string(p.threads_per_worker), p.ensemble_alg, string(p.particles),
                string(p.steps), @sprintf("%.6f", p.median_s),
                @sprintf("%.3f", throughput(p.median_s)),
            )
            println(io, join(row, ','))
        end
    end
    return path
end

"""
    read_points(path) -> Vector{ScalingPoint}

Parse a CSV written by [`write_points`](@ref). A file that holds only a header yields an
empty vector.
"""
function read_points(path::AbstractString)
    data = readdlm(path, ',')
    size(data, 1) < 1 && return ScalingPoint[]
    header = String.(data[1, :])
    header == CSV_HEADER || error(
        "$path has columns $(join(header, ',')); expected $(join(CSV_HEADER, ','))"
    )
    return [parse_point(data[i, :]) for i in 2:size(data, 1)]
end

parse_point(row) = ScalingPoint(
    String(row[1]), String(row[2]), Int(Float64(row[3])), Int(Float64(row[4])),
    Int(Float64(row[5])), String(row[6]), Int(Float64(row[7])), Int(Float64(row[8])),
    Float64(row[9]),
)

end # module ScalingCommon

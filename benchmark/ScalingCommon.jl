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

Each point records the median of its repetitions together with the fastest and slowest one,
so a figure can show how much the repetitions disagreed.
"""
module ScalingCommon

using TestParticle
using StaticArrays
using Statistics
using Printf
using DelimitedFiles

export ScalingPoint, CSV_HEADER, N_PARTICLES, N_STEPS, N_SAMPLES, throughput,
    make_problem, measure, speedups, write_points, read_points, merge_points, progress

# An exported-but-empty variable counts as unset rather than as a parse error.
env_int(name, default) = (v = strip(get(ENV, name, "")); isempty(v) ? default : parse(Int, v))

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
problem sizes stay distinguishable. `min_s` and `max_s` are the fastest and slowest
repetition behind `median_s`.
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
    min_s::Float64
    max_s::Float64
end

const CSV_HEADER = [
    "mode", "units", "count", "workers", "threads_per_worker", "ensemble_alg",
    "particles", "steps", "median_s", "min_s", "max_s", "msteps_per_s",
]

# Written before the sample range was recorded. Still accepted on read, with `min_s` and
# `max_s` set to the median, so a figure drawn from older measurements carries no bars
# rather than refusing to draw.
const LEGACY_CSV_HEADER = [
    "mode", "units", "count", "workers", "threads_per_worker", "ensemble_alg",
    "particles", "steps", "median_s", "msteps_per_s",
]

"""
    throughput(p::ScalingPoint) -> Float64

Millions of particle-steps per second for a recorded point, taken from the workload that
the point carries. Use this for stored data: `throughput(t_s)` uses the configured
`N_PARTICLES` and `N_STEPS`, so it is only right for a point measured in this session, and
a collected sweep from a run at a different size would report the wrong rate.
"""
throughput(p::ScalingPoint) = p.particles * p.steps / (p.median_s * 1.0e6)

# Progress lines should appear promptly even when SLURM redirects stdout to a file.
function progress(msg::AbstractString)
    println(msg)
    flush(stdout)
    return nothing
end

"""
    measure(ensemblealg, prob) -> (; median, min, max)

Timing over `N_SAMPLES` runs: the median wall time in seconds plus the fastest and slowest
repetition. Two untimed runs come first: a small
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
    progress(
        @sprintf(
            "  median: %.3f s (range %.3f to %.3f)",
            t, minimum(samples), maximum(samples)
        )
    )
    return (; median = t, min = minimum(samples), max = maximum(samples))
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
                @sprintf("%.6f", p.min_s), @sprintf("%.6f", p.max_s),
                @sprintf("%.3f", throughput(p)),
            )
            println(io, join(row, ','))
        end
    end
    return path
end

"""
    read_points(path) -> Vector{ScalingPoint}

Parse a CSV written by [`write_points`](@ref). A file that holds only a header yields an
empty vector. A file written before the sample range was recorded is still accepted, with
`min_s` and `max_s` set to the median, so it draws without error bars.
"""
function read_points(path::AbstractString)
    data = readdlm(path, ',')
    size(data, 1) < 1 && return ScalingPoint[]
    header = String.(data[1, :])
    legacy = header == LEGACY_CSV_HEADER
    header == CSV_HEADER || legacy || error(
        "$path has columns $(join(header, ',')); expected $(join(CSV_HEADER, ','))"
    )
    return [parse_point(data[i, :], legacy) for i in 2:size(data, 1)]
end

function parse_point(row, legacy)
    median_s = Float64(row[9])
    min_s = legacy ? median_s : Float64(row[10])
    max_s = legacy ? median_s : Float64(row[11])
    return ScalingPoint(
        String(row[1]), String(row[2]), Int(Float64(row[3])), Int(Float64(row[4])),
        Int(Float64(row[5])), String(row[6]), Int(Float64(row[7])), Int(Float64(row[8])),
        median_s, min_s, max_s,
    )
end

"""
    merge_points(existing, fresh) -> Vector{ScalingPoint}

Fold freshly measured points into the table already on disk. A stored point whose mode,
count and ensemble algorithm match a fresh point is replaced and the rest are kept, so
re-measuring one suspicious count refreshes that row instead of discarding the sweep.
Points stored at a different workload are dropped, since one table describes one problem
size.
"""
function merge_points(existing, fresh)
    isempty(existing) && return fresh
    workload = (first(fresh).particles, first(fresh).steps)
    key(p) = (p.mode, p.count, p.ensemble_alg)
    fresh_keys = Set(key(p) for p in fresh)
    kept = filter(existing) do p
        (p.particles, p.steps) == workload && key(p) ∉ fresh_keys
    end
    dropped = count(p -> (p.particles, p.steps) != workload, existing)
    dropped > 0 && progress(
        "Dropped $dropped stored point(s) measured at a different workload; " *
            "the table now holds $(workload[1]) particles, $(workload[2]) steps."
    )
    return vcat(kept, fresh)
end

end # module ScalingCommon

# Scaling plots
#
# Reads the CSV files produced by run_scaling.jl and renders the corresponding figures.
# All three benchmark modes share one CSV schema, so every figure is drawn from the same
# reader. Select figures by name, or use `all`:
#
#   julia --project=<env-with-CairoMakie> benchmark/plot_scaling.jl             # cpu
#   julia --project=<env-with-CairoMakie> benchmark/plot_scaling.jl threads
#   julia --project=<env-with-CairoMakie> benchmark/plot_scaling.jl distributed
#   julia --project=<env-with-CairoMakie> benchmark/plot_scaling.jl multinode
#   julia --project=<env-with-CairoMakie> benchmark/plot_scaling.jl all
#
# Options:
#   --out=<path>      output figure (threads and distributed only)
#   --file=<csv>      input CSV, relative to the benchmark directory (single figure only)
#   --machine=<name>  machine name shown in the title
#
# Figures:
#   cpu          every CPU sweep on one axis: aggregate throughput against the cores used,
#                which is the only axis all three sweeps share
#   threads      EnsembleThreads alone
#   distributed  EnsembleDistributed alone
#   multinode    EnsembleSplitThreads against EnsembleDistributed across nodes
#
# From Julia, the same figures are available as
#
#   include("benchmark/plot_scaling.jl")
#   plot_scaling(:cpu)
#
# Measurement scripts do not depend on CairoMakie, so run this with an environment that
# provides it.

using CairoMakie
using DelimitedFiles

const BENCH_DIR = @__DIR__
# The docs figure belongs to the docs in the repository; a standalone copy of these
# scripts keeps it next to the measurements instead.
const FIG_DIR = isdir(joinpath(BENCH_DIR, "..", "docs")) ?
    joinpath(BENCH_DIR, "..", "docs", "src", "figures") : BENCH_DIR

const DEFAULT_MACHINE = "Perlmutter"

# Mirrors ScalingCommon.CSV_HEADER. Checked on read so that a CSV from an older
# run_scaling.jl fails loudly instead of being silently misread.
const CSV_COLUMNS = [
    "mode", "units", "count", "workers", "threads_per_worker", "ensemble_alg",
    "particles", "steps", "median_s", "min_s", "max_s", "msteps_per_s",
]

# Written before the sample range was recorded; still accepted, with `min_s` and `max_s`
# equal to the median, so older measurements draw without error bars.
const LEGACY_CSV_COLUMNS = [
    "mode", "units", "count", "workers", "threads_per_worker", "ensemble_alg",
    "particles", "steps", "median_s", "msteps_per_s",
]

const INPUT_FILES = Dict(
    :threads => "threads_scaling.csv",
    :distributed => "distributed_scaling.csv",
    :multinode => "multinode_scaling.csv",
)

const OUTPUT_FILES = Dict(
    :threads => joinpath(BENCH_DIR, "threads_scaling.png"),
    :distributed => joinpath(BENCH_DIR, "distributed_scaling.png"),
    :multinode => joinpath(BENCH_DIR, "multinode_scaling.png"),
    :cpu => joinpath(FIG_DIR, "scaling_cpu.png"),
)

"""
    read_points(file) -> NamedTuple

Read a scaling CSV written by run_scaling.jl, with `file` relative to the benchmark
directory. Rows are returned in file order, which run_scaling.jl already sorts by the
swept count. An input that holds only a header yields empty vectors. A file written before
the sample range was recorded is accepted with `min_s` and `max_s` equal to `median_s`, so
it draws without error bars.
"""
function read_points(file::AbstractString)
    data = readdlm(joinpath(BENCH_DIR, file), ',')
    size(data, 1) < 2 && return empty_points()
    names = String.(data[1, :])
    legacy = names == LEGACY_CSV_COLUMNS
    names == CSV_COLUMNS || legacy || error(
        "$file was written by an older run_scaling.jl; columns $(names), " *
            "expected $(CSV_COLUMNS). Regenerate the measurements."
    )
    rows = data[2:end, :]
    median_s = Float64.(rows[:, 9])
    min_s = legacy ? median_s : Float64.(rows[:, 10])
    max_s = legacy ? median_s : Float64.(rows[:, 11])
    return (
        mode = String.(rows[:, 1]),
        units = String.(rows[:, 2]),
        count = Int.(Float64.(rows[:, 3])),
        workers = Int.(Float64.(rows[:, 4])),
        threads_per_worker = Int.(Float64.(rows[:, 5])),
        ensemble_alg = String.(rows[:, 6]),
        particles = Int(Float64(rows[1, 7])),
        steps = Int(Float64(rows[1, 8])),
        median_s = median_s,
        min_s = min_s,
        max_s = max_s,
    )
end

empty_points() = (
    mode = String[],
    units = String[],
    count = Int[],
    workers = Int[],
    threads_per_worker = Int[],
    ensemble_alg = String[],
    particles = 0,
    steps = 0,
    median_s = Float64[],
    min_s = Float64[],
    max_s = Float64[],
)

function ideal_line!(ax, counts)
    x = [minimum(counts), maximum(counts)]
    lines!(
        ax, x, Float64.(x);
        color = :black, linestyle = :dash, linewidth = 2, label = "Ideal Scaling"
    )
    return nothing
end

"""
    with_errorbars!(ax, x, y, y_low, y_high)

Draw asymmetric error bars spanning `y_low` to `y_high` at `x`. Nothing is drawn when both
bounds coincide with `y`, which is what a CSV measured before the sample range was recorded
produces, so those figures look exactly as they did.
"""
function with_errorbars!(ax, x, y, y_low, y_high)
    all(y_low .== y) && all(y_high .== y) && return nothing
    errorbars!(
        ax, x, y, y .- y_low, y_high .- y;
        color = :black, whiskerwidth = 10, linewidth = 1.5
    )
    return nothing
end

function mode_xlabel(mode::Symbol, data)
    mode == :threads && return "Number of Threads"
    mode == :distributed && return "Number of Workers"
    return "Number of Nodes ($(data.threads_per_worker[1]) threads per node)"
end

function mode_title(mode::Symbol)
    mode == :threads && return "Thread Strong Scaling"
    mode == :distributed && return "Distributed Strong Scaling"
    return "Multi-node Strong Scaling"
end

"""
    plot_mode(mode::Symbol; kwargs...)

Plot the speedup of every ensemble algorithm in one mode's CSV against the swept count.
`mode` selects the default input and output files. The thread and worker sweeps use a
log2 x-axis, the node sweep stays linear.

# Keywords
- `file`: CSV table to read, relative to the benchmark directory.
- `out_path`: file to save the figure to.
- `machine`: machine name shown in the title.
"""
function plot_mode(
        mode::Symbol;
        file = INPUT_FILES[mode], out_path = OUTPUT_FILES[mode], machine = DEFAULT_MACHINE
    )
    data = read_points(file)
    isempty(data.count) && error("$file holds no measurements.")
    title = mode_title(mode)
    log_x = mode != :multinode

    fig = Figure(; size = (900, 520), fontsize = 20)
    ax = Axis(
        fig[1, 1];
        xscale = log_x ? log2 : identity,
        yscale = log2,
        xlabel = mode_xlabel(mode, data),
        ylabel = "Speedup",
        title = "$(title) ($(data.particles) Particles, $(data.steps) Steps on $(machine))",
        xticks = log_x ? data.count : sort(unique(data.count)),
        xminorticksvisible = true,
        yminorticksvisible = true,
    )

    for alg in unique(data.ensemble_alg)
        sel = data.ensemble_alg .== alg
        base = data.median_s[sel][1]
        y = base ./ data.median_s[sel]
        scatterlines!(ax, data.count[sel], y; label = alg, linewidth = 3, markersize = 14)
        with_errorbars!(
            ax, data.count[sel], y, base ./ data.max_s[sel], base ./ data.min_s[sel]
        )
    end
    ideal_line!(ax, data.count)

    axislegend(ax; position = :lt)

    mkpath(dirname(out_path))
    save(out_path, fig)
    println("$(title) plot saved to: ", out_path)
    return out_path
end

"""
    plot_cpu(; kwargs...)

Plot every CPU sweep on one axis: aggregate throughput against the number of cores used.
Throughput is a rate, so the single node sweeps and the multi-node sweep share an axis even
though they were measured separately. Each multi-node point contributes `threads_per_worker`
cores per node, so that series extends the axis past a single node.

A dashed line shows ideal scaling at the measured single-core rate for comparison.

# Keywords
- `out_path`: file to save the figure to.
- `machine`: machine name shown in the title.
"""
function plot_cpu(; out_path = OUTPUT_FILES[:cpu], machine = DEFAULT_MACHINE)
    for (key, f) in INPUT_FILES
        isfile(joinpath(BENCH_DIR, f)) ||
            error("$(f) is missing; run `run_scaling.jl $(key)` first.")
    end
    threads = read_points(INPUT_FILES[:threads])
    dist = read_points(INPUT_FILES[:distributed])
    multi = read_points(INPUT_FILES[:multinode])
    isempty(threads.count) && error("threads_scaling.csv holds no measurements.")
    isempty(dist.count) && error("distributed_scaling.csv holds no measurements.")
    isempty(multi.count) && error("multinode_scaling.csv holds no measurements.")

    multi_cores = multi.count .* multi.threads_per_worker[1]
    series = (
        ("Multithreading, $(threads.particles) particles", threads.count, threads),
        ("Distributed, $(dist.particles) particles", dist.count, dist),
        ("Multi-node SplitThreads, $(multi.particles) particles", multi_cores, multi),
    )
    xticks = sort(unique(vcat([s[2] for s in series]...)))

    fig = Figure(; size = (1000, 620), fontsize = 19)
    ax = Axis(
        fig[1, 1];
        xscale = log2,
        yscale = log2,
        xlabel = "CPU Cores Used",
        ylabel = "Particle-steps per Second (millions)",
        title = "CPU Scaling Throughput on $(machine)",
        xticks = xticks,
        xminorticksvisible = true,
        yminorticksvisible = true,
    )

    for (label, x, data) in series
        total = data.particles * data.steps
        y = total ./ (data.median_s .* 1.0e6)
        scatterlines!(ax, x, y; label = label, linewidth = 3, markersize = 12)
        with_errorbars!(
            ax, x, y, total ./ (data.max_s .* 1.0e6), total ./ (data.min_s .* 1.0e6)
        )
    end

    base = threads.particles * threads.steps / (threads.median_s[1] * 1.0e6)
    ideal_x = [minimum(xticks), maximum(xticks)]
    lines!(
        ax, ideal_x, base .* Float64.(ideal_x);
        color = :black, linestyle = :dash, linewidth = 2, label = "Ideal (from 1 core)"
    )

    axislegend(ax; position = :lt)

    mkpath(dirname(out_path))
    save(out_path, fig)
    println("CPU throughput plot saved to: ", out_path)
    return out_path
end

plot_threads(; kwargs...) = plot_mode(:threads; kwargs...)

plot_distributed(; kwargs...) = plot_mode(:distributed; kwargs...)

plot_multinode(; kwargs...) = plot_mode(:multinode; kwargs...)

const PLOT_METHODS = Dict(
    :threads => plot_threads,
    :distributed => plot_distributed,
    :multinode => plot_multinode,
    :cpu => plot_cpu,
)

const METHOD_NAMES = ["cpu", "threads", "distributed", "multinode"]

"""
    plot_scaling(method::Symbol; kwargs...)

Draw one of the scaling figures, forwarding keyword arguments to its plotting function.
Use `:cpu` for all three sweeps on one axis, or `:threads`, `:distributed` and
`:multinode` for a single mode.
"""
function plot_scaling(method::Symbol; kwargs...)
    plot = get(PLOT_METHODS, method) do
        error("Unknown plot method $(method); choose one of $(join(METHOD_NAMES, ", ")).")
    end
    return plot(; kwargs...)
end

function parse_args(args)
    methods = Symbol[]
    options = Dict{Symbol, String}()
    for arg in args
        if arg in METHOD_NAMES
            push!(methods, Symbol(arg))
        elseif arg == "all"
            push!(methods, :all)
        elseif startswith(arg, "--out=")
            options[:out_path] = split(arg, "="; limit = 2)[2]
        elseif startswith(arg, "--file=")
            options[:file] = split(arg, "="; limit = 2)[2]
        elseif startswith(arg, "--machine=")
            options[:machine] = split(arg, "="; limit = 2)[2]
        else
            error("Unrecognized argument: $(arg). See the header of this file.")
        end
    end

    if isempty(methods)
        push!(methods, :cpu)
    elseif :all in methods
        haskey(options, :out_path) &&
            error("--out is ambiguous for `all`; plot one method at a time.")
        haskey(options, :file) &&
            error("--file is ambiguous for `all`; plot one method at a time.")
        methods = Symbol.(METHOD_NAMES)
    end
    return methods, options
end

function main(args = ARGS)
    methods, options = parse_args(args)
    for method in methods
        plot_scaling(method; options...)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

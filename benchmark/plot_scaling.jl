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
#   julia --project=<env-with-CairoMakie> benchmark/plot_scaling.jl gpu
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
#   gpu          GPU pusher throughput on uniform and 3D grid fields
#
# From Julia, the same figures are available as
#
#   include("benchmark/plot_scaling.jl")
#   plot_scaling(:cpu)
#   plot_scaling(:gpu)
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
    :gpu => "perlmutter-benchmark-result-gpu.md",
)

const OUTPUT_FILES = Dict(
    :threads => joinpath(BENCH_DIR, "threads_scaling.png"),
    :distributed => joinpath(BENCH_DIR, "distributed_scaling.png"),
    :multinode => joinpath(BENCH_DIR, "multinode_scaling.png"),
    :cpu => joinpath(FIG_DIR, "scaling_cpu.png"),
    :gpu => joinpath(FIG_DIR, "scaling_gpu.png"),
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

    fig = Figure(; size = (1000, 450), fontsize = 19)
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

function parse_gpu_number(s::AbstractString)
    m = match(r"([0-9]+(?:\.[0-9]+)?)", s)
    return m === nothing ? NaN : parse(Float64, m.match)
end

function read_gpu_benchmark(file::AbstractString)
    filepath = isabspath(file) ? file : joinpath(BENCH_DIR, file)
    lines = readlines(filepath)

    gpu_device = "NVIDIA A100-SXM4-40GB"
    cpu_threads = "128"
    for line in lines[1:min(10, length(lines))]
        if occursin("GPU Device:", line)
            gpu_device = strip(split(line, ":"; limit = 2)[2])
        elseif occursin("CPU Threads:", line)
            cpu_threads = strip(split(line, ":"; limit = 2)[2])
        end
    end

    analytical_n = Float64[]
    cpu_ser = Float64[]
    cpu_thr = Float64[]
    ka_cpu = Float64[]
    gpu_fp64 = Float64[]
    gpu_fp32 = Float64[]
    ka_cpu_raw = Float64[]
    gpu_fp64_raw = Float64[]
    gpu_fp32_raw = Float64[]

    grid_n = Float64[]
    grid_cpu_thr = Float64[]
    grid_ka_cpu = Float64[]
    grid_gpu = Float64[]
    grid_ka_cpu_raw = Float64[]
    grid_gpu_raw = Float64[]

    large_n = Float64[]
    large_gpu_fp64 = Float64[]
    large_gpu_fp32 = Float64[]
    large_gpu_fp32_raw = Float64[]

    mode = :none
    for line in lines
        s = strip(line)
        if occursin("[1] Analytical", s)
            mode = :s1_base
        elseif occursin("[2] 3D Grid", s)
            mode = :s2_base
        elseif occursin("[3] GPU-only large-scale", s)
            mode = :s3
        elseif occursin("[4] GPU kernel isolation", s)
            mode = :s4
        elseif occursin("Raw output", s)
            if mode in (:s1_base, :s1_other)
                mode = :s1_raw
            elseif mode in (:s2_base, :s2_other)
                mode = :s2_raw
            end
        elseif occursin("Speedup relative to", s) || occursin("Effective throughput", s)
            if mode in (:s1_base, :s1_raw)
                mode = :s1_other
            elseif mode in (:s2_base, :s2_raw)
                mode = :s2_other
            elseif mode == :s3
                mode = :s3_other
            end
        elseif occursin("|", s) && !occursin("Particles", s) && !occursin("Steps", s) &&
                !occursin("---", s)
            parts = strip.(split(s, "|"))
            if mode == :s1_base && length(parts) >= 6
                push!(analytical_n, parse_gpu_number(parts[1]))
                push!(cpu_ser, parse_gpu_number(parts[2]))
                push!(cpu_thr, parse_gpu_number(parts[3]))
                push!(ka_cpu, parse_gpu_number(parts[4]))
                push!(gpu_fp64, parse_gpu_number(parts[5]))
                push!(gpu_fp32, parse_gpu_number(parts[6]))
            elseif mode == :s1_raw && length(parts) >= 4
                push!(ka_cpu_raw, parse_gpu_number(parts[2]))
                push!(gpu_fp64_raw, parse_gpu_number(parts[3]))
                push!(gpu_fp32_raw, parse_gpu_number(parts[4]))
            elseif mode == :s2_base && length(parts) >= 4
                push!(grid_n, parse_gpu_number(parts[1]))
                push!(grid_cpu_thr, parse_gpu_number(parts[2]))
                push!(grid_ka_cpu, parse_gpu_number(parts[3]))
                push!(grid_gpu, parse_gpu_number(parts[4]))
            elseif mode == :s2_raw && length(parts) >= 3
                push!(grid_ka_cpu_raw, parse_gpu_number(parts[2]))
                push!(grid_gpu_raw, parse_gpu_number(parts[3]))
            elseif mode == :s3 && length(parts) >= 4
                push!(large_n, parse_gpu_number(parts[1]))
                push!(large_gpu_fp64, parse_gpu_number(parts[2]))
                push!(large_gpu_fp32, parse_gpu_number(parts[3]))
                push!(large_gpu_fp32_raw, parse_gpu_number(parts[4]))
            end
        end
    end

    return (
        gpu_device = gpu_device,
        cpu_threads = cpu_threads,
        analytical = (
            n = analytical_n, steps = 1000,
            cpu_ser = cpu_ser, cpu_thr = cpu_thr,
            ka_cpu = ka_cpu, ka_cpu_raw = ka_cpu_raw,
            gpu_fp64 = gpu_fp64, gpu_fp64_raw = gpu_fp64_raw,
            gpu_fp32 = gpu_fp32, gpu_fp32_raw = gpu_fp32_raw,
        ),
        grid = (
            n = grid_n, steps = 200,
            cpu_thr = grid_cpu_thr,
            ka_cpu = grid_ka_cpu, ka_cpu_raw = grid_ka_cpu_raw,
            gpu = grid_gpu, gpu_raw = grid_gpu_raw,
        ),
        large = (
            n = large_n, steps = 1000,
            gpu_fp64 = large_gpu_fp64,
            gpu_fp32 = large_gpu_fp32,
            gpu_fp32_raw = large_gpu_fp32_raw,
        ),
    )
end

"""
    plot_gpu(; kwargs...)

Plot GPU pusher throughput on uniform and 3D grid fields against particle count.

# Keywords
- `file`: markdown file holding benchmark output, relative to the benchmark directory.
- `out_path`: file to save the figure to.
- `machine`: machine name shown in the title.
"""
function plot_gpu(;
        file = INPUT_FILES[:gpu], out_path = OUTPUT_FILES[:gpu], machine = DEFAULT_MACHINE
    )
    data = read_gpu_benchmark(file)
    fig = Figure(; size = (1000, 450), fontsize = 18)

    a = data.analytical
    lg = data.large
    tp(n, t_ms, steps) = (n .* steps) ./ (t_ms .* 1000.0)

    ax1 = Axis(
        fig[1, 1];
        xscale = log10,
        yscale = log10,
        xlabel = "Number of Particles",
        ylabel = "Throughput (M particle-steps/s)",
        title = "Analytical Field (1,000 Steps)",
        xticks = ([1e3, 1e4, 1e5, 1e6, 5e6], ["1k", "10k", "100k", "1M", "5M"]),
        xminorticksvisible = true,
        yminorticksvisible = true,
    )

    scatterlines!(
        ax1, a.n, tp(a.n, a.cpu_ser, a.steps);
        label = "CPU Serial (FP64)", color = :gray50, linestyle = :dot,
        linewidth = 2, markersize = 8
    )
    scatterlines!(
        ax1, a.n, tp(a.n, a.cpu_thr, a.steps);
        label = "CPU Threads ($(data.cpu_threads) thr)", color = :darkorange,
        linewidth = 2.5, markersize = 10
    )
    scatterlines!(
        ax1, a.n, tp(a.n, a.ka_cpu_raw, a.steps);
        label = "KA CPU ($(data.cpu_threads) thr, Raw)", color = :chocolate3,
        linestyle = :dash, linewidth = 2.5, markersize = 10
    )

    if !isempty(lg.n)
        fp64_n = vcat(a.n[1:2], lg.n)
        fp64_t = vcat(a.gpu_fp64[1:2], lg.gpu_fp64)
        scatterlines!(
            ax1, fp64_n, tp(fp64_n, fp64_t, 1000);
            label = "GPU FP64", color = :steelblue, linewidth = 2.5, markersize = 10
        )

        fp32_n = vcat(a.n[1:2], lg.n)
        fp32_t = vcat(a.gpu_fp32[1:2], lg.gpu_fp32)
        scatterlines!(
            ax1, fp32_n, tp(fp32_n, fp32_t, 1000);
            label = "GPU FP32", color = :teal, linewidth = 2.5, markersize = 10
        )

        fp32_raw_n = vcat(a.n[1:2], lg.n)
        fp32_raw_t = vcat(a.gpu_fp32_raw[1:2], lg.gpu_fp32_raw)
        scatterlines!(
            ax1, fp32_raw_n, tp(fp32_raw_n, fp32_raw_t, 1000);
            label = "GPU FP32 (Raw)", color = :forestgreen, linewidth = 3, markersize = 11
        )
    else
        scatterlines!(
            ax1, a.n, tp(a.n, a.gpu_fp64, a.steps);
            label = "GPU FP64", color = :steelblue, linewidth = 2.5, markersize = 10
        )
        scatterlines!(
            ax1, a.n, tp(a.n, a.gpu_fp32, a.steps);
            label = "GPU FP32", color = :teal, linewidth = 2.5, markersize = 10
        )
        scatterlines!(
            ax1, a.n, tp(a.n, a.gpu_fp32_raw, a.steps);
            label = "GPU FP32 (Raw)", color = :forestgreen, linewidth = 3, markersize = 11
        )
    end

    axislegend(ax1; position = :rb, labelsize = 12, backgroundcolor = (:white, 0.85))

    ax2 = Axis(
        fig[1, 2];
        xscale = log10,
        yscale = log10,
        xlabel = "Number of Particles",
        ylabel = "Throughput (M particle-steps/s)",
        title = "Numerical Field (32³ Grid, 200 Steps)",
        xticks = ([2e3, 1e4, 5e4], ["2k", "10k", "50k"]),
        yticks = [1e2, 1e3, 1e4],
        xminorticksvisible = true,
        yminorticksvisible = true,
    )

    g = data.grid
    scatterlines!(
        ax2, g.n, tp(g.n, g.cpu_thr, g.steps);
        label = "CPU Threads ($(data.cpu_threads) thr)", color = :darkorange,
        linewidth = 2.5, markersize = 10
    )
    scatterlines!(
        ax2, g.n, tp(g.n, g.ka_cpu, g.steps);
        label = "KA CPU ($(data.cpu_threads) thr)", color = :chocolate3,
        linewidth = 2.5, markersize = 10
    )
    scatterlines!(
        ax2, g.n, tp(g.n, g.ka_cpu_raw, g.steps);
        label = "KA CPU ($(data.cpu_threads) thr, Raw)", color = :chocolate3,
        linestyle = :dash, linewidth = 2.5, markersize = 10
    )
    scatterlines!(
        ax2, g.n, tp(g.n, g.gpu, g.steps);
        label = "GPU (FP32)", color = :teal, linewidth = 2.5, markersize = 10
    )
    scatterlines!(
        ax2, g.n, tp(g.n, g.gpu_raw, g.steps);
        label = "GPU (FP32 Raw)", color = :forestgreen, linewidth = 3, markersize = 11
    )

    ylims!(ax2, 70, 25000)
    axislegend(ax2; position = :lt, labelsize = 12, backgroundcolor = (:white, 0.85))

    Label(
        fig[0, :], "GPU Throughput on $(machine) ($(data.gpu_device))";
        fontsize = 20, font = :bold
    )

    mkpath(dirname(out_path))
    save(out_path, fig)
    println("GPU throughput plot saved to: ", out_path)
    return out_path
end

const PLOT_METHODS = Dict(
    :threads => plot_threads,
    :distributed => plot_distributed,
    :multinode => plot_multinode,
    :cpu => plot_cpu,
    :gpu => plot_gpu,
)

const METHOD_NAMES = ["cpu", "threads", "distributed", "multinode", "gpu"]

"""
    plot_scaling(method::Symbol; kwargs...)

Draw one of the scaling figures, forwarding keyword arguments to its plotting function.
Use `:cpu` for CPU core scaling, `:gpu` for GPU pusher scaling, or `:threads`,
`:distributed` and `:multinode` for a single ensemble mode.
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

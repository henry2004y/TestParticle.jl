# Scaling plots
#
# Reads the CSV files produced by the run_scaling_*.jl drivers and renders the
# corresponding strong scaling speedup figures. Every figure has its own plotting
# function; select which one to draw with a keyword (or the matching CLI flag):
#
#   julia --project=benchmark benchmark/plot_scaling.jl                # combined
#   julia --project=benchmark benchmark/plot_scaling.jl multinode
#   julia --project=benchmark benchmark/plot_scaling.jl all
#   julia --project=benchmark benchmark/plot_scaling.jl --out=/tmp/scaling.png
#
# From Julia, the same figures are available as
#
#   include("benchmark/plot_scaling.jl")
#   plot_scaling(:multinode)
#   plot_combined(; out_path = "figures/combined.png", particles = 65_536)
#
# CairoMakie is not part of the benchmark project, so run this with an environment
# that provides it.

using CairoMakie
using DelimitedFiles

const BENCH_DIR = @__DIR__

const DEFAULT_PARTICLES = 16_384
const DEFAULT_MACHINE = "Perlmutter"

"""
    read_scaling(file)

Read a two-column `counts,times` CSV and return the processor counts together with
the speedup relative to the first entry.
"""
function read_scaling(file::AbstractString)
    data = readdlm(joinpath(BENCH_DIR, file), ',')
    counts = Int.(Float64.(data[:, 1]))
    times = Float64.(data[:, 2])
    return counts, times[1] ./ times
end

"""
    read_multinode(file)

Read a `multinode_scaling.csv` table (with a header row) and return the node counts,
ensemble algorithms, per-algorithm speedups, and the run configuration.
"""
function read_multinode(file::AbstractString)
    rows = readdlm(joinpath(BENCH_DIR, file), ',')[2:end, :]
    return (
        nodes = Int.(Float64.(rows[:, 1])),
        algs = String.(rows[:, 4]),
        speedups = Float64.(rows[:, 8]),
        threads_per_worker = Int(Float64(rows[1, 3])),
        particles = Int(Float64(rows[1, 5])),
    )
end

"""
    plot_combined(; kwargs...)

Plot the strong scaling speedup of the multithreading and distributed backends of a
single node run against each other.

# Keywords
- `out_path`: file to save the figure to.
- `particles`: particle count shown in the title.
- `machine`: machine name shown in the title.
"""
function plot_combined(;
    out_path = joinpath(BENCH_DIR, "..", "docs", "src", "figures", "scaling_combined.png"),
    particles = DEFAULT_PARTICLES,
    machine = DEFAULT_MACHINE,
)
    t_counts, threads_speedup = read_scaling("threads_scaling.csv")
    d_counts, dist_speedup = read_scaling("distributed_scaling.csv")
    if t_counts != d_counts
        error("threads_scaling.csv and distributed_scaling.csv cover different counts.")
    end
    counts = t_counts

    fig = Figure(; size = (1000, 500), fontsize = 20)
    ax = Axis(
        fig[1, 1];
        xscale = log2,
        yscale = log2,
        xlabel = "Number of Threads / Workers",
        ylabel = "Speedup",
        title = "Strong Scaling ($(particles) Particles on $(machine))",
        xticks = counts,
        yticks = counts,
        xminorticksvisible = true,
        yminorticksvisible = true,
    )

    scatterlines!(
        ax, counts, threads_speedup; label = "Multithreading", linewidth = 3
    )
    scatterlines!(
        ax, counts, dist_speedup; label = "Distributed", linewidth = 3
    )
    lines!(
        ax, counts, Float64.(counts);
        color = :black, linestyle = :dash, label = "Ideal Scaling", linewidth = 2
    )

    axislegend(ax; position = :lt)

    mkpath(dirname(out_path))
    save(out_path, fig)
    println("Combined scaling plot saved to: ", out_path)
    return out_path
end

"""
    plot_multinode(; kwargs...)

Plot the speedup of each ensemble algorithm against the number of nodes.

# Keywords
- `out_path`: file to save the figure to.
- `file`: CSV table to read, relative to the benchmark directory.
"""
function plot_multinode(;
    out_path = joinpath(BENCH_DIR, "multinode_scaling.png"),
    file = "multinode_scaling.csv",
)
    data = read_multinode(file)

    fig = Figure(; size = (900, 520), fontsize = 20)
    ax = Axis(
        fig[1, 1];
        xlabel = "Number of Nodes ($(data.threads_per_worker) threads per node)",
        ylabel = "Speedup",
        title = "Multi-node Strong Scaling ($(data.particles) Particles)",
        xticks = sort(unique(data.nodes)),
        xminorticksvisible = true,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(5),
    )

    for alg in unique(data.algs)
        sel = data.algs .== alg
        scatterlines!(
            ax, data.nodes[sel], data.speedups[sel];
            label = alg, linewidth = 3, markersize = 14
        )
    end

    nmax = maximum(data.nodes)
    lines!(
        ax, [1, nmax], [1, nmax];
        color = :black, linestyle = :dash, linewidth = 2, label = "Ideal Scaling"
    )

    axislegend(ax; position = :lt)

    mkpath(dirname(out_path))
    save(out_path, fig)
    println("Multi-node scaling plot saved to: ", out_path)
    return out_path
end

const PLOT_METHODS = Dict(:combined => plot_combined, :multinode => plot_multinode)

"""
    plot_scaling(method::Symbol; kwargs...)

Draw one of the scaling figures, forwarding keyword arguments to its plotting
function. Use `:combined` for the single node comparison and `:multinode` for the
node sweep.
"""
function plot_scaling(method::Symbol; kwargs...)
    plot = get(PLOT_METHODS, method) do
        error(
            "Unknown plot method $(method); choose one of " *
            join(sort!(collect(keys(PLOT_METHODS)); by = string), ", "),
        )
    end
    return plot(; kwargs...)
end

function parse_args(args)
    methods = Symbol[]
    options = Dict{Symbol, String}()
    for arg in args
        if arg in ("combined", "multinode", "all")
            push!(methods, Symbol(arg))
        elseif startswith(arg, "--out=")
            options[:out_path] = split(arg, "="; limit = 2)[2]
        else
            error(
                "Unrecognized argument: $(arg). " *
                "Choose from combined, multinode, all, --out=PATH."
            )
        end
    end

    isempty(methods) && push!(methods, :combined)
    if :all in methods
        haskey(options, :out_path) &&
            error("--out is ambiguous for `all`; plot one method at a time.")
        methods = [:combined, :multinode]
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
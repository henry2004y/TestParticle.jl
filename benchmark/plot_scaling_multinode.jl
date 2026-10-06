# Multi-node scaling plot
#
# Reads multinode_scaling.csv produced by run_scaling_multinode.jl and plots the speedup
# of each ensemble algorithm against the number of nodes.
#
# CairoMakie is not part of the benchmark project, so run this with an environment that
# provides it, e.g.
#   julia --project benchmark/plot_scaling_multinode.jl

using CairoMakie
using DelimitedFiles

const BENCH_DIR = @__DIR__

data = readdlm(joinpath(BENCH_DIR, "multinode_scaling.csv"), ',')
rows = data[2:end, :]

nodes_all = Int.(Float64.(rows[:, 1]))
algs = unique(String.(rows[:, 4]))
threads_per_worker = Int(Float64(rows[1, 3]))

fig = Figure(; size = (900, 520), fontsize = 20)
ax = Axis(
    fig[1, 1];
    xlabel = "Number of Nodes ($threads_per_worker threads per node)",
    ylabel = "Speedup",
    title = "Multi-node Strong Scaling ($(Int(Float64(rows[1, 5]))) Particles)",
    xticks = sort(unique(nodes_all)),
    xminorticksvisible = true,
    yminorticksvisible = true,
    yminorticks = IntervalsBetween(5),
)

for alg in algs
    sel = String.(rows[:, 4]) .== alg
    scatterlines!(
        ax, nodes_all[sel], Float64.(rows[sel, 8]);
        label = alg, linewidth = 3, markersize = 14
    )
end

nmax = maximum(nodes_all)
lines!(
    ax, [1, nmax], [1, nmax];
    color = :black, linestyle = :dash, linewidth = 2, label = "Ideal Scaling"
)

axislegend(ax; position = :lt)

plot_path = joinpath(BENCH_DIR, "multinode_scaling.png")
save(plot_path, fig)
println("Saved multi-node scaling plot to: ", plot_path)

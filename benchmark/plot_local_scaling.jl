# Local strong-scaling plot (Boris, 16,384 particles).
# Distributed data: official benchmark/distributed_scaling.csv (8 workers, 1 thread each).
# Thread data: clean EnsembleThreads measurement (2,048 particles) x8 to match the
# 16,384-particle workload of the distributed run.
using CairoMakie
using DelimitedFiles

dist_data = readdlm(joinpath(@__DIR__, "distributed_scaling.csv"), ',')
counts = Int.(dist_data[:, 1])
dist_times = dist_data[:, 2]            # seconds, 16384 particles

# Clean EnsembleThreads measurement on 2048 particles; scale to 16384.
thread_raw = Dict(1 => 11.718, 2 => 7.111, 4 => 4.508, 8 => 2.799, 9 => 3.021)
thread_counts = sort(collect(keys(thread_raw)))
thread_times = [thread_raw[c] * 8.0 for c in thread_counts]

dist_speedup = dist_times[1] ./ dist_times
thread_speedup = thread_times[1] ./ thread_times

fig = Figure(size = (1000, 560), fontsize = 20)
ax = Axis(
    fig[1, 1],
    xscale = log2,
    yscale = log2,
    xlabel = "Number of Threads / Workers",
    ylabel = "Speedup",
    title = "Boris Strong Scaling (16,384 particles, local 8-core machine)",
    xticks = [1, 2, 4, 8, 9],
    yticks = [1, 2, 4, 8, 9],
    xminorticksvisible = true,
    yminorticksvisible = true,
)

scatterlines!(ax, thread_counts, thread_speedup; label = "EnsembleThreads", linewidth = 3, markersize = 10)
scatterlines!(ax, counts, dist_speedup; label = "EnsembleDistributed", linewidth = 3, markersize = 10)
lines!(ax, [1, 9], [1.0, 9.0]; color = :black, linestyle = :dash, label = "Ideal", linewidth = 2)

axislegend(ax, position = :lt)
save(joinpath(@__DIR__, "parallel_scaling_local.png"), fig)
println("Saved parallel_scaling_local.png")

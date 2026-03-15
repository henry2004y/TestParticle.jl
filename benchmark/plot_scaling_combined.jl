# Combined Strong Scaling Plot
using CairoMakie
using DelimitedFiles

# Load data from CSV files
threads_data = readdlm(joinpath(@__DIR__, "threads_scaling.csv"), ',')
dist_data = readdlm(joinpath(@__DIR__, "distributed_scaling.csv"), ',')

counts = Int.(threads_data[:, 1])

# Multithreading median times (converted from ms to s for consistency)
threads_times = threads_data[:, 2] ./ 1000.0

# Distributed median times in seconds
dist_times = dist_data[:, 2]

# Calculate speedups
threads_speedup = threads_times[1] ./ threads_times
dist_speedup = dist_times[1] ./ dist_times

# --- Plotting ---
fig = Figure(size = (1000, 500), fontsize = 20)
ax = Axis(
    fig[1, 1],
    xscale = log2,
    yscale = log2,
    xlabel = "Number of Threads / Workers",
    ylabel = "Speedup",
    title = "Strong Scaling (16,384 Particles on Perlmutter)",
    xticks = counts,
    yticks = counts,
    xminorticksvisible = true,
    yminorticksvisible = true,
)

# Multithreading
scatterlines!(ax, counts, threads_speedup, label = "Multithreading", linewidth = 3)

# Distributed
scatterlines!(ax, counts, dist_speedup, label = "Distributed", linewidth = 3)

# Ideal reference
lines!(
    ax, counts, Float64.(counts),
    color = :black, linestyle = :dash, label = "Ideal Scaling", linewidth = 2
)

axislegend(ax, position = :lt)

plot_path = joinpath(@__DIR__, "..", "docs", "src", "figures", "scaling_combined.png")
save(plot_path, fig)
println("Combined scaling plot saved to: ", plot_path)

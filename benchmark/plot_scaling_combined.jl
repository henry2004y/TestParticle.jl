# Combined Strong Scaling Plot
using CairoMakie

# Thread/Worker counts
counts = [2^i for i in 0:8] # 1, 2, 4, 8, 16, 32, 64, 128, 256

# Multithreading median times (converted from ms to s for consistency)
threads_times = [
    151531.9724, 75669.774758, 37841.402016, 19013.7110955,
    9746.676036, 5334.494392, 2685.373244, 1437.751099, 1084.204129,
] ./ 1000.0

# Distributed median times in seconds
dist_times = [
    152.689, 76.496, 38.926, 20.007, 10.517, 5.979, 3.401, 2.647, 2.234,
]

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

plot_path = joinpath(@__DIR__, "scaling_combined.png")
save(plot_path, fig)
println("Combined scaling plot saved to: ", plot_path)

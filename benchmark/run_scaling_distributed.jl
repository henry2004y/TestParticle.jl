# Distributed Strong Scaling Driver
#
# Measures strong scaling of EnsembleDistributed by varying the number of
# worker processes. Each benchmark point starts fresh workers, runs the
# simulation, then removes them to ensure clean isolation.
#
# Usage:
#   julia --project benchmark/run_scaling_distributed.jl

using CairoMakie
using Statistics
using Printf
using Distributed
using TestParticle
using StaticArrays

const N_PARTICLES = 16384
const N_SAMPLES = 3

# Worker counts to test
max_procs = Sys.CPU_THREADS
proc_counts = [2^i for i in 0:floor(Int, log2(max_procs))]
if !(max_procs in proc_counts)
    push!(proc_counts, max_procs)
end

println("Starting distributed scaling benchmark on worker counts: ", proc_counts)
println("Particles: $N_PARTICLES, samples per point: $N_SAMPLES")

uniform_B = x -> SA[0.0, 0.0, 1.0e-8]
uniform_E = x -> SA[0.0, 0.0, 0.0]
param = prepare(uniform_E, uniform_B; species = Proton)

x0 = [0.0, 0.0, 0.0]
v0 = [1.0e5, 0.0, 0.0]
stateinit = [x0..., v0...]
tspan = (0.0, 1.0e-5)
dt = 1.0e-9
savestepinterval = 10000

function make_prob(param, stateinit, tspan)
    function prob_func(prob, i, repeat)
        return remake(
            prob;
            u0 = [prob.u0[1], prob.u0[2], prob.u0[3], (i / 1000.0) * 1.0e5, 0.0, 0.0],
        )
    end
    return TraceProblem(stateinit, tspan, param; prob_func = prob_func)
end

times = Float64[]

for nw in proc_counts
    println("\n" * "-"^50)
    println("Benchmarking with $nw worker(s)...")

    pids = addprocs(nw)
    try
        @everywhere pids using TestParticle
        @everywhere pids using StaticArrays

        prob = make_prob(param, stateinit, tspan)

        # Warmup
        TestParticle.solve(
            prob, EnsembleDistributed(); trajectories = 10, dt, savestepinterval, batch_size = 1
        )

        # Timed samples
        sample_times = Float64[]
        for _ in 1:N_SAMPLES
            t = @elapsed TestParticle.solve(
                prob, EnsembleDistributed(); trajectories = N_PARTICLES, dt, savestepinterval, batch_size = max(1, N_PARTICLES ÷ nw)
            )
            push!(sample_times, t)
        end

        t_med = median(sample_times)
        push!(times, t_med)
        println(@sprintf("  Median time: %.3f s", t_med))
    finally
        rmprocs(pids)
    end
end

speedups = times[1] ./ times

# --- Plotting ---
println("\nGenerating plot...")
fig = Figure(; fontsize = 20)
ax = Axis(
    fig[1, 1];
    xlabel = "Number of Workers", ylabel = "Speedup",
    title = "Distributed Strong Scaling ($N_PARTICLES Particles)",
    xminorticksvisible = true, xminorticks = IntervalsBetween(5),
    yminorticksvisible = true, yminorticks = IntervalsBetween(5),
)

scatterlines!(ax, proc_counts, speedups; label = "Measured Speedup")
lines!(
    ax, proc_counts, Float64.(proc_counts);
    color = :black, linestyle = :dash, label = "Ideal Linear Scaling",
)

axislegend(ax; position = :lt)

plot_path = joinpath(@__DIR__, "distributed_scaling.png")
save(plot_path, fig)
println("Saved scaling plot to: ", plot_path)

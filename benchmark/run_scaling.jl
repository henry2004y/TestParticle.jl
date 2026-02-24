# Parallel Strong Scaling Driver

# This script runs `bench_parallel.jl` with different thread counts and plots the result
using CairoMakie
using Statistics
using Printf


# Thread counts to test (powers of 2 up to the system's logical cores)
max_threads = Sys.CPU_THREADS
threads_to_test = [2^i for i in 0:floor(Int, log2(max_threads))]
if !(max_threads in threads_to_test)
    push!(threads_to_test, max_threads)
end

println("Starting scaling benchmark on threads: ", threads_to_test)

times = Float64[]

for t in threads_to_test
    println("\n" * "-"^50)
    println("Benchmarking with $t threads...")

    # We write a small runner script to extract the median time of the threads benchmark
    runner_code = """
    using BenchmarkTools, TestParticle, StaticArrays, Printf, SciMLBase
    n_particles = 10_000
    uniform_B(x) = SA[0.0, 0.0, 1.0e-8]
    uniform_E(x) = SA[0.0, 0.0, 0.0]
    param = prepare(uniform_E, uniform_B; species = Proton)
    x0 = [0.0, 0.0, 0.0]; v0 = [1.0e5, 0.0, 0.0]; stateinit = [x0..., v0...]
    tspan = (0.0, 1.0e-5); dt = 1.0e-9
    prob_func(prob, i, repeat) = SciMLBase.remake(prob; u0 = [prob.u0[1:3]..., (i / 1000.0) * 1.0e5, 0.0, 0.0])
    prob_multi = TraceProblem(stateinit, tspan, param; prob_func = prob_func)

    # Warmup
    TestParticle.solve(prob_multi, EnsembleThreads(); trajectories = 10, dt=dt, savestepinterval=10000)

    bench_threads = @benchmark TestParticle.solve(\$prob_multi, EnsembleThreads(); trajectories = \$n_particles, dt = \$dt, savestepinterval = 10000) samples=5 seconds=30
    time_ms = median(bench_threads).time / 1.0e6
    println("RESULT_TIME_MS: \$time_ms")
    """

    # Write to temp file and execute
    tmp_file = tempname() * ".jl"
    write(tmp_file, runner_code)

    cmd = `julia --project=$(joinpath(@__DIR__, "..")) -t $t $tmp_file`
    output = read(cmd, String)
    println(output)

    # Parse output
    m = match(r"RESULT_TIME_MS:\s*([0-9.]+)", output)
    if m !== nothing
        push!(times, parse(Float64, m.captures[1]))
    else
        error("Could not parse result for $t threads")
    end
end

speedups = times[1] ./ times

# --- Plotting ---
println("\nGenerating plot...")
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Number of Threads", ylabel = "Speedup", title = "Strong Scaling (10,000 Particles)")

# Measured speedup
scatterlines!(ax, threads_to_test, speedups, color = :blue, label = "Measured Speedup")

# Ideal linear scaling
lines!(ax, threads_to_test, threads_to_test, color = :black, linestyle = :dash, label = "Ideal Linear Scaling")

axislegend(ax, position = :lt)

plot_path = joinpath(@__DIR__, "parallel_scaling.png")
save(plot_path, fig)
println("\nSaved scaling plot to: ", plot_path)

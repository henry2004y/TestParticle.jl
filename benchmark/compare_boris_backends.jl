# Performance Comparison: Boris CPU Serial vs CPU Threads vs GPU/KernelAbstractions
#
# Run locally with:
#   julia --project=test -t auto benchmark/compare_boris_backends.jl
#
# If you have an NVIDIA/AMD/Apple GPU, install the backend package (CUDA, AMDGPU, Metal)
# and run with that package loaded, e.g.:
#   julia --project=test -t auto -e "using CUDA; include(\"benchmark/compare_boris_backends.jl\")"

using TestParticle
using StaticArrays
using Statistics: median
import TestParticle: solve
using Printf
using KernelAbstractions: CPU
import KernelAbstractions as KA

const N_REPEATS = 3
const NSTEPS_ANALYTICAL = 1_000
const NSTEPS_GRID = 200
const ANALYTICAL_COUNTS = [1_000, 10_000, 100_000]
const GRID_COUNTS = [2_000, 10_000, 50_000]
# GPU-only. At these sizes the kernel, not the host, sets the time.
const LARGE_COUNTS = [100_000, 500_000, 1_000_000, 5_000_000]
# GPU-only. The sweep must reach where b*steps exceeds a, which is around 10^4 steps.
const KERNEL_SWEEP_N = 1_000_000
const KERNEL_SWEEP_STEPS = [100, 1_000, 10_000, 50_000, 200_000]

"""
    measure_time(f::Function, repeats::Int = N_REPEATS) -> Float64

Call `f` once as a warmup, then time it `repeats` times with a `GC.gc()` before each
run and return the median elapsed wall time in milliseconds.
"""
function measure_time(f::Function, repeats::Int = N_REPEATS)
    f()
    samples = Vector{Float64}(undef, repeats)
    for i in eachindex(samples)
        GC.gc()
        samples[i] = @elapsed f()
    end
    return median(samples) * 1000
end

function detect_gpu_backend()
    for (pkg, backend_expr) in [
            (:CUDA, "CUDA.CUDABackend()"),
            (:AMDGPU, "AMDGPU.ROCBackend()"),
            (:Metal, "Metal.MetalBackend()"),
            (:oneAPI, "oneAPI.oneAPIBackend()"),
        ]
        if isdefined(Main, pkg)
            try
                backend = eval(Meta.parse(backend_expr))
                return (string(pkg), backend)
            catch e
                @warn "Failed to initialize $pkg backend: $e"
            end
        end
    end
    return ("None", nothing)
end

function gpu_device_name()
    if isdefined(Main, :CUDA)
        try
            return string(CUDA.name(CUDA.device()))
        catch
        end
    end
    if isdefined(Main, :AMDGPU)
        try
            return string(AMDGPU.device())
        catch
        end
    end
    return "unknown"
end

# `solve` already synchronizes internally; the extra call guards against backends that
# only queue the trailing host-side copies.
function solve_gpu(prob, alg, backend, N, dt; kwargs...)
    sol = solve(prob, alg, backend; trajectories = N, dt, kwargs...)
    KA.synchronize(backend)
    return sol
end

function run_benchmark_analytical(
        prob64, prob32, N, dt64, dt32;
        saveat = (), save_everystep = false
    )
    # 1. CPU Serial (Float64)
    t_serial = measure_time() do
        solve(
            prob64, Boris(), EnsembleSerial();
            trajectories = N, dt = dt64, saveat, save_everystep
        )
        return nothing
    end

    # 2. CPU Threads (Float64)
    t_threads = measure_time() do
        solve(
            prob64, Boris(), EnsembleThreads();
            trajectories = N, dt = dt64, saveat, save_everystep
        )
        return nothing
    end

    # 3. KernelAbstractions CPU Batch (Float64)
    t_ka_cpu = measure_time() do
        solve(prob64, Boris(), CPU(); trajectories = N, dt = dt64, saveat, save_everystep)
        KA.synchronize(CPU())
        return nothing
    end

    # 4. GPU (Float64 and Float32)
    _, gpu_backend = detect_gpu_backend()
    t_gpu64 = nothing
    t_gpu32 = nothing
    if gpu_backend !== nothing
        try
            t_gpu64 = measure_time() do
                solve_gpu(prob64, Boris(), gpu_backend, N, dt64; saveat, save_everystep)
                return nothing
            end

            t_gpu32 = measure_time() do
                solve_gpu(prob32, Boris(), gpu_backend, N, dt32; saveat, save_everystep)
                return nothing
            end
        catch err
            @warn "GPU execution failed: $err"
        end
    end

    return (; N, t_serial, t_threads, t_ka_cpu, t_gpu64, t_gpu32)
end

function setup_grid_problem(; grid_res = 32, dt = 1.0f-8, nsteps = NSTEPS_GRID)
    x = range(0.0f0, 1.0f0, length = grid_res)
    y = range(0.0f0, 1.0f0, length = grid_res)
    z = range(0.0f0, 1.0f0, length = grid_res)

    B_data = [
        SA[0.0f0, 0.0f0, 1.0f-8] for _ in 1:grid_res, _ in 1:grid_res, _ in 1:grid_res
    ]
    E_data = [
        SA[0.0f0, 0.0f0, 0.0f0] for _ in 1:grid_res, _ in 1:grid_res, _ in 1:grid_res
    ]

    param32 = prepare(x, y, z, E_data, B_data; species = Proton, type = Float32)
    prob_func(prob, ctx) = remake(
        prob;
        u0 = Float32[
            0.1f0 + 0.8f0 * Float32(mod(ctx.sim_id, 100)) / 100.0f0,
            0.1f0 + 0.8f0 * Float32(mod(div(ctx.sim_id, 100), 100)) / 100.0f0,
            0.5f0,
            1.0f4, 0.0f0, 0.0f0,
        ]
    )
    tspan32 = (0.0f0, Float32(dt * nsteps))
    prob32 = TraceProblem(Float32[0.1, 0.1, 0.5, 1.0e4, 0, 0], tspan32, param32; prob_func)

    return prob32
end

function run_benchmark_grid(prob32, N, dt; saveat = (), save_everystep = false)
    # 1. CPU multithreading grid (one ODESolution per particle)
    t_cpu = measure_time() do
        solve(
            prob32, Boris(), EnsembleThreads();
            trajectories = N, dt, saveat, save_everystep
        )
        return nothing
    end

    # 2. KernelAbstractions CPU batch grid
    t_ka_cpu = measure_time() do
        solve(prob32, Boris(), CPU(); trajectories = N, dt, saveat, save_everystep)
        KA.synchronize(CPU())
        return nothing
    end

    _, gpu_backend = detect_gpu_backend()
    t_gpu = nothing
    t_gpu_sorted = nothing
    if gpu_backend !== nothing
        try
            t_gpu = measure_time() do
                solve_gpu(prob32, Boris(), gpu_backend, N, dt; saveat, save_everystep)
                return nothing
            end

            t_gpu_sorted = measure_time() do
                solve_gpu(
                    prob32, Boris(), gpu_backend, N, dt;
                    saveat, save_everystep, sort_particles = true
                )
                return nothing
            end
        catch err
            @warn "GPU Grid execution failed: $err"
        end
    end

    return (; N, t_cpu, t_ka_cpu, t_gpu, t_gpu_sorted)
end

"""
    initial_state_matrix(prob, N)

Build the `(N, 6)` matrix of initial states by calling the problem's own `prob_func`, so
the `u0` shortcut describes exactly the same ensemble as the default path.
"""
function initial_state_matrix(prob, N)
    states = Matrix{eltype(prob.u0)}(undef, N, 6)
    for i in 1:N
        states[i, :] .= prob.prob_func(prob, (sim_id = i, repeat = 1)).u0
    end
    return states
end

"""
    raw_kwargs(prob, N)

Keywords that switch a run to the configuration with the lowest host cost: a bulk `u0`
matrix instead of a per-particle `prob_func`, raw output instead of one `ODESolution` per
particle, and no saved initial state, which lets the raw output alias the state buffer.
The pusher kernel is untouched, so only the host side of the timing changes.
"""
raw_kwargs(prob, N) = (;
    u0 = initial_state_matrix(prob, N), raw_output = true, save_start = false,
)

"""
    run_benchmark_gpu_scaling(prob64, prob32, counts, dt64, dt32)

Time the GPU backend alone over `counts` particles. `nothing` is returned when no GPU
backend is available. A failure at some N (typically out of memory) stops the sweep and
keeps the points measured so far. With `raw = true` the same sweep is repeated with
`raw_kwargs`, which isolates how much of the wall time is host-side bookkeeping.
"""
function run_benchmark_gpu_scaling(
        prob64, prob32, counts, dt64, dt32;
        raw::Bool = false, saveat = (), save_everystep = false
    )
    _, gpu_backend = detect_gpu_backend()
    gpu_backend === nothing && return nothing

    Ns = Int[]
    times64 = Float64[]
    times32 = Float64[]
    times32raw = Float64[]
    for N in counts
        try
            kw = raw ? raw_kwargs(prob32, N) : (;)
            t64 = measure_time() do
                solve_gpu(prob64, Boris(), gpu_backend, N, dt64; saveat, save_everystep)
                return nothing
            end
            t32 = measure_time() do
                solve_gpu(prob32, Boris(), gpu_backend, N, dt32; saveat, save_everystep)
                return nothing
            end
            t32raw = measure_time() do
                solve_gpu(
                    prob32, Boris(), gpu_backend, N, dt32; saveat, save_everystep, kw...
                )
                return nothing
            end
            push!(Ns, N)
            push!(times64, t64)
            push!(times32, t32)
            push!(times32raw, t32raw)
        catch err
            @warn "GPU large-scale run stopped at $N particles" exception = err
            break
        end
    end
    return (; Ns, times64, times32, times32raw)
end

"""
    run_benchmark_kernel_sweep(prob64, prob32, N, dt64, dt32, steps)

Time the GPU backend at a fixed particle count while sweeping the number of Boris steps.
The host cost does not depend on the step count, so a fit of `t = a + b*steps` separates
the host cost `a` from the kernel `b`.
"""
function run_benchmark_kernel_sweep(
        prob64, prob32, N, dt64, dt32, steps;
        raw::Bool = false, saveat = (), save_everystep = false
    )
    _, gpu_backend = detect_gpu_backend()
    gpu_backend === nothing && return nothing

    kw64 = raw ? raw_kwargs(prob64, N) : (;)
    kw32 = raw ? raw_kwargs(prob32, N) : (;)

    used = Int[]
    times64 = Float64[]
    times32 = Float64[]
    times64raw = Float64[]
    times32raw = Float64[]
    for nt in steps
        try
            prob64_nt = remake(prob64; tspan = (zero(dt64), dt64 * nt))
            prob32_nt = remake(prob32; tspan = (zero(dt32), dt32 * nt))
            t64 = measure_time() do
                solve_gpu(prob64_nt, Boris(), gpu_backend, N, dt64; saveat, save_everystep)
                return nothing
            end
            t32 = measure_time() do
                solve_gpu(prob32_nt, Boris(), gpu_backend, N, dt32; saveat, save_everystep)
                return nothing
            end
            t64raw = measure_time() do
                solve_gpu(
                    prob64_nt, Boris(), gpu_backend, N, dt64;
                    saveat, save_everystep, kw64...
                )
                return nothing
            end
            t32raw = measure_time() do
                solve_gpu(
                    prob32_nt, Boris(), gpu_backend, N, dt32;
                    saveat, save_everystep, kw32...
                )
                return nothing
            end
            push!(used, nt)
            push!(times64, t64)
            push!(times32, t32)
            push!(times64raw, t64raw)
            push!(times32raw, t32raw)
        catch err
            @warn "Kernel sweep stopped at $nt steps" exception = err
            break
        end
    end
    return (; N, steps = used, times64, times32, times64raw, times32raw)
end

# Millions of particle-steps per second.
throughput(N, nsteps, t_ms) = N * nsteps / (t_ms * 1.0e3)

"""
    linear_fit(xs, ys) -> (; intercept, slope)

Least-squares fit of `y = intercept + slope * x`.
"""
function linear_fit(xs, ys)
    n = length(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    slope = sum((x - mx) * (y - my) for (x, y) in zip(xs, ys)) /
        sum((x - mx)^2 for x in xs)
    return (; intercept = my - slope * mx, slope)
end

"""
    print_speedup_table(ylabel, row_names, col_names, values, baselines)

Print `baselines[i] / values[i][j]` for every row and column. A `nothing` entry, either
in the values or in the baseline, prints `-`.
"""
function print_speedup_table(ylabel, row_names, col_names, values, baselines)
    width = 13
    @printf("%-9s", ylabel)
    for c in col_names
        @printf(" | %-*s", width, c)
    end
    println()
    println("-"^(10 + (width + 3) * length(col_names)))
    for (i, row) in enumerate(row_names)
        @printf("%-9s", row)
        for j in eachindex(col_names)
            v = values[i][j]
            s = (v === nothing || baselines[i] === nothing) ? "-" :
                @sprintf("%.2fx", baselines[i] / v)
            @printf(" | %-*s", width, s)
        end
        println()
    end
    return nothing
end

function main()
    println("="^88)
    println("Boris Pusher Benchmark: CPU Serial vs CPU Threads vs GPU Backends")
    gpu_name, _ = detect_gpu_backend()
    println("CPU Threads: ", Threads.nthreads())
    println("GPU Backend: ", gpu_name)
    if gpu_name != "None"
        println("GPU Device:  ", gpu_device_name())
    end
    @printf(
        "Timing: median of %d runs after one warmup; save_everystep = false\n", N_REPEATS
    )
    println("="^88)

    if Threads.nthreads() == 1
        @warn "Single Julia thread: restart with -t auto (or -t N) to measure " *
            "multithreaded CPU performance."
    end

    # 1. Analytical Field Benchmark
    println("\n[1] Analytical Uniform Field Tracing (", NSTEPS_ANALYTICAL, " steps)")
    println("-"^88)

    B64(x) = SA[0.0, 0.0, 1.0e-8]
    E64(x) = SA[0.0, 0.0, 0.0]
    param64 = prepare(E64, B64; species = Proton)
    dt64 = 1.0e-8
    tspan64 = (0.0, dt64 * NSTEPS_ANALYTICAL)
    prob_func64(prob, ctx) = remake(
        prob; u0 = [prob.u0[1:3]..., ctx.sim_id * 1.0e4, 0.0, 0.0]
    )
    prob64 = TraceProblem(
        [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan64, param64; prob_func = prob_func64
    )

    B32(x) = SA[0.0f0, 0.0f0, 1.0f-8]
    E32(x) = SA[0.0f0, 0.0f0, 0.0f0]
    param32 = prepare(E32, B32; species = Proton, type = Float32)
    dt32 = 1.0f-8
    tspan32 = (0.0f0, dt32 * NSTEPS_ANALYTICAL)
    prob_func32(prob, ctx) = remake(
        prob; u0 = Float32[prob.u0[1:3]..., Float32(ctx.sim_id) * 1.0f4, 0.0f0, 0.0f0]
    )
    prob32 = TraceProblem(
        Float32[0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan32, param32;
        prob_func = prob_func32
    )

    if gpu_name != "None"
        @printf(
            "%-9s | %-10s | %-10s | %-10s | %-10s | %-10s\n",
            "Particles", "CPU Ser", "CPU Thr", "KA CPU", "GPU (FP64)", "GPU (FP32)"
        )
    else
        @printf(
            "%-9s | %-13s | %-13s | %-13s\n",
            "Particles", "CPU Serial", "CPU Threads", "KA CPU"
        )
    end
    println("-"^88)

    analytical = Vector{NamedTuple}(undef, length(ANALYTICAL_COUNTS))
    for (k, N) in enumerate(ANALYTICAL_COUNTS)
        res = run_benchmark_analytical(prob64, prob32, N, dt64, dt32)
        analytical[k] = res
        if res.t_gpu32 !== nothing
            @printf(
                "%-9d | %7.2f ms | %7.2f ms | %7.2f ms | %7.2f ms | %7.2f ms\n",
                N, res.t_serial, res.t_threads, res.t_ka_cpu, res.t_gpu64, res.t_gpu32
            )
        else
            @printf(
                "%-9d | %8.2f ms  | %8.2f ms  | %8.2f ms\n",
                N, res.t_serial, res.t_threads, res.t_ka_cpu
            )
        end
    end

    println("\n  Speedup relative to CPU Serial:")
    print_speedup_table(
        "Particles",
        string.(ANALYTICAL_COUNTS),
        ["CPU Thr", "KA CPU", "GPU (FP64)", "GPU (FP32)"],
        [[r.t_threads, r.t_ka_cpu, r.t_gpu64, r.t_gpu32] for r in analytical],
        [r.t_serial for r in analytical]
    )

    let res = analytical[end], N = res.N
        println("\n  Effective throughput at $N particles ($(NSTEPS_ANALYTICAL) steps):")
        @printf(
            "    CPU serial %8.1f | CPU threads %8.1f | KA CPU %8.1f",
            throughput(N, NSTEPS_ANALYTICAL, res.t_serial),
            throughput(N, NSTEPS_ANALYTICAL, res.t_threads),
            throughput(N, NSTEPS_ANALYTICAL, res.t_ka_cpu)
        )
        if res.t_gpu32 !== nothing
            @printf(
                " | GPU FP64 %8.1f | GPU FP32 %8.1f",
                throughput(N, NSTEPS_ANALYTICAL, res.t_gpu64),
                throughput(N, NSTEPS_ANALYTICAL, res.t_gpu32)
            )
        end
        println(" M particle-steps/s")
    end

    # 2. 3D Grid Field Benchmark
    println("\n[2] 3D Grid Field Interpolation (32³ Grid, $(NSTEPS_GRID) steps, Float32)")
    println("-"^88)

    dt_grid = 1.0f-8
    grid_prob = setup_grid_problem(; dt = dt_grid)

    if gpu_name != "None"
        @printf(
            "%-9s | %-15s | %-15s | %-15s | %-15s\n",
            "Particles", "CPU Thr", "KA CPU", "GPU Unsort", "GPU Morton"
        )
    else
        @printf("%-9s | %-15s | %-15s\n", "Particles", "CPU Threads", "KA CPU")
    end
    println("-"^88)

    grid = Vector{NamedTuple}(undef, length(GRID_COUNTS))
    for (k, N) in enumerate(GRID_COUNTS)
        gres = run_benchmark_grid(grid_prob, N, dt_grid)
        grid[k] = gres
        if gres.t_gpu !== nothing
            @printf(
                "%-9d | %12.2f ms | %12.2f ms | %12.2f ms | %12.2f ms\n",
                N, gres.t_cpu, gres.t_ka_cpu, gres.t_gpu, gres.t_gpu_sorted
            )
        else
            @printf("%-9d | %12.2f ms | %12.2f ms\n", N, gres.t_cpu, gres.t_ka_cpu)
        end
    end

    println(
        "\n  `CPU Thr` builds one ODESolution per particle, which with grid fields\n" *
            "  costs far more than the pusher itself; `KA CPU` is the fair CPU baseline."
    )
    println("\n  Speedup relative to CPU Threads:")
    print_speedup_table(
        "Particles",
        string.(GRID_COUNTS),
        ["KA CPU", "GPU Unsort", "GPU Morton"],
        [[g.t_ka_cpu, g.t_gpu, g.t_gpu_sorted] for g in grid],
        [g.t_cpu for g in grid]
    )

    # 3. GPU-only large-scale scaling
    if gpu_name != "None"
        println("\n[3] GPU-only large-scale throughput (", NSTEPS_ANALYTICAL, " steps)")
        println("-"^88)
        @printf(
            "%-9s | %-13s | %-13s | %-13s | %-12s | %-12s\n",
            "Particles", "GPU (FP64)", "GPU (FP32)", "FP32 raw", "FP32 M/s", "Raw M/s"
        )
        println("-"^88)

        large = run_benchmark_gpu_scaling(
            prob64, prob32, LARGE_COUNTS, dt64, dt32; raw = true
        )
        for (i, N) in enumerate(large.Ns)
            @printf(
                "%-9d | %10.2f ms | %10.2f ms | %10.2f ms | %12.1f | %12.1f\n",
                N, large.times64[i], large.times32[i], large.times32raw[i],
                throughput(N, NSTEPS_ANALYTICAL, large.times32[i]),
                throughput(N, NSTEPS_ANALYTICAL, large.times32raw[i])
            )
        end

        println("\n  Speedup relative to GPU (FP64):")
        print_speedup_table(
            "Particles", string.(large.Ns), ["GPU (FP32)", "FP32 raw"],
            [[large.times32[i], large.times32raw[i]] for i in eachindex(large.Ns)],
            large.times64
        )
        println(
            "  `FP32 raw` uses a bulk `u0` matrix and `raw_output = true`, so no\n" *
                "  ODESolution is built. It changes the host bookkeeping only."
        )

        # 4. Kernel isolation: sweep the step count at fixed particle count.
        println(
            "\n[4] GPU-only kernel isolation ($KERNEL_SWEEP_N particles, step sweep)"
        )
        println("-"^88)
        @printf(
            "%-9s | %-13s | %-13s | %-13s | %-13s\n",
            "Steps", "GPU (FP64)", "GPU (FP32)", "FP64 raw", "FP32 raw"
        )
        println("-"^88)

        sweep = run_benchmark_kernel_sweep(
            prob64, prob32, KERNEL_SWEEP_N, dt64, dt32, KERNEL_SWEEP_STEPS; raw = true
        )
        for (i, nt) in enumerate(sweep.steps)
            @printf(
                "%-9d | %10.2f ms | %10.2f ms | %10.2f ms | %10.2f ms\n",
                nt, sweep.times64[i], sweep.times32[i],
                sweep.times64raw[i], sweep.times32raw[i]
            )
        end

        println("\n  Speedup relative to GPU (FP64):")
        print_speedup_table(
            "Steps", string.(sweep.steps), ["GPU (FP32)", "FP32 raw"],
            [[sweep.times32[i], sweep.times32raw[i]] for i in eachindex(sweep.steps)],
            sweep.times64
        )

        if length(sweep.steps) >= 2
            fits = (
                ("baseline", "FP64", linear_fit(Float64.(sweep.steps), sweep.times64)),
                ("baseline", "FP32", linear_fit(Float64.(sweep.steps), sweep.times32)),
                ("raw", "FP64", linear_fit(Float64.(sweep.steps), sweep.times64raw)),
                ("raw", "FP32", linear_fit(Float64.(sweep.steps), sweep.times32raw)),
            )
            println("\n  Fit of t = a + b*steps at $(sweep.N) particles")
            @printf(
                "    %-10s %-5s | %-13s | %-17s | %s\n",
                "variant", "prec", "host cost a", "kernel per step b", "kernel M steps/s"
            )
            println("    " * "-"^70)
            for (variant, prec, fit) in fits
                @printf(
                    "    %-10s %-5s | %10.2f ms | %14.4f ms | %15.1f\n",
                    variant, prec, fit.intercept, fit.slope,
                    sweep.N / (fit.slope * 1.0e3)
                )
            end
            println(
                "    `a` is the host cost (initialization and ODESolution assembly);\n" *
                    "    `b` is the pusher kernel alone. A negative `a` means the model\n" *
                    "    does not hold and the fit cannot be trusted. `b` must agree\n" *
                    "    between variants, since the raw path changes the host side only."
            )
        end
    end

    println("\n" * "="^88)
    println("Tips:")
    println(" - To run on NVIDIA GPUs:")
    println(
        "     julia --project=test -t auto -e " *
            "\"using CUDA; include(\\\"benchmark/compare_boris_backends.jl\\\")\""
    )
    println(" - To run with multiple CPU threads:")
    println("     julia --project=test -t auto benchmark/compare_boris_backends.jl")
    println(" - Sections [1] and [2] include host-side ODESolution assembly. Sections [3]")
    println("   and [4] show what that costs by repeating the GPU runs with a bulk")
    println("   `u0` matrix and `raw_output = true`.")
    return println("="^88)
end

if isempty(PROGRAM_FILE) || abspath(PROGRAM_FILE) == @__FILE__
    main()
end

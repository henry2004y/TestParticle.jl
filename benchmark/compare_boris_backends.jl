# Performance Comparison: Boris CPU Serial vs CPU Threads vs GPU/KernelAbstractions
#
# Run locally with:
#   julia --project=test benchmark/compare_boris_backends.jl
#
# If you have an NVIDIA/AMD/Apple GPU, install the backend package (CUDA, AMDGPU, Metal)
# and run with that package loaded, e.g.:
#   julia --project=test -e "using CUDA; include(\"benchmark/compare_boris_backends.jl\")"

using TestParticle
using StaticArrays
using KernelAbstractions
import TestParticle: solve
using Printf

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

function run_benchmark_analytical(
        prob64, prob32, N, dt64, dt32;
        saveat = (), save_everystep = false
    )
    # Warmup
    solve(
        prob64, Boris(), EnsembleSerial();
        trajectories = 1, dt = dt64, saveat, save_everystep
    )
    solve(
        prob64, Boris(), EnsembleThreads();
        trajectories = 1, dt = dt64, saveat, save_everystep
    )
    solve(prob64, Boris(), CPU(); trajectories = 1, dt = dt64, saveat, save_everystep)

    # 1. CPU Serial (Float64)
    GC.gc()
    b_cpu_serial = @timed solve(
        prob64, Boris(), EnsembleSerial();
        trajectories = N, dt = dt64, saveat, save_everystep
    )

    # 2. CPU Threads (Float64)
    GC.gc()
    b_cpu_threads = @timed solve(
        prob64, Boris(), EnsembleThreads();
        trajectories = N, dt = dt64, saveat, save_everystep
    )

    # 3. KernelAbstractions CPU Batch (Float64)
    GC.gc()
    b_ka_cpu = @timed solve(
        prob64, Boris(), CPU();
        trajectories = N, dt = dt64, saveat, save_everystep
    )

    # 4. GPU (Float64 and Float32)
    gpu_name, gpu_backend = detect_gpu_backend()
    t_gpu64 = nothing
    t_gpu32 = nothing
    if gpu_backend !== nothing
        try
            solve(
                prob64, Boris(), gpu_backend;
                trajectories = 1, dt = dt64, saveat, save_everystep
            )
            GC.gc()
            b_gpu64 = @timed solve(
                prob64, Boris(), gpu_backend;
                trajectories = N, dt = dt64, saveat, save_everystep
            )
            t_gpu64 = b_gpu64.time * 1000

            # Float32 GPU
            solve(
                prob32, Boris(), gpu_backend;
                trajectories = 1, dt = dt32, saveat, save_everystep
            )
            GC.gc()
            b_gpu32 = @timed solve(
                prob32, Boris(), gpu_backend;
                trajectories = N, dt = dt32, saveat, save_everystep
            )
            t_gpu32 = b_gpu32.time * 1000
        catch err
            @warn "GPU execution failed: $err"
        end
    end

    return (;
        N,
        t_serial = b_cpu_serial.time * 1000,
        t_threads = b_cpu_threads.time * 1000,
        t_ka_cpu = b_ka_cpu.time * 1000,
        gpu_name,
        t_gpu64,
        t_gpu32,
    )
end

function run_benchmark_grid(N, dt; grid_res = 32)
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
            1.0f4, 0.0f0, 0.0f0
        ]
    )
    tspan32 = (0.0f0, Float32(dt * 200))
    prob32 = TraceProblem(Float32[0.1, 0.1, 0.5, 1e4, 0, 0], tspan32, param32; prob_func)

    # CPU multithreading grid
    GC.gc()
    b_cpu = @timed solve(
        prob32, Boris(), EnsembleThreads();
        trajectories = N, dt
    )
    t_cpu = b_cpu.time * 1000

    gpu_name, gpu_backend = detect_gpu_backend()
    t_gpu = nothing
    t_gpu_sorted = nothing
    if gpu_backend !== nothing
        try
            solve(prob32, Boris(), gpu_backend; trajectories = 1, dt)
            GC.gc()
            b_gpu = @timed solve(prob32, Boris(), gpu_backend; trajectories = N, dt)
            t_gpu = b_gpu.time * 1000

            GC.gc()
            b_sorted = @timed solve(
                prob32, Boris(), gpu_backend;
                trajectories = N, dt, sort_particles = true
            )
            t_gpu_sorted = b_sorted.time * 1000
        catch err
            @warn "GPU Grid execution failed: $err"
        end
    end

    return (; N, t_cpu, gpu_name, t_gpu, t_gpu_sorted)
end

function main()
    println("="^88)
    println("Boris Pusher Benchmark: CPU Serial vs CPU Threads vs GPU Backends")
    gpu_name, _ = detect_gpu_backend()
    println("CPU Threads: ", Threads.nthreads())
    println("GPU Backend: ", gpu_name)
    println("="^88)

    # 1. Analytical Field Benchmark
    println("\n[1] Analytical Uniform Field Tracing (1,000 steps, save end state)")
    println("-"^88)

    B64(x) = SA[0.0, 0.0, 1.0e-8]
    E64(x) = SA[0.0, 0.0, 0.0]
    param64 = prepare(E64, B64; species = Proton)
    tspan64 = (0.0, 1.0e-5)
    dt64 = 1.0e-8
    prob_func64(prob, ctx) = remake(
        prob; u0 = [prob.u0[1:3]..., ctx.sim_id * 1.0e4, 0.0, 0.0]
    )
    prob64 = TraceProblem(
        [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan64, param64; prob_func = prob_func64
    )

    B32(x) = SA[0.0f0, 0.0f0, 1.0f-8]
    E32(x) = SA[0.0f0, 0.0f0, 0.0f0]
    param32 = prepare(E32, B32; species = Proton, type = Float32)
    tspan32 = (0.0f0, 1.0f-5)
    dt32 = 1.0f-8
    prob_func32(prob, ctx) = remake(
        prob; u0 = Float32[prob.u0[1:3]..., Float32(ctx.sim_id) * 1.0f4, 0.0f0, 0.0f0]
    )
    prob32 = TraceProblem(
        Float32[0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan32, param32;
        prob_func = prob_func32
    )

    analytical_counts = [100, 1000, 5000, 10000]

    if gpu_name != "None"
        @printf(
            "%-7s | %-10s | %-10s | %-10s | %-10s | %-10s | %-8s\n",
            "N", "CPU Ser", "CPU Thr", "KA CPU", "GPU (FP64)", "GPU (FP32)", "F32/Ser"
        )
        println("-"^88)
    else
        @printf(
            "%-7s | %-12s | %-12s | %-12s | %-12s\n",
            "N", "CPU Serial", "CPU Threads", "KA CPU", "KA/Serial"
        )
        println("-"^70)
    end

    for N in analytical_counts
        res = run_benchmark_analytical(prob64, prob32, N, dt64, dt32)
        if res.t_gpu32 !== nothing
            speedup_str = @sprintf("%.1fx", res.t_serial / res.t_gpu32)
            @printf(
                "%-7d | %7.2f ms | %7.2f ms | %7.2f ms | %7.2f ms  | %7.2f ms  | %-8s\n",
                N, res.t_serial, res.t_threads, res.t_ka_cpu,
                res.t_gpu64, res.t_gpu32, speedup_str
            )
        else
            speedup_str = @sprintf("%.2fx", res.t_serial / res.t_ka_cpu)
            @printf(
                "%-7d | %8.2f ms  | %8.2f ms  | %8.2f ms  | %-12s\n",
                N, res.t_serial, res.t_threads, res.t_ka_cpu, speedup_str
            )
        end
    end

    # 2. 3D Grid Field Benchmark
    println("\n[2] 3D Grid Field Interpolation (32³ Grid, 200 steps, Float32)")
    println("-"^88)
    grid_counts = [500, 2000, 10000]

    if gpu_name != "None"
        @printf(
            "%-7s | %-14s | %-14s | %-14s | %-10s\n",
            "N", "CPU Threads", "GPU Unsorted", "GPU Morton", "Speedup"
        )
        println("-"^88)
        for N in grid_counts
            gres = run_benchmark_grid(N, 1.0f-8)
            speedup_str = @sprintf("%.1fx", gres.t_cpu / gres.t_gpu_sorted)
            @printf(
                "%-7d | %9.2f ms   | %9.2f ms   | %9.2f ms   | %-10s\n",
                N, gres.t_cpu, gres.t_gpu, gres.t_gpu_sorted, speedup_str
            )
        end
    else
        @printf("%-7s | %-14s\n", "N", "CPU Threads")
        println("-"^30)
        for N in grid_counts
            gres = run_benchmark_grid(N, 1.0f-8)
            @printf("%-7d | %9.2f ms\n", N, gres.t_cpu)
        end
    end

    println("\n" * "="^88)
    println("Tips:")
    println(" - To run on NVIDIA GPUs:")
    println(
        "     julia --project=test -e " *
            "\"using CUDA; include(\\\"benchmark/compare_boris_backends.jl\\\")\""
    )
    println(" - To run with multiple CPU threads:")
    println("     julia --project=test -t auto benchmark/compare_boris_backends.jl")
    return println("="^88)
end

if isempty(PROGRAM_FILE) || abspath(PROGRAM_FILE) == @__FILE__
    main()
end

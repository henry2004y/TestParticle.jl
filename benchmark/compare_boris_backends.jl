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

function run_benchmark(
        prob, N, dt, mode_desc;
        saveat = (), save_everystep = false
    )
    # Warmup / JIT
    solve(
        prob, Boris(), EnsembleSerial();
        trajectories = 1, dt, saveat, save_everystep
    )
    solve(
        prob, Boris(), EnsembleThreads();
        trajectories = 1, dt, saveat, save_everystep
    )
    solve(
        prob, Boris(), CPU();
        trajectories = 1, dt, saveat, save_everystep
    )

    # 1. CPU Serial
    GC.gc()
    b_cpu_serial = @timed solve(
        prob, Boris(), EnsembleSerial();
        trajectories = N, dt, saveat, save_everystep
    )

    # 2. CPU Threads
    GC.gc()
    b_cpu_threads = @timed solve(
        prob, Boris(), EnsembleThreads();
        trajectories = N, dt, saveat, save_everystep
    )

    # 3. KernelAbstractions CPU Batch
    GC.gc()
    b_ka_cpu = @timed solve(
        prob, Boris(), CPU();
        trajectories = N, dt, saveat, save_everystep
    )

    # 4. GPU (if loaded)
    gpu_name, gpu_backend = detect_gpu_backend()
    b_gpu = if gpu_backend !== nothing
        try
            # Warmup
            solve(
                prob, Boris(), gpu_backend;
                trajectories = 1, dt, saveat, save_everystep
            )
            GC.gc()
            @timed solve(
                prob, Boris(), gpu_backend;
                trajectories = N, dt, saveat, save_everystep
            )
        catch err
            @warn "GPU execution failed: $err"
            nothing
        end
    else
        nothing
    end

    return (;
        N,
        t_serial = b_cpu_serial.time * 1000,
        m_serial = b_cpu_serial.bytes / 1024,
        t_threads = b_cpu_threads.time * 1000,
        m_threads = b_cpu_threads.bytes / 1024,
        t_ka_cpu = b_ka_cpu.time * 1000,
        m_ka_cpu = b_ka_cpu.bytes / 1024,
        gpu_name,
        t_gpu = b_gpu === nothing ? nothing : b_gpu.time * 1000,
        m_gpu = b_gpu === nothing ? nothing : b_gpu.bytes / 1024,
    )
end

function main()
    println("="^88)
    println("Boris Pusher Benchmark: CPU Serial vs CPU Threads vs KernelAbstractions")
    gpu_name, _ = detect_gpu_backend()
    println("CPU Threads: ", Threads.nthreads())
    println("GPU Backend: ", gpu_name)
    println("="^88)

    # Problem setup: uniform B along z, zero E
    B(x) = SA[0.0, 0.0, 1.0e-8]
    E(x) = SA[0.0, 0.0, 0.0]
    param = prepare(E, B; species = Proton)

    # 1,000 steps
    tspan = (0.0, 1.0e-5)
    dt = 1.0e-8
    saveat_interval = 1.0e-7 # 100 steps interval

    # Ensemble initial perturbation
    prob_func(prob, ctx) = remake(
        prob; u0 = [prob.u0[1:3]..., ctx.sim_id * 1.0e4, 0.0, 0.0]
    )
    prob = TraceProblem(
        [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan, param; prob_func
    )

    trajectory_counts = [10, 100, 1000, 5000]

    for (scen_idx, (scen_name, saveat, save_everystep)) in enumerate([
        ("Scenario A: Save End State Only (Pusher Compute)", (), false),
        ("Scenario B: Save Every 100 Steps (Dense Output)", saveat_interval, false),
    ])
        println("\n" * "-"^88)
        println(" $scen_name")
        println("-"^88)

        if gpu_name != "None"
            @printf(
                "%-7s | %-12s | %-12s | %-12s | %-12s | %-10s\n",
                "N", "CPU Serial", "CPU Threads", "KA CPU", "GPU ($gpu_name)", "Speedup"
            )
            println("-"^88)
        else
            @printf(
                "%-7s | %-12s | %-12s | %-12s | %-12s\n",
                "N", "CPU Serial", "CPU Threads", "KA CPU", "KA/Serial"
            )
            println("-"^70)
        end

        for N in trajectory_counts
            res = run_benchmark(
                prob, N, dt, scen_name;
                saveat, save_everystep
            )

            if res.t_gpu !== nothing
                speedup_str = @sprintf("%.2fx", res.t_serial / res.t_gpu)
                @printf(
                    "%-7d | %8.2f ms  | %8.2f ms  | %8.2f ms  | %8.2f ms  | %-10s\n",
                    N, res.t_serial, res.t_threads, res.t_ka_cpu, res.t_gpu, speedup_str
                )
            else
                speedup_str = @sprintf("%.2fx", res.t_serial / res.t_ka_cpu)
                @printf(
                    "%-7d | %8.2f ms  | %8.2f ms  | %8.2f ms  | %-12s\n",
                    N, res.t_serial, res.t_threads, res.t_ka_cpu, speedup_str
                )
            end
        end
    end

    println("\n" * "="^88)
    println("Tips:")
    println(" - To run on NVIDIA GPUs:")
    println("     julia --project=test -e " *
            "\"using CUDA; include(\\\"benchmark/compare_boris_backends.jl\\\")\"")
    println(" - To run with multiple CPU threads:")
    println("     julia --project=test -t auto benchmark/compare_boris_backends.jl")
    println("="^88)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

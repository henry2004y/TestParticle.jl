# Precompiling workloads

@setup_workload begin
    @compile_workload begin
        # numerical field parameters
        x = range(-10, 10, length = 4)
        y = range(-10, 10, length = 6)
        z = range(-10, 10, length = 8)
        B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
        E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]

        B[3, :, :, :] .= 10.0e-9
        E[3, :, :, :] .= 5.0e-10

        mesh = CartesianGrid(
            (first(x), first(y), first(z)), (last(x), last(y), last(z));
            dims = (length(x) - 1, length(y) - 1, length(z) - 1)
        )

        vdf = Maxwellian([0.0, 0.0, 0.0], 1.0e-9, 1.0e6)
        v = rand(vdf)
        vdf = BiMaxwellian([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0e-9, 1.0e-9, 1.0e6)
        v = rand(vdf)
        # numerical field
        param = prepare(x, y, z, E, B)
        param = prepare(mesh, E, B)
        # analytical field
        param = prepare(DipoleField())
        # Boris pusher
        stateinit = let x0 = [0.0, 0.0, 0.0], v0 = [0.0, 1.0e5, 0.0]
            [x0..., v0...]
        end
        tspan = (0.0, 1.0)
        dt = 0.5
        prob = TraceProblem(stateinit, tspan, param)
        sol = solve(prob; dt, savestepinterval = 100)
        sol = solve(prob, EnsembleThreads(); dt, savestepinterval = 100)

        # Adaptive Boris
        alg_adaptive = AdaptiveBoris(dtmax = 1.0)
        sol_adaptive = solve(prob, alg_adaptive)

        sol_adaptive_sf = solve(prob, alg_adaptive; save_fields = true, save_work = true)
        E, B = get_fields(sol_adaptive)
        work = get_work(sol_adaptive)

        # Multistep Boris
        sol_multistep = solve(prob; dt, n = 2)
        sol_multistep_sf = solve(prob; dt, n = 2, save_fields = true)
        # guiding center
        gc_init = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        stateinit_gc, param_gc = prepare_gc(gc_init, ZeroField(), DipoleField())
        stateinit_gc, param_gc = prepare_gc(gc_init, x, y, z, E, B)
        # native GC solvers
        prob_native = TraceGCProblem(stateinit_gc, tspan, param_gc)
        sol_native_rk4 = solve(prob_native; dt = 0.1, alg = :rk4)
        sol_native_rk45 = solve(prob_native; alg = :rk45)

        # hybrid solver
        prob_hybrid = TraceHybridProblem(stateinit, tspan, param)
        alg_hybrid = AdaptiveHybrid(threshold = 0.1, dtmax = 1.0)
        sol_hybrid = solve(prob_hybrid, alg_hybrid)
    end
end

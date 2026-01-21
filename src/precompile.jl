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
        # guiding center
        gc_init = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        stateinit_gc, param_gc = prepare_gc(gc_init, ZeroField(), DipoleField())
        stateinit_gc, param_gc = prepare_gc(gc_init, x, y, z, E, B)
        # native GC
        prob_gc = TraceGCProblem(stateinit_gc, tspan, param_gc)
        sol = solve(prob_gc; dt, alg = :rk4)
        sol = solve(prob_gc; dt, alg = :rk45)
    end
end

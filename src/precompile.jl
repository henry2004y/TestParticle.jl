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

        # numerical field
        param = prepare(x, y, z, E, B)
        # analytical field
        B_analytic(r) = SA[0.0, 0.0, 1.0e-9]
        param = prepare(B_analytic)
        stateinit = let x0 = [0.0, 0.0, 0.0], v0 = [0.0, 1.0e5, 0.0]
            [x0..., v0...]
        end
        tspan = (0.0, 1.0)
        prob = TraceProblem(stateinit, tspan, param)

        # guiding center
        gc_init = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        stateinit_gc, param_gc = prepare_gc(gc_init; B = B_analytic)
        stateinit_gc, param_gc = prepare_gc(gc_init, x, y, z, E, B)
        # native GC solvers
        prob_native = TraceGCProblem(stateinit_gc, tspan, param_gc)
        sol_native_rk4 = solve(prob_native; dt = 0.1, alg = :rk4)
        sol_native_rk45 = solve(prob_native; alg = :rk45)

        # explicit trace functions
        t = 0.0
        out6, out4 = zeros(6), zeros(4)
        # non-relativistic
        trace(stateinit, param, t)
        trace!(out6, stateinit, param, t)
        # relativistic
        trace_relativistic(stateinit, param, t)
        trace_relativistic!(out6, stateinit, param, t)
        # guiding center
        trace_gc(stateinit_gc, param_gc, t)
        trace_gc!(out4, stateinit_gc, param_gc, t)
    end
end

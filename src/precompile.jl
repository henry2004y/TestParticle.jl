# Precompiling workloads
#
# Analytic fields are user functions, each its own type, so code specialized on one
# defined here is never reused; only interpolated fields are worth caching.

@setup_workload begin
    x = range(-10, 10, length = 4)
    y = range(-10, 10, length = 6)
    z = range(-10, 10, length = 8)
    B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
    E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]
    B[3, :, :, :] .= 10.0e-9
    E[3, :, :, :] .= 5.0e-10

    stateinit = [0.0, 0.0, 0.0, 0.0, 1.0e5, 0.0]
    stateinit_s = SVector{6}(stateinit)
    gc_init = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    tspan = (0.0, 1.0)
    t = 0.0
    out6, out4 = zeros(6), zeros(4)

    @compile_workload begin
        param = prepare(x, y, z; E, B)
        prob = TraceProblem(stateinit, tspan, param)
        trace(stateinit_s, param, t)
        trace!(out6, stateinit, param, t)
        trace_relativistic(stateinit_s, param, t)
        trace_relativistic!(out6, stateinit, param, t)

        stateinit_gc, param_gc = prepare_gc(gc_init, x, y, z; E, B)
        prob_gc = TraceGCProblem(stateinit_gc, tspan, param_gc)
        solve(prob_gc, GCRK4(); dt = 0.1)
        solve(prob_gc, GCRK45())
        trace_gc(stateinit_gc, param_gc, t)
        trace_gc!(out4, stateinit_gc, param_gc, t)
    end
end

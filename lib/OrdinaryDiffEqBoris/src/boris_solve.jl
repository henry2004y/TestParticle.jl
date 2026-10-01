# Fast-path solve! for fixed-step AbstractBoris on ODEIntegrator.

@inline function _can_fastpath_boris(integrator::ODEIntegrator{<:AbstractBoris})
    opts = integrator.opts
    t0, t1 = integrator.sol.prob.tspan
    dt = integrator.dt
    iszero(dt) && return false

    !opts.adaptive || return false
    isempty(opts.callback.continuous_callbacks) || return false
    isempty(opts.callback.discrete_callbacks) || return false
    opts.isoutofdomain === ODE_DEFAULT_ISOUTOFDOMAIN ||
        return false
    opts.save_idxs === nothing || return false
    isempty(opts.saveat) || return false
    length(opts.tstops) <= 1 || return false
    if !isempty(opts.tstops)
        first(opts.tstops) == integrator.tdir * t1 || return false
    end

    nt = round(Int, (t1 - t0) / dt)
    nt > 0 || return false
    isapprox(t0 + nt * dt, t1; atol = 100 * eps(typeof(dt))) || return false

    return true
end

function SciMLBase.solve!(integrator::ODEIntegrator{<:AbstractBoris})
    if !_can_fastpath_boris(integrator)
        return invoke(SciMLBase.solve!, Tuple{ODEIntegrator}, integrator)
    end

    dt = integrator.dt
    t0, t1 = integrator.sol.prob.tspan
    p = integrator.p
    alg = integrator.alg
    u = integrator.u
    opts = integrator.opts
    r = SVector(u[1], u[2], u[3])

    nt = round(Int, (t1 - t0) / dt)

    cache = integrator.cache
    v_half = cache.v_half
    fields = cache.fields

    save_everystep = opts.save_everystep && opts.save_on

    if !save_everystep
        t = t0
        r, v_half = advance_boris(v_half, r, dt, t, fields, alg)
        t += dt
        for _ in 2:nt
            r, v_half = advance_boris(v_half, r, dt, t, p, alg)
            t += dt
        end
        v_node = update_velocity_node(v_half, r, dt, t, p, alg)
        u_final = u isa SVector ? vcat(r, v_node) :
            typeof(u)([r[1], r[2], r[3], v_node[1], v_node[2], v_node[3]])
        integrator.u = u_final
        integrator.t = t
        if opts.save_end
            push!(integrator.sol.u, u_final)
            push!(integrator.sol.t, t)
        end
    else
        sol_u = integrator.sol.u
        sol_t = integrator.sol.t
        idx = length(sol_u)
        nout = idx + nt
        resize!(sol_u, nout)
        resize!(sol_t, nout)
        t = t0
        for _ in 1:nt
            r, v_half = advance_boris(v_half, r, dt, t, fields, alg)
            t += dt
            fields = _node_fields(p, r, t)
            v_node = update_velocity_node(v_half, r, dt, t, fields, alg)
            u_node = u isa SVector ? vcat(r, v_node) :
                typeof(u)([r[1], r[2], r[3], v_node[1], v_node[2], v_node[3]])
            idx += 1
            sol_u[idx] = u_node
            sol_t[idx] = t
        end
        integrator.u = sol_u[end]
        integrator.t = t
        cache.fields = fields
    end

    integrator.sol.stats.naccept = nt
    integrator.sol = SciMLBase.solution_new_retcode(integrator.sol, ReturnCode.Success)
    return integrator.sol
end

# Boris particle pusher

get_EField(p::AbstractODEProblem) = get_EField(p.p)
get_BField(p::AbstractODEProblem) = get_BField(p.p)

get_BField(sol::AbstractODESolution) = get_BField(sol.prob)
get_EField(sol::AbstractODESolution) = get_EField(sol.prob)

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false

"""
    boris_velocity_update(v, E, B, qdt_2m)

Update velocity using the Boris method, returning the new velocity as an SVector.
This is the core logic shared between the standard solver and the kernel solver.
"""
@inline @muladd function boris_velocity_update(v, E, B, qdt_2m)
    t_rotate = qdt_2m * B
    t_mag2 = sum(abs2, t_rotate)
    s_rotate = 2 * t_rotate / (1 + t_mag2)

    v⁻ = v + qdt_2m * E
    v′ = v⁻ + (v⁻ × t_rotate)
    v⁺ = v⁻ + (v′ × s_rotate)

    v_new = v⁺ + qdt_2m * E

    return v_new
end

"""
    update_velocity(v, r, dt, t, param)

Update velocity using the Boris method, returning the new velocity as an SVector.
"""
@inline @muladd function update_velocity(v, r, dt, t, param)
    q2m, _, Efunc, Bfunc, _ = param
    E = Efunc(r, t)
    B = Bfunc(r, t)
    qdt_2m = q2m * 0.5 * dt

    return boris_velocity_update(v, E, B, qdt_2m)
end

function _default_batch_size(ensemblealg, trajectories)
    if ensemblealg isa EnsembleDistributed || ensemblealg isa EnsembleSplitThreads
        return max(1, trajectories ÷ nworkers())
    end
    return 1
end

function _get_sol_type(prob, dt, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    u0 = prob.u0
    tspan = prob.tspan
    T_t = typeof(tspan[1] + dt)
    T = eltype(u0)

    n_vars = 6
    if SaveFields
        n_vars += 6
    end
    if SaveWork
        n_vars += 4
    end

    u = SVector{n_vars, T}[]
    interp = LinearInterpolation(T_t[], u)
    alg = :boris

    sol = build_solution(prob, alg, T_t[], u; interp = interp)
    return typeof(sol)
end

@inline function _prepare_saved_data(xv, p, t, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    data = xv

    # Pre-declare variables to share between blocks if both are true
    local E_field, magnetic_props

    if SaveFields
        r = get_x(xv)
        T = eltype(xv)
        E_field = SVector{3, T}(get_EField(p)(r, t))

        # We need magnetic properties for work, so if SaveWork is also true, compute them now
        if SaveWork
            q2m, m, Efunc, Bfunc, _ = p
            # get_magnetic_properties returns (B, ∇B, κ, b̂, Bmag)
            magnetic_props = get_magnetic_properties(r, t, Bfunc)
            B_vec = SVector{3, T}(magnetic_props[1])
            data = vcat(data, E_field, B_vec)
        else
            B_vec = SVector{3, T}(get_BField(p)(r, t))
            data = vcat(data, E_field, B_vec)
        end
    end

    if SaveWork
        # If SaveFields was true, we already computed E_field and magnetic_props
        if SaveFields
            work = get_work_rates(xv, p, t, magnetic_props, E_field)
        else
            work = get_work_rates(xv, p, t)
        end
        data = vcat(data, work)
    end

    return data
end

"""
    get_fields(sol::AbstractODESolution)

Return the electric and magnetic fields from the solution `sol`.
"""
function get_fields(sol::AbstractODESolution)
    Efunc, Bfunc = _get_field_funcs(sol.prob)

    E = map((u, t) -> Efunc(get_x(u), t), sol.u, sol.t)
    B = map((u, t) -> Bfunc(get_x(u), t), sol.u, sol.t)

    return E, B
end

function _get_field_funcs(prob::TraceGCProblem)
    # p = (q, q2m, μ, Efunc, Bfunc)
    p = prob.p
    return p[4], p[5]
end

function _get_field_funcs(prob)
    p = prob.p
    return get_EField(p), get_BField(p)
end

"""
    get_work(sol::AbstractODESolution)

Return the work done by the electric field from the solution `sol`.
"""
function get_work(sol::AbstractODESolution)
    return _get_work(sol, sol.prob)
end

function _get_work(sol, prob::TraceGCProblem)
    p = prob.p
    return map((u, t) -> get_work_rates_gc(u, p, t), sol.u, sol.t)
end

function _get_work(sol, prob)
    p = prob.p
    return map((u, t) -> get_work_rates(u, p, t), sol.u, sol.t)
end

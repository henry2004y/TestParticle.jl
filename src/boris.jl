# Native particle pusher

struct TraceProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <:
    AbstractODEProblem{uType, tType, isinplace}
    f::F
    "initial condition"
    u0::uType
    "time span"
    tspan::tType
    "(q2m, m, E, B, F)"
    p::P
    "function for setting initial conditions"
    prob_func::PF
end

const TN_MAG_THRESHOLD = 1.0e-4

get_EField(p::AbstractODEProblem) = get_EField(p.p)
get_BField(p::AbstractODEProblem) = get_BField(p.p)

get_BField(sol::AbstractODESolution) = get_BField(sol.prob)
get_EField(sol::AbstractODESolution) = get_EField(sol.prob)

DEFAULT_PROB_FUNC(prob, i, repeat) = prob

function TraceProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    _f = ODEFunction{true, DEFAULT_SPECIALIZATION}(x -> nothing) # dummy func
    return TraceProblem{
        typeof(u0), typeof(tspan), true, typeof(p), typeof(_f),
        typeof(prob_func),
    }(
        _f,
        u0,
        tspan,
        p,
        prob_func
    )
end
# For remake
function TraceProblem{iip}(; f, u0, tspan, p, prob_func) where {iip}
    return TraceProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f),
        typeof(prob_func),
    }(
        f,
        u0,
        tspan,
        p,
        prob_func
    )
end

struct BorisMethod{TV}
    # intermediate variables used in the solver
    v⁻::TV
    v′::TV
    v⁺::TV
    t_rotate::TV
    s_rotate::TV
    v⁻_cross_t::TV
    v′_cross_s::TV
end

function BorisMethod(T::Type{<:AbstractFloat} = Float64)
    v⁻ = MVector{3, T}(undef)
    v′ = MVector{3, T}(undef)
    v⁺ = MVector{3, T}(undef)
    t_rotate = MVector{3, T}(undef)
    s_rotate = MVector{3, T}(undef)
    v⁻_cross_t = MVector{3, T}(undef)
    v′_cross_s = MVector{3, T}(undef)

    return BorisMethod{typeof(v⁻)}(v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s)
end

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

    v_minus = v + qdt_2m * E
    v_prime = v_minus + (v_minus × t_rotate)
    v_plus = v_minus + (v_prime × s_rotate)

    v_new = v_plus + qdt_2m * E

    return v_new
end

"""
    update_velocity(v, r, param, dt, t)

Update velocity using the Boris method, returning the new velocity as an SVector.
"""
@inline @muladd function update_velocity(v, r, param::P, dt, t) where {P}
    q2m, _, Efunc, Bfunc, _ = param
    E = Efunc(r, t)
    B = Bfunc(r, t)
    qdt_2m = q2m * 0.5 * dt

    return boris_velocity_update(v, E, B, qdt_2m)
end

"""
    update_velocity!(xv, paramBoris, param, dt, t)

Update velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation.
Reference: [DTIC](https://apps.dtic.mil/sti/citations/ADA023511)
"""
@muladd function update_velocity!(xv, paramBoris, param, dt, t)
    (; v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s) = paramBoris
    q2m, _, Efunc, Bfunc, _ = param
    E = Efunc(xv, t)
    B = Bfunc(xv, t)
    # t vector
    for dim in 1:3
        t_rotate[dim] = q2m * B[dim] * 0.5 * dt
    end
    t_mag2 = sum(abs2, t_rotate)
    # s vector
    for dim in 1:3
        s_rotate[dim] = 2 * t_rotate[dim] / (1 + t_mag2)
    end
    # v-
    for dim in 1:3
        v⁻[dim] = xv[dim + 3] + q2m * E[dim] * 0.5 * dt
    end
    # v′
    cross!(v⁻, t_rotate, v⁻_cross_t)
    for dim in 1:3
        v′[dim] = v⁻[dim] + v⁻_cross_t[dim]
    end
    # v+
    cross!(v′, s_rotate, v′_cross_s)
    for dim in 1:3
        v⁺[dim] = v⁻[dim] + v′_cross_s[dim]
    end
    # v[n+1/2]
    for dim in 1:3
        xv[dim + 3] = v⁺[dim] + q2m * E[dim] * 0.5 * dt
    end

    return
end

"""
Update location in one timestep `dt`.
"""
@muladd function update_location!(xv, dt)
    xv[1] = xv[1] + xv[4] * dt
    xv[2] = xv[2] + xv[5] * dt
    xv[3] = xv[3] + xv[6] * dt

    return
end

"""
In-place cross product.
"""
@muladd function cross!(v1, v2, vout)
    vout[1] = v1[2] * v2[3] - v1[3] * v2[2]
    vout[2] = v1[3] * v2[1] - v1[1] * v2[3]
    vout[3] = v1[1] * v2[2] - v1[2] * v2[1]

    return
end

"""
    solve(prob::TraceProblem; trajectories::Int=1, dt::AbstractFloat,
        savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN,
        n::Int=1, save_start::Bool=true, save_end::Bool=true, save_everystep::Bool=true,
        save_fields::Bool=false, save_work::Bool=false)

Trace particles using the Boris method with specified `prob`.

# keywords

  - `trajectories::Int`: number of trajectories to trace.
  - `dt::AbstractFloat`: time step.
  - `savestepinterval::Int`: saving output interval.
  - `isoutofdomain::Function`: a function with input of position and velocity vector `xv` that determines whether to stop tracing.
  - `n::Int`: number of substeps for the Multistep Boris method. Default is 1 (standard Boris).
  - `save_start::Bool`: save the initial condition. Default is `true`.
  - `save_end::Bool`: save the final condition. Default is `true`.
  - `save_everystep::Bool`: save the state at every `savestepinterval`. Default is `true`.
  - `save_fields::Bool`: save the electric and magnetic fields. Default is `false`.
  - `save_work::Bool`: save the work done by the electric field. Default is `false`.
  - `batch_size::Int`: the number of trajectories to process per worker in `EnsembleDistributed`. Default is `max(1, trajectories ÷ Threads.nthreads())`.
"""
@inline function solve(
        prob::TraceProblem, ensemblealg::EA = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1, dt::AbstractFloat,
        isoutofdomain::F = ODE_DEFAULT_ISOUTOFDOMAIN, n::Int = 1,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        save_fields::Bool = false, save_work::Bool = false, maxiters::Int = 1_000_000,
        batch_size::Int = max(1, trajectories ÷ Threads.nthreads())
    ) where {EA <: BasicEnsembleAlgorithm, F}
    return _solve(
        ensemblealg, prob, trajectories, dt, savestepinterval, isoutofdomain, n,
        save_start, save_end, save_everystep, Val(save_fields), Val(save_work), maxiters,
        batch_size
    )
end

function _dispatch_boris!(
        sols, prob::TraceProblem, irange, savestepinterval, dt, nt, nout, isoutofdomain::F, n,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork, F}
    return if n == 1
        _boris!(
            sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
            save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
        )
    else
        _multistep_boris!(
            sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
            save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
        )
    end
end

@inline function _solve(
        ::EnsembleSerial, prob::TraceProblem, trajectories, dt, savestepinterval, isoutofdomain::F, n,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}, maxiters,
        batch_size
    ) where {SaveFields, SaveWork, F}
    sols, nt,
        nout = _prepare(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork), maxiters
    )
    irange = 1:trajectories
    _dispatch_boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
    )

    return sols
end

@inline function _solve(
        ::EnsembleThreads, prob::TraceProblem, trajectories, dt, savestepinterval, isoutofdomain::F, n,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}, maxiters,
        batch_size
    ) where {SaveFields, SaveWork, F}
    sols, nt,
        nout = _prepare(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork), maxiters
    )

    nchunks = Threads.nthreads()
    Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        _dispatch_boris!(
            sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
            save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
        )
    end

    return sols
end

"""
    _solve_single_boris(prob, i, ...)

Solve a single trajectory `i` of `prob` for use with `EnsembleDistributed`.

`_generic_boris!` uses the loop index `i` for two purposes simultaneously:
applying `prob_func(prob, i, false)` to select per-particle initial conditions,
and storing the result at `sols[i]`. For a distributed worker handling only one
trajectory at a time, a 1-element `local_sols` would be out of bounds if `i > 1`
were passed directly. Decoupling the two uses would require threading a storage
offset through the entire `_dispatch_boris!` → `_generic_boris!` call chain.

Instead, we pre-apply `prob_func` to get the correct IC for trajectory `i` and
wrap the result in a fresh `TraceProblem` with the default (identity) `prob_func`,
so `_generic_boris!` can safely iterate `1:1` without applying `prob_func` again.
The `TraceProblem` construction is a negligible struct copy relative to the
simulation cost and the serialization overhead inherent in `pmap`.
"""
function _solve_single_boris(
        prob::TraceProblem, i, savestepinterval, dt, nt, nout, isoutofdomain::F, n,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork, F}
    new_prob = prob.prob_func(prob, i, false)
    single_prob = TraceProblem(new_prob.u0, new_prob.tspan, new_prob.p)
    sol_type = _get_sol_type(single_prob, dt, Val(SaveFields), Val(SaveWork))
    local_sols = Vector{sol_type}(undef, 1)
    _dispatch_boris!(
        local_sols, single_prob, 1:1, savestepinterval, dt, nt, nout,
        isoutofdomain, n, save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork)
    )
    return local_sols[1]
end

@inline function _solve(
        ::EnsembleDistributed, prob::TraceProblem, trajectories, dt, savestepinterval,
        isoutofdomain::F, n, save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}, maxiters, batch_size
    ) where {SaveFields, SaveWork, F}
    _, nt, nout = _prepare(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork), maxiters
    )
    return pmap(1:trajectories; batch_size = batch_size) do i
        _solve_single_boris(
            prob, i, savestepinterval, dt, nt, nout, isoutofdomain, n,
            save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
        )
    end
end

function _get_sol_type(prob, dt, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    u0 = prob.u0
    tspan = prob.tspan
    T_t = typeof(tspan[1] + dt)
    t = Vector{T_t}(undef, 0)
    # Force u to be Vector{SVector{6, T}} as used in _boris!
    T = eltype(u0)

    n_vars = 6
    if SaveFields
        n_vars += 6
    end
    if SaveWork
        n_vars += 4
    end

    u = Vector{SVector{n_vars, T}}(undef, 0)
    interp = LinearInterpolation(t, u)
    alg = :boris

    sol = build_solution(prob, alg, t, u; interp = interp)
    return typeof(sol)
end

"""
Prepare for advancing.
"""
function _prepare(
        prob::TraceProblem, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}, maxiters
    ) where {SaveFields, SaveWork}
    if abs(dt) < 10 * eps(typeof(dt))
        throw(ArgumentError("time step dt is too small, violating min_dt = 10 * eps(typeof(dt))"))
    end
    ttotal = prob.tspan[2] - prob.tspan[1]
    nt = round(Int, ttotal / dt) |> abs
    if nt > maxiters
        throw(ArgumentError("number of iterations nt ($nt) exceeds maxiters ($maxiters)"))
    end

    nout = 0
    if save_start
        nout += 1
    end

    if save_everystep
        steps = nt ÷ savestepinterval
        last_is_step = (nt > 0) && (nt % savestepinterval == 0)
        nout += steps
        if !save_end && last_is_step
            nout -= 1
        end
        if save_end && !last_is_step
            nout += 1
        end
    elseif save_end
        nout += 1
    end

    sol_type = _get_sol_type(prob, dt, Val(SaveFields), Val(SaveWork))
    sols = Vector{sol_type}(undef, trajectories)

    return sols, nt, nout
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
Apply Boris method for particles with index in `irange`.
"""
@inline function _boris_loop!(
        traj, tsave, iout, r, v, p, dt, nt, tspan,
        savestepinterval, save_everystep, isoutofdomain::F1, velocity_updater::F2,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {F1, F2, SaveFields, SaveWork}
    it = 1
    while it <= nt
        v_prev = v
        t = (it - 0.5) * dt
        v = velocity_updater(v, r, p, dt, t)

        if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
            iout += 1
            if iout <= length(traj)
                t_current = tspan[1] + (it - 1) * dt
                v_save = velocity_updater(v_prev, r, p, 0.5 * dt, t_current)
                data = vcat(r, v_save)
                traj[iout] = _prepare_saved_data(data, p, t_current, Val(SaveFields), Val(SaveWork))
                tsave[iout] = t_current
            end
        end

        r += v * dt
        isoutofdomain(vcat(r, v), p, it * dt) && break
        it += 1
    end
    return it, iout, r, v
end

@inline @muladd function _generic_boris!(
        sols, prob::TraceProblem, irange, savestepinterval, dt, nt, nout, isoutofdomain::F1,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        velocity_updater::F2, alg_name
    ) where {SaveFields, SaveWork, F1, F2}
    (; tspan, p, u0) = prob
    T = eltype(u0)

    vars_dim = 6
    if SaveFields
        vars_dim += 6
    end
    if SaveWork
        vars_dim += 4
    end

    @inbounds for i in irange
        traj = Vector{SVector{vars_dim, T}}(undef, nout)
        tsave = Vector{typeof(tspan[1] + dt)}(undef, nout)

        # set initial conditions for each trajectory i
        iout = 0
        new_prob = prob.prob_func(prob, i, false)
        u0_i = SVector{6, T}(new_prob.u0)
        r = u0_i[SVector(1, 2, 3)]
        v = u0_i[SVector(4, 5, 6)]

        if save_start
            iout += 1
            traj[iout] = _prepare_saved_data(u0_i, p, tspan[1], Val(SaveFields), Val(SaveWork))
            tsave[iout] = tspan[1]
        end

        # push velocity back in time by 1/2 dt
        v = velocity_updater(v, r, p, -0.5 * dt, tspan[1])

        it, iout, r, v = _boris_loop!(
            traj, tsave, iout, r, v, p, dt, nt, tspan,
            savestepinterval, save_everystep, isoutofdomain, velocity_updater,
            Val(SaveFields), Val(SaveWork)
        )

        final_step = min(it, nt)
        should_save_final = false
        if save_end
            should_save_final = true
        elseif save_everystep && (final_step > 0) && (final_step % savestepinterval == 0)
            should_save_final = true
        end

        if iout < nout && should_save_final
            iout += 1
            t_final = final_step == nt ? tspan[2] : tspan[1] + final_step * dt
            dt_final = t_final - (tspan[1] + (final_step - 0.5) * dt)
            v_final = velocity_updater(v, r, p, dt_final, t_final)

            data = vcat(r, v_final)
            traj[iout] = _prepare_saved_data(
                data, p, t_final, Val(SaveFields), Val(SaveWork)
            )
            tsave[iout] = t_final
        end

        if iout < nout
            resize!(traj, iout)
            resize!(tsave, iout)
        end
        traj_save = traj
        t = tsave

        alg = alg_name
        interp = LinearInterpolation(t, traj_save)
        retcode = ReturnCode.Default
        stats = nothing

        sols[i] = build_solution(
            prob, alg, t, traj_save; interp = interp, retcode = retcode, stats = stats
        )
    end

    return
end

"""
Apply Boris method for particles with index in `irange`.
"""
@inline @muladd function _boris!(
        sols, prob::TraceProblem, irange, savestepinterval, dt, nt, nout, isoutofdomain::F,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork, F}

    _generic_boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
        update_velocity, :boris
    )

    return
end

"""
    update_velocity_multistep(v, r, param, dt, t, n)

Update velocity using the Multistep Boris method, returning the new velocity as an SVector.
Reference: [Zenitani & Kato 2025](https://arxiv.org/abs/2505.02270)
"""
@muladd function update_velocity_multistep(v, r, param, dt, t, n::Int)
    q2m, _, Efunc, Bfunc, _ = param
    E = Efunc(r, t)
    B = Bfunc(r, t)

    # t_n and e_n vectors
    factor = q2m * dt / (2 * n)

    t_n = factor * B # (q/m * dt/(2n)) * B
    e_n = factor * E # (q/m * dt/(2n)) * E

    t_n_mag2 = sum(abs2, t_n)
    t_n_mag = sqrt(t_n_mag2)

    # Calculate coefficients
    # Check for small t_n to avoid division by zero or precision loss
    if t_n_mag < TN_MAG_THRESHOLD
        # Taylor expansion limits as t_n -> 0
        c_n1 = 1 - 2 * n * n * t_n_mag2
        c_n2 = 2 * n - (4 / 3) * n * n * n * t_n_mag2
        c_n3 = 2 * n * n
        c_n6 = (4 / 3) * n * n * n
    else
        alpha_n = atan(t_n_mag)
        n_alpha_n = n * alpha_n
        sin_n_alpha, cos_n_alpha = sincos(n_alpha_n)
        sin_2n_alpha = 2 * sin_n_alpha * cos_n_alpha
        cos_2n_alpha = cos_n_alpha * cos_n_alpha - sin_n_alpha * sin_n_alpha

        c_n1 = cos_2n_alpha
        c_n2 = sin_2n_alpha / t_n_mag
        c_n3 = 2 * sin_n_alpha * sin_n_alpha / t_n_mag2
        c_n6 = (2 * n - c_n2) / t_n_mag2
    end

    c_n4 = c_n2
    c_n5 = c_n3

    v_dot_t = v ⋅ t_n
    e_dot_t = e_n ⋅ t_n

    v_cross_t = v × t_n
    e_cross_t = e_n × t_n

    # Update velocity
    # Equation 39:
    # v_new = c_n1*v + c_n2*(v x t_n) + c_n3*(v . t_n)t_n + c_n4*e_n + c_n5*(e_n x t_n) + c_n6*(e_n . t_n)t_n
    v_new = c_n1 * v +
        c_n2 * v_cross_t +
        c_n3 * v_dot_t * t_n +
        c_n4 * e_n +
        c_n5 * e_cross_t +
        c_n6 * e_dot_t * t_n

    return v_new
end

@inline @muladd function _multistep_boris!(
        sols, prob::TraceProblem, irange, savestepinterval, dt, nt, nout, isoutofdomain::F, n_steps::Int,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork, F}

    velocity_updater = (v, r, p, dt, t) ->
    update_velocity_multistep(v, r, p, dt, t, n_steps)

    _generic_boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
        velocity_updater, :multistep_boris
    )

    return
end

"""
    get_fields(sol::AbstractODESolution)

Return the electric and magnetic fields from the solution `sol`.
"""
function get_fields(sol::AbstractODESolution)
    Efunc, Bfunc = _get_field_funcs(sol.prob)
    T = eltype(sol.u[1])

    E = [SVector{3, T}(Efunc(get_x(u), t)) for (u, t) in zip(sol.u, sol.t)]
    B = [SVector{3, T}(Bfunc(get_x(u), t)) for (u, t) in zip(sol.u, sol.t)]

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
    return [get_work_rates_gc(u, p, t) for (u, t) in zip(sol.u, sol.t)]
end

function _get_work(sol, prob)
    p = prob.p
    return [get_work_rates(u, p, t) for (u, t) in zip(sol.u, sol.t)]
end

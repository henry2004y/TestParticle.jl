# Native particle pusher

"""
    Boris(; safety=0.0)

The standard Boris method for particle pushing in electric and magnetic fields.
When `safety > 0.0`, it uses adaptive time stepping based on local gyroperiod.
"""
struct Boris{T}
    safety::T
end
Boris(; safety = 0.0) = Boris(safety)

"""
    AdaptiveBoris(; safety=0.1)

The adaptive Boris method with adaptive time stepping based on local gyroperiod.
Returns a `Boris` instance with the specified `safety` factor.
"""
struct AdaptiveBoris{T}
    safety::T
end
function AdaptiveBoris(; safety = 0.1)
    @warn "AdaptiveBoris is deprecated. Use Boris(safety=$safety) instead." maxlog=1
    return AdaptiveBoris(safety)
end

"""
    MultistepBoris{N}(; n=1, safety=0.0)

The Multistep/Hyper Boris method of order `N`.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order: 2 (standard), 4, or 6 (Hyper-Boris).
When `safety > 0.0`, it uses adaptive time stepping.
"""
struct MultistepBoris{N, T}
    n::Int
    safety::T
end
@inline function MultistepBoris{N}(; n::Int = 1, safety = 0.0) where {N}
    if N ∉ (2, 4, 6)
        throw(ArgumentError("Multistep Boris order N must be 2, 4, or 6."))
    end
    return MultistepBoris{N, typeof(safety)}(n, safety)
end

const MultistepBoris2 = MultistepBoris{2}
const MultistepBoris4 = MultistepBoris{4}
const MultistepBoris6 = MultistepBoris{6}


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
    solve(prob::TraceProblem, alg::Union{Boris, MultistepBoris}=Boris();
        trajectories::Int=1, dt=nothing,
        savestepinterval::Int=1, isoutside::Function=ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool=true, save_end::Bool=true, save_everystep::Bool=true,
        save_fields::Bool=false, save_work::Bool=false)

Trace particles using the Boris method with specified `prob`.

# keywords

  - `trajectories::Int`: number of trajectories to trace.
  - `dt::AbstractFloat`: time step for fixed-step Boris. If not provided, `alg.safety` must be > 0.
  - `savestepinterval::Int`: saving output interval.
  - `isoutside::Function`: pinpointing impact or checking boundaries.
  - `save_start::Bool=true`: save the initial condition.
  - `save_end::Bool=true`: save the final condition.
  - `save_everystep::Bool=true`: save the state at every `savestepinterval`.
  - `save_fields::Bool=false`: save the electric and magnetic fields.
  - `save_work::Bool=false`: save the work done by the electric field.
  - `batch_size::Int=max(1, trajectories ÷ nworkers())`: the number of trajectories to process per worker in `EnsembleDistributed` and `EnsembleSplitThreads`.
"""
@inline function solve(
        prob::TraceProblem, alg::Union{Boris, MultistepBoris, AdaptiveBoris} = Boris(),
        ensemblealg::EA = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1, dt = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        save_fields::Bool = false, save_work::Bool = false, maxiters::Int = 1_000_000,
        batch_size::Int = (ensemblealg isa EnsembleDistributed || ensemblealg isa EnsembleSplitThreads) ?
            max(1, trajectories ÷ nworkers()) : 1,
        n::Int = 1, N::Int = 2
    ) where {EA <: BasicEnsembleAlgorithm, F}

    # Backward compatibility: promote AdaptiveBoris to Boris
    if alg isa AdaptiveBoris
        alg = Boris(alg.safety)
    end

    # Backward compatibility: promote Boris to MultistepBoris if n > 1 or N != 2
    if alg isa Boris && (n > 1 || N != 2)
        if N ∉ (2, 4, 6)
            throw(ArgumentError("Multistep Boris order N must be 2, 4, or 6."))
        end
        alg = MultistepBoris{N}(n = n, safety = alg.safety)
    end

    if isnothing(dt) && alg.safety <= 0.0
        throw(ArgumentError("Time step dt must be provided for fixed-step Boris solver."))
    end

    if !isnothing(dt)
        if dt < eps(eltype(dt)) * 100
            throw(ArgumentError("Time step dt is too small."))
        end
        ttotal = prob.tspan[2] - prob.tspan[1]
        if abs(ttotal / dt) > maxiters
            throw(ArgumentError("Total steps exceed maxiters. Increase maxiters or dt."))
        end
    end

    return _solve(
        ensemblealg, prob, trajectories, alg, dt, savestepinterval, isoutside,
        save_start, save_end, save_everystep, Val(save_fields), Val(save_work), maxiters,
        batch_size
    )
end

function solve(prob::TraceProblem, ensemblealg::BasicEnsembleAlgorithm; kwargs...)
    return solve(prob, Boris(), ensemblealg; kwargs...)
end

function _dispatch_boris!(
        sols, prob::TraceProblem, irange, alg, dt, savestepinterval, isoutside::F,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        maxiters
    ) where {SaveFields, SaveWork, F}
    return _generic_boris!(
        sols, prob, irange, alg, dt, savestepinterval, isoutside,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
        maxiters
    )
end

@inline function _solve(
        ::EnsembleSerial, prob::TraceProblem, trajectories, alg, dt, savestepinterval,
        isoutside::F, save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}, maxiters, batch_size
    ) where {SaveFields, SaveWork, F}
    sols = _prepare(
        prob, trajectories, alg, dt, save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork)
    )
    irange = 1:trajectories
    _dispatch_boris!(
        sols, prob, irange, alg, dt, savestepinterval, isoutside,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
        maxiters
    )

    return sols
end

@inline function _solve(
        ::EnsembleThreads, prob::TraceProblem, trajectories, alg, dt, savestepinterval,
        isoutside::F, save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}, maxiters, batch_size
    ) where {SaveFields, SaveWork, F}
    sols = _prepare(
        prob, trajectories, alg, dt, save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork)
    )

    nchunks = Threads.nthreads()
    Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        _dispatch_boris!(
            sols, prob, irange, alg, dt, savestepinterval, isoutside,
            save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
            maxiters
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
        prob::TraceProblem, i, alg, dt, savestepinterval, isoutside::F,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        maxiters
    ) where {SaveFields, SaveWork, F}
    new_prob = prob.prob_func(prob, i, false)
    single_prob = TraceProblem(new_prob.u0, new_prob.tspan, new_prob.p)
    sol_type = _get_sol_type(
        single_prob, isnothing(dt) ? zero(eltype(single_prob.tspan)) : dt,
        Val(SaveFields), Val(SaveWork)
    )
    local_sols = Vector{sol_type}(undef, 1)
    _dispatch_boris!(
        local_sols, single_prob, 1:1, alg, dt, savestepinterval,
        isoutside, save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork), maxiters
    )
    return local_sols[1]
end

@inline function _solve(
        ::EnsembleDistributed, prob::TraceProblem, trajectories, alg, dt, savestepinterval,
        isoutside::F, save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}, maxiters, batch_size
    ) where {SaveFields, SaveWork, F}
    return pmap(1:trajectories; batch_size = batch_size) do i
        _solve_single_boris(
            prob, i, alg, dt, savestepinterval, isoutside,
            save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
            maxiters
        )
    end
end

@inline function _solve(
        ::EnsembleSplitThreads, prob::TraceProblem, trajectories, alg, dt, savestepinterval,
        isoutside::F, save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}, maxiters, batch_size
    ) where {SaveFields, SaveWork, F}
    # _solve_single_boris wraps each trajectory in a fresh TraceProblem(u0, tspan, p)
    # with DEFAULT_PROB_FUNC. We get a sample problem from prob_func to ensure the
    # u0 type (uType) matches.
    sample_prob = prob.prob_func(prob, 1, false)
    dummy_prob = TraceProblem(sample_prob.u0, sample_prob.tspan, sample_prob.p)
    sol_type = _get_sol_type(dummy_prob, zero(eltype(dummy_prob.tspan)), Val(SaveFields), Val(SaveWork))
    ichunks = index_chunks(1:trajectories; size = batch_size)
    results = pmap(ichunks) do irange
        local_sols = Vector{sol_type}(undef, length(irange))
        Threads.@threads for k in eachindex(irange)
            i = irange[k]
            local_sols[k] = _solve_single_boris(
                prob, i, alg, dt, savestepinterval, isoutside,
                save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
                maxiters
            )
        end
        local_sols
    end
    return reduce(vcat, results)
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
        prob::TraceProblem, trajectories, alg, dt,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    if !isnothing(dt) && abs(dt) < 10 * eps(typeof(dt))
        throw(ArgumentError("time step dt is too small, violating min_dt = 10 * eps(typeof(dt))"))
    end

    sol_type = _get_sol_type(
        prob, isnothing(dt) ? zero(eltype(prob.tspan)) : dt,
        Val(SaveFields), Val(SaveWork)
    )
    sols = Vector{sol_type}(undef, trajectories)

    return sols
end

@inline function _prepare_saved_data(xv, p, t, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    data = xv
    if SaveFields || SaveWork
        r = get_x(xv)
        T = eltype(xv)
        Bfunc = get_BField(p)
        E_field = SVector{3, T}(get_EField(p)(r, t))

        if SaveWork
            # get_magnetic_properties returns (B, ∇B, κ, b̂, Bmag)
            magnetic_props = get_magnetic_properties(r, t, Bfunc)
            if SaveFields
                B_vec = SVector{3, T}(magnetic_props[1])
                data = vcat(data, E_field, B_vec)
            end
            work = get_work_rates(xv, p, t, magnetic_props, E_field)
            data = vcat(data, work)
        elseif SaveFields
            B_vec = SVector{3, T}(Bfunc(r, t))
            data = vcat(data, E_field, B_vec)
        end
    end

    return data
end

@inline @muladd function _boris_loop!(
        traj, tsave, r, v, p, dt, tspan,
        savestepinterval, save_everystep, isoutside::F1, velocity_updater::F2,
        ::Val{SaveFields}, ::Val{SaveWork}, alg, maxiters
    ) where {F1, F2, SaveFields, SaveWork}
    t = tspan[1]
    ttotal = tspan[2] - tspan[1]
    it = 1
    retcode = ReturnCode.Success

    C = zero(eltype(v))
    if alg.safety > 0.0
        q2m, _, _, Bfunc, _ = p
        C = (2π * alg.safety * sign(ttotal)) / abs(q2m)
        Bmag = norm(Bfunc(r, t))
        dt = C / Bmag
    end

    # push velocity back in time by 1/2 dt
    v = velocity_updater(v, r, -0.5 * dt, t, p)

    while abs(t - tspan[1]) < abs(ttotal)
        if abs(t + dt - tspan[1]) >= abs(ttotal) - 100 * eps(eltype(v)(abs(ttotal)))
            dt_step = tspan[2] - t
            if dt_step != dt
                # Resync v from `t - 0.5*dt` to `t - 0.5*dt_step`
                v = velocity_updater(v, r, 0.5 * dt, t, p)
                v = velocity_updater(v, r, -0.5 * dt_step, t, p)
                dt = dt_step
            end
        end

        if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
            # Advance to t to get v_n
            v_save = velocity_updater(v, r, 0.5 * dt, t, p)
            data = vcat(r, v_save)
            push!(traj, _prepare_saved_data(data, p, t, Val(SaveFields), Val(SaveWork)))
            push!(tsave, t)
        end

        v_new = velocity_updater(v, r, dt, t + 0.5 * dt, p)
        r_next = r + v_new * dt
        t_next = t + dt

        # NaN check
        if any(isnan, r_next) || any(isnan, v_new)
            retcode = ReturnCode.Unstable
            break
        end

        if isoutside(vcat(r_next, v_new), p, t_next)
            retcode = ReturnCode.Terminated
            break
        end

        r, v, t = r_next, v_new, t_next

        if alg.safety > 0.0 && abs(t - tspan[1]) < abs(ttotal)
            q2m, _, _, Bfunc, _ = p
            Bmag = norm(Bfunc(r, t))
            dt_new = C / Bmag
            # Resync v_{n+1/2}(dt) to v_{n+1/2}(dt_new)
            v = velocity_updater(v, r, 0.5 * dt, t, p)
            v = velocity_updater(v, r, -0.5 * dt_new, t, p)
            dt = dt_new
        end
        if it > maxiters
            retcode = ReturnCode.MaxIters
            break
        end
        it += 1
    end

    return r, v, t, dt, retcode
end

@inline @muladd function _generic_boris!(
        sols, prob::TraceProblem, irange, alg, dt, savestepinterval, isoutside::F1,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        maxiters
    ) where {SaveFields, SaveWork, F1}
    (; tspan, p, u0) = prob
    T = eltype(u0)

    vars_dim = 6
    if SaveFields
        vars_dim += 6
    end
    if SaveWork
        vars_dim += 4
    end

    algorithm_name = alg isa Boris ? :boris : :multistep_boris

    velocity_updater = if alg isa Boris
        update_velocity
    else
        (v, r, dt, t, p) -> update_velocity_multistep(v, r, dt, t, alg.n, _get_val_N(alg), p)
    end

    # Calculate exact nout for fixed-step to avoid memory regression
    ttotal = tspan[2] - tspan[1]
    nout_fixed = if !isnothing(dt)
        nt = round(Int, abs(ttotal / dt))
        nsteps = 0
        if save_start
            nsteps += 1
        end
        if save_everystep
            steps = nt ÷ savestepinterval
            last_is_step = (nt > 0) && (nt % savestepinterval == 0)
            nsteps += steps
            if !save_end && last_is_step
                nsteps -= 1
            end
            if save_end && !last_is_step
                nsteps += 1
            end
        elseif save_end
            nsteps += 1
        end
        nsteps
    else
        1000 # Default capacity for adaptive Boris
    end

    @inbounds for i in irange
        traj = Vector{SVector{vars_dim, T}}(undef, 0)
        tsave = Vector{typeof(tspan[1] + (isnothing(dt) ? 0.0 : dt))}(undef, 0)
        sizehint!(traj, nout_fixed)
        sizehint!(tsave, nout_fixed)

        # set initial conditions for each trajectory i
        new_prob = prob.prob_func(prob, i, false)
        u0_i = SVector{6, T}(new_prob.u0)
        r = u0_i[SVector(1, 2, 3)]
        v = u0_i[SVector(4, 5, 6)]

        # If dt is not provided, it must be adaptive (safety > 0)
        _dt = isnothing(dt) ? zero(T) : T(dt)

        if save_start
            push!(traj, _prepare_saved_data(u0_i, p, tspan[1], Val(SaveFields), Val(SaveWork)))
            push!(tsave, tspan[1])
        end

        r, v, t, _dt_final, retcode = _boris_loop!(
            traj, tsave, r, v, p, _dt, tspan,
            savestepinterval, save_everystep, isoutside, velocity_updater,
            Val(SaveFields), Val(SaveWork), alg, maxiters
        )

        should_save_final = save_end
        if should_save_final && (isempty(tsave) || tsave[end] != t)
            # v is at t - 0.5*_dt_final. To get v at t, advance by 0.5*_dt_final
            v_save = velocity_updater(v, r, 0.5 * _dt_final, t, p)
            data = vcat(r, v_save)
            push!(traj, _prepare_saved_data(data, p, t, Val(SaveFields), Val(SaveWork)))
            push!(tsave, t)
        end

        alg_sol = algorithm_name
        interp = LinearInterpolation(tsave, traj)
        stats = nothing

        sols[i] = build_solution(prob, alg_sol, tsave, traj; interp, retcode, stats)
    end

    return
end

_get_val_N(::MultistepBoris{N}) where {N} = Val{N}()

"""
    update_velocity_multistep(v, r, dt, t, n, N, param)
    update_velocity_multistep(v, r, dt, t, n, ::Val{N}, param)

Update velocity using the Multistep/Hyper Boris method, returning the new velocity as an SVector.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order. When N=2, it corresponds to the Multicycle solver. When N=4 or N=6, it is the Hyper Boris solver.
Reference: [Zenitani & Kato 2025](https://arxiv.org/abs/2505.02270)
"""
@inline function update_velocity_multistep(v, r, dt, t, n::Int, N::Int, param)
    return update_velocity_multistep(v, r, dt, t, n, Val(N), param)
end

@muladd function update_velocity_multistep(v, r, dt, t, n::Int, ::Val{N}, param) where {N}
    q2m, _, Efunc, Bfunc, _ = param
    E = Efunc(r, t)
    B = Bfunc(r, t)

    # t_n and e_n vectors
    factor = q2m * dt / (2 * n)

    t_n = factor * B # (q/m * dt/(2n)) * B
    e_n = factor * E # (q/m * dt/(2n)) * E

    # Hyper Boris N-th order gyrophase correction
    if N != 2
        t_mag2 = sum(abs2, t_n)
        if N == 4
            f_N = 1 + t_mag2 / 3
            e_corr_factor = -1 / 3
        else # N == 6
            f_N = 1 + t_mag2 / 3 + 2 * t_mag2 * t_mag2 / 15
            e_corr_factor = -1 / 3 - 2 * t_mag2 / 15
        end

        e_dot_t = e_n ⋅ t_n
        e_n = f_N * e_n + (e_corr_factor * e_dot_t) * t_n
        t_n = f_N * t_n
    end

    t_n_mag2 = sum(abs2, t_n)
    t_n_mag = sqrt(t_n_mag2)

    # Calculate coefficients
    # Check for small t_n to avoid division by zero or precision loss
    if t_n_mag < TN_MAG_THRESHOLD
        # Taylor expansion limits as t_n -> 0
        c_n1 = 1 - 2 * n * n * t_n_mag2

        n_term1 = 2 * n
        n_term3 = 4 * n * n * n

        c_n2 = n_term1 - (n_term1 + n_term3) / 3 * t_n_mag2
        c_n3 = 2 * n * n - (4 * n * n + 2 * n * n * n * n) / 3 * t_n_mag2
        c_n6 = (n_term1 + n_term3) / 3
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

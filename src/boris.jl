# Native particle pusher

struct Boris{Adaptive} <: AbstractBoris
    safety::Float64
end

Boris() = Boris{false}(0.0)

"""
    MultistepBoris{N}(; n=1)

The Multistep/Hyper Boris method of order `N`.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order:
2 (standard), 4, or 6 (Hyper-Boris).
"""
struct MultistepBoris{N} <: AbstractBoris
    n::Int
end

@inline function MultistepBoris{N}(; n::Int = 1) where {N}
    if N ∉ (2, 4, 6)
        throw(
            ArgumentError(
                "Multistep Boris order N must be 2, 4, or 6."
            )
        )
    end
    return MultistepBoris{N}(n)
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
    solve(prob::TraceProblem,
        alg::Union{Boris{false}, MultistepBoris},
        ensemblealg=EnsembleSerial(); dt, kwargs...)

Trace particles using the fixed-step Boris or Multistep/Hyper Boris method.

# Keywords
  - `trajectories::Int=1`: number of trajectories.
  - `dt::AbstractFloat`: time step.
  - `savestepinterval::Int=1`: saving output interval.
  - `isoutside`: boundary check function.
  - `save_start::Bool=true`: save initial condition.
  - `save_end::Bool=true`: save final condition.
  - `save_everystep::Bool=true`: save at intervals.
  - `save_fields::Bool=false`: save E and B fields.
  - `save_work::Bool=false`: save work rates.
  - `maxiters::Int=1_000_000`: maximum iterations.
"""
@inline function solve(
        prob::TraceProblem,
        alg::Union{Boris{false}, MultistepBoris},
        ensemblealg::EA = EnsembleSerial();
        trajectories::Int = 1,
        savestepinterval::Int = 1,
        dt::AbstractFloat,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true,
        save_end::Bool = true,
        save_everystep::Bool = true,
        save_fields::Bool = false,
        save_work::Bool = false,
        maxiters::Int = 1_000_000,
        batch_size::Int = _default_batch_size(
            ensemblealg, trajectories
        ),
    ) where {EA <: BasicEnsembleAlgorithm, F}
    return _solve(
        ensemblealg, prob, alg, trajectories, dt,
        savestepinterval, isoutside,
        save_start, save_end, save_everystep,
        Val(save_fields), Val(save_work),
        maxiters, batch_size
    )
end

function _default_batch_size(ensemblealg, trajectories)
    if ensemblealg isa EnsembleDistributed ||
            ensemblealg isa EnsembleSplitThreads
        return max(1, trajectories ÷ nworkers())
    end
    return 1
end

function _dispatch_boris!(
        sols, prob::TraceProblem, irange,
        savestepinterval, dt, nt, nout, isoutside::F,
        ::Boris{false}, save_start, save_end,
        save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork, F}
    return _boris!(
        sols, prob, irange, savestepinterval, dt, nt,
        nout, isoutside,
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork)
    )
end

function _dispatch_boris!(
        sols, prob::TraceProblem, irange,
        savestepinterval, dt, nt, nout, isoutside::F,
        alg::MultistepBoris{N}, save_start, save_end,
        save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {N, SaveFields, SaveWork, F}
    return _multistep_boris!(
        sols, prob, irange, savestepinterval, dt, nt,
        nout, isoutside, alg.n, Val(N),
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork)
    )
end

@inline function _solve(
        ::EnsembleSerial, prob::TraceProblem,
        alg::AbstractBoris, trajectories, dt,
        savestepinterval, isoutside::F,
        save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork},
        maxiters, batch_size
    ) where {SaveFields, SaveWork, F}
    sols, nt,
        nout = _prepare(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork), maxiters
    )
    irange = 1:trajectories
    _dispatch_boris!(
        sols, prob, irange, savestepinterval, dt, nt,
        nout, isoutside, alg,
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork)
    )

    return sols
end

@inline function _solve(
        ::EnsembleThreads, prob::TraceProblem,
        alg::AbstractBoris, trajectories, dt,
        savestepinterval, isoutside::F,
        save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork},
        maxiters, batch_size
    ) where {SaveFields, SaveWork, F}
    sols, nt,
        nout = _prepare(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork), maxiters
    )

    nchunks = Threads.nthreads()
    Threads.@threads for irange in index_chunks(
            1:trajectories; n = nchunks
        )
        _dispatch_boris!(
            sols, prob, irange, savestepinterval,
            dt, nt, nout, isoutside, alg,
            save_start, save_end, save_everystep,
            Val(SaveFields), Val(SaveWork)
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
        prob::TraceProblem, i, savestepinterval,
        dt, nt, nout, isoutside::F,
        alg::AbstractBoris,
        save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork, F}
    new_prob = prob.prob_func(prob, i, false)
    single_prob = TraceProblem(
        new_prob.u0, new_prob.tspan, new_prob.p
    )
    sol_type = _get_sol_type(
        single_prob, dt, Val(SaveFields), Val(SaveWork)
    )
    local_sols = Vector{sol_type}(undef, 1)
    _dispatch_boris!(
        local_sols, single_prob, 1:1,
        savestepinterval, dt, nt, nout,
        isoutside, alg,
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork)
    )
    return local_sols[1]
end

@inline function _solve(
        ::EnsembleDistributed, prob::TraceProblem,
        alg::AbstractBoris, trajectories, dt,
        savestepinterval, isoutside::F,
        save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork},
        maxiters, batch_size
    ) where {SaveFields, SaveWork, F}
    _, nt, nout = _prepare(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork), maxiters
    )
    return pmap(
        1:trajectories; batch_size = batch_size
    ) do i
        _solve_single_boris(
            prob, i, savestepinterval, dt, nt, nout,
            isoutside, alg,
            save_start, save_end, save_everystep,
            Val(SaveFields), Val(SaveWork)
        )
    end
end

@inline function _solve(
        ::EnsembleSplitThreads, prob::TraceProblem,
        alg::AbstractBoris, trajectories, dt,
        savestepinterval, isoutside::F,
        save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork},
        maxiters, batch_size
    ) where {SaveFields, SaveWork, F}
    _, nt, nout = _prepare(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork), maxiters
    )
    sample_prob = prob.prob_func(prob, 1, false)
    dummy_prob = TraceProblem(
        sample_prob.u0, sample_prob.tspan,
        sample_prob.p
    )
    sol_type = _get_sol_type(
        dummy_prob, dt, Val(SaveFields), Val(SaveWork)
    )
    ichunks = index_chunks(
        1:trajectories; size = batch_size
    )
    results = pmap(ichunks) do irange
        local_sols = Vector{sol_type}(
            undef, length(irange)
        )
        Threads.@threads for k in eachindex(irange)
            i = irange[k]
            local_sols[k] = _solve_single_boris(
                prob, i, savestepinterval,
                dt, nt, nout, isoutside, alg,
                save_start, save_end,
                save_everystep,
                Val(SaveFields), Val(SaveWork)
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
@inline @muladd function _boris_loop!(
        traj, tsave, iout, r, v, p, dt, nt, tspan,
        savestepinterval, save_everystep, isoutside::F1, velocity_updater::F2,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {F1, F2, SaveFields, SaveWork}
    it = 1
    t = tspan[1] - 0.5 * dt
    while it <= nt
        v_prev = v
        t += dt
        v = velocity_updater(v, r, dt, t, p)

        r_next = r + v * dt
        t_next = t + 0.5 * dt
        if isoutside(vcat(r_next, v), p, t_next)
            return it - 1, iout, r, v_prev
        end

        if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
            iout += 1
            if iout <= length(traj)
                t_current = t - 0.5 * dt
                v_save = velocity_updater(v_prev, r, 0.5 * dt, t_current, p)
                data = vcat(r, v_save)
                traj[iout] = _prepare_saved_data(data, p, t_current, Val(SaveFields), Val(SaveWork))
                tsave[iout] = t_current
            end
        end

        r = r_next
        it += 1
    end
    return it - 1, iout, r, v
end

@inline @muladd function _generic_boris!(
        sols, prob::TraceProblem, irange, savestepinterval, dt, nt, nout, isoutside::F1,
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
        v = velocity_updater(v, r, -0.5 * dt, tspan[1], p)

        it, iout, r, v = _boris_loop!(
            traj, tsave, iout, r, v, p, dt, nt, tspan,
            savestepinterval, save_everystep, isoutside, velocity_updater,
            Val(SaveFields), Val(SaveWork)
        )

        final_step = it
        should_save_final = false
        if save_end
            should_save_final = true
        elseif save_everystep && (final_step > 0) && (final_step % savestepinterval == 0)
            should_save_final = true
        end

        if iout < nout && should_save_final
            t_final = final_step == nt ? tspan[2] : tspan[1] + final_step * dt
            if iout == 0 || tsave[iout] < t_final
                iout += 1
                dt_final = t_final - (tspan[1] + (final_step - 0.5) * dt)
                v_final = velocity_updater(v, r, dt_final, t_final, p)

                data = vcat(r, v_final)
                traj[iout] = _prepare_saved_data(
                    data, p, t_final, Val(SaveFields), Val(SaveWork)
                )
                tsave[iout] = t_final
            end
        end

        if iout < nout
            resize!(traj, iout)
            resize!(tsave, iout)
            retcode = ReturnCode.Terminated
        else
            retcode = ReturnCode.Success
        end
        traj_save = traj
        t = tsave

        alg = alg_name
        interp = LinearInterpolation(t, traj_save)
        stats = nothing

        sols[i] = build_solution(prob, alg, t, traj_save; interp, retcode, stats)
    end

    return
end

"""
Apply Boris method for particles with index in `irange`.
"""
@inline @muladd function _boris!(
        sols, prob::TraceProblem, irange, savestepinterval, dt, nt, nout, isoutside::F,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork, F}

    _generic_boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutside,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
        update_velocity, :boris
    )

    return
end

"""
    update_velocity_multistep(v, r, dt, t, n, ::Val{N}, param)

Update velocity using the Multistep/Hyper Boris method, returning the new velocity as an SVector.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order. When N=2, it corresponds to the Multicycle solver. When N=4 or N=6, it is the Hyper Boris solver.
Reference: [Zenitani & Kato 2025](https://arxiv.org/abs/2505.02270)
"""
@muladd function update_velocity_multistep(
        v, r, dt, t, n::Int, ::Val{N}, param
    ) where {N}
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

@inline @muladd function _multistep_boris!(
        sols, prob::TraceProblem, irange, savestepinterval, dt, nt, nout, isoutside::F,
        n_steps::Int, ::Val{N_order}, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {N_order, SaveFields, SaveWork, F}

    velocity_updater = (v, r, dt, t, p) ->
    update_velocity_multistep(
        v, r, dt, t, n_steps, Val(N_order), p
    )

    _generic_boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutside,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
        velocity_updater, :multistep_boris
    )

    return
end

# --- Adaptive Boris solve entry point and ensemble methods ---

"""
    solve(prob::TraceProblem, alg::Boris{true},
        ensemblealg=EnsembleSerial(); kwargs...)

Trace particles using the adaptive Boris method.
The time step is determined by `dt = safety * 2π / |qB/m|`.

# Keywords
  - `trajectories::Int=1`: number of trajectories.
  - `savestepinterval::Int=1`: saving output interval.
  - `isoutside`: boundary check function.
  - `save_start::Bool=true`: save initial condition.
  - `save_end::Bool=true`: save final condition.
  - `save_everystep::Bool=true`: save at intervals.
  - `save_fields::Bool=false`: save E and B fields.
  - `save_work::Bool=false`: save work rates.
  - `batch_size::Int`: batch size for distributed.
"""
function solve(
        prob::TraceProblem, alg::Boris{true},
        ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1,
        savestepinterval::Int = 1,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true,
        save_end::Bool = true,
        save_everystep::Bool = true,
        save_fields::Bool = false,
        save_work::Bool = false,
        batch_size::Int = _default_batch_size(ensemblealg, trajectories),
    ) where {F}
    return _solve_adaptive(
        ensemblealg, prob, trajectories, alg, savestepinterval, isoutside, save_start,
        save_end, save_everystep, Val(save_fields), Val(save_work), batch_size
    )
end

function _solve_adaptive(
        ::EnsembleSerial, prob, trajectories, alg::Boris{true}, savestepinterval, isoutside,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}, batch_size
    ) where {SaveFields, SaveWork}
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)), Val(SaveFields), Val(SaveWork))
    sols = Vector{sol_type}(undef, trajectories)
    irange = 1:trajectories

    _adaptive_boris!(
        sols, prob, irange, alg, savestepinterval, isoutside, save_start, save_end,
        save_everystep, Val(SaveFields), Val(SaveWork)
    )

    return sols
end

function _solve_adaptive(
        ::EnsembleThreads, prob, trajectories, alg::Boris{true}, savestepinterval,
        isoutside, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        batch_size
    ) where {SaveFields, SaveWork}
    sol_type = _get_sol_type(
        prob, zero(eltype(prob.tspan)),
        Val(SaveFields), Val(SaveWork)
    )
    sols = Vector{sol_type}(undef, trajectories)

    nchunks = Threads.nthreads()
    Threads.@threads for irange in index_chunks(
            1:trajectories; n = nchunks
        )
        _adaptive_boris!(
            sols, prob, irange, alg, savestepinterval, isoutside, save_start, save_end,
            save_everystep, Val(SaveFields), Val(SaveWork)
        )
    end

    return sols
end

"See `_solve_single_boris` for rationale."
function _solve_single_adaptive_boris(
        prob, i, alg::Boris{true}, savestepinterval, isoutside, save_start, save_end,
        save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    new_prob = prob.prob_func(prob, i, false)
    single_prob = TraceProblem(
        new_prob.u0, new_prob.tspan, new_prob.p
    )
    sol_type = _get_sol_type(
        single_prob,
        zero(eltype(single_prob.tspan)),
        Val(SaveFields), Val(SaveWork)
    )
    local_sols = Vector{sol_type}(undef, 1)
    _adaptive_boris!(
        local_sols, single_prob, 1:1, alg,
        savestepinterval, isoutside,
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork)
    )
    return local_sols[1]
end

function _solve_adaptive(
        ::EnsembleDistributed, prob, trajectories, alg::Boris{true}, savestepinterval,
        isoutside, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        batch_size
    ) where {SaveFields, SaveWork}
    return pmap(
        1:trajectories; batch_size = batch_size
    ) do i
        _solve_single_adaptive_boris(
            prob, i, alg, savestepinterval,
            isoutside, save_start, save_end,
            save_everystep,
            Val(SaveFields), Val(SaveWork)
        )
    end
end

function _solve_adaptive(
        ::EnsembleSplitThreads, prob, trajectories, alg::Boris{true}, savestepinterval,
        isoutside, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        batch_size
    ) where {SaveFields, SaveWork}
    sample_prob = prob.prob_func(prob, 1, false)
    dummy_prob = TraceProblem(
        sample_prob.u0, sample_prob.tspan,
        sample_prob.p
    )
    sol_type = _get_sol_type(
        dummy_prob,
        zero(eltype(dummy_prob.tspan)),
        Val(SaveFields), Val(SaveWork)
    )
    ichunks = index_chunks(
        1:trajectories; size = batch_size
    )
    results = pmap(ichunks) do irange
        local_sols = Vector{sol_type}(
            undef, length(irange)
        )
        Threads.@threads for k in eachindex(irange)
            i = irange[k]
            local_sols[k] =
                _solve_single_adaptive_boris(
                prob, i, alg, savestepinterval, isoutside, save_start, save_end,
                save_everystep, Val(SaveFields), Val(SaveWork)
            )
        end
        local_sols
    end
    return reduce(vcat, results)
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

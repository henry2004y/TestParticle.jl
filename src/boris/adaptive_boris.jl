# Adaptive Boris solver

"""
    solve(prob::TraceProblem, alg::Union{Boris{true}, MultistepBoris{N, true} where N},
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
@inline function solve(
        prob::TraceProblem, alg::Union{Boris{true}, MultistepBoris{N, true} where {N}},
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

@inline function _solve_adaptive(
        ::EnsembleSerial, prob, trajectories, alg::AbstractBoris, savestepinterval, isoutside,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}, batch_size
    ) where {SaveFields, SaveWork}
    # Use array comprehension for EnsembleSerial to avoid _get_sol_type allocations
    elapsed_time = @elapsed sols = [
        _adaptive_boris_single(
                prob, i, savestepinterval, isoutside, save_start, save_end,
                save_everystep, Val(SaveFields), Val(SaveWork), alg
            ) for i in 1:trajectories
    ]
    return EnsembleSolution(sols, elapsed_time, true)
end

function _solve_adaptive(
        ::EnsembleThreads, prob, trajectories, alg::AbstractBoris, savestepinterval,
        isoutside, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        batch_size
    ) where {SaveFields, SaveWork}
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)), Val(SaveFields), Val(SaveWork))
    sols = Vector{sol_type}(undef, trajectories)

    nchunks = Threads.nthreads()
    elapsed_time = @elapsed Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        _adaptive_boris!(
            sols, prob, irange, savestepinterval, isoutside, save_start, save_end,
            save_everystep, Val(SaveFields), Val(SaveWork), alg
        )
    end

    return EnsembleSolution(sols, elapsed_time, true)
end

"See `_solve_single_boris` for rationale."
function _solve_single_adaptive_boris(
        prob, i, savestepinterval, isoutside, save_start, save_end,
        save_everystep, ::Val{SaveFields}, ::Val{SaveWork}, alg::AbstractBoris
    ) where {SaveFields, SaveWork}
    new_prob = prob.prob_func(prob, (i = i, repeat = false))
    single_prob = TraceProblem(
        new_prob.u0, new_prob.tspan, new_prob.p
    )
    sol_type = _get_sol_type(
        single_prob, zero(eltype(single_prob.tspan)),
        Val(SaveFields), Val(SaveWork)
    )
    local_sols = Vector{sol_type}(undef, 1)
    _adaptive_boris!(
        local_sols, single_prob, 1:1, savestepinterval, isoutside,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork), alg
    )
    return local_sols[1]
end

function _solve_adaptive(
        ::EnsembleDistributed, prob, trajectories, alg::AbstractBoris, savestepinterval,
        isoutside, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        batch_size
    ) where {SaveFields, SaveWork}
    start_time = time()
    sols = pmap(1:trajectories; batch_size = batch_size) do i
        _solve_single_adaptive_boris(
            prob, i, savestepinterval,
            isoutside, save_start, save_end,
            save_everystep,
            Val(SaveFields), Val(SaveWork), alg
        )
    end
    end_time = time()
    return EnsembleSolution(sols, end_time - start_time, true)
end

function _solve_adaptive(
        ::EnsembleSplitThreads, prob, trajectories, alg::AbstractBoris, savestepinterval,
        isoutside, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        batch_size
    ) where {SaveFields, SaveWork}
    sample_prob = prob.prob_func(prob, (i = 1, repeat = false))
    dummy_prob = TraceProblem(sample_prob.u0, sample_prob.tspan, sample_prob.p)
    sol_type = _get_sol_type(
        dummy_prob, zero(eltype(dummy_prob.tspan)), Val(SaveFields), Val(SaveWork)
    )
    ichunks = index_chunks(1:trajectories; size = batch_size)
    start_time = time()
    results = pmap(ichunks) do irange
        local_sols = Vector{sol_type}(undef, length(irange))
        Threads.@threads for k in eachindex(irange)
            i = irange[k]
            local_sols[k] =
                _solve_single_adaptive_boris(
                prob, i, savestepinterval, isoutside, save_start, save_end,
                save_everystep, Val(SaveFields), Val(SaveWork), alg
            )
        end
        local_sols
    end
    end_time = time()
    sols = reduce(vcat, results)
    return EnsembleSolution(sols, end_time - start_time, true)
end

@inline get_velocity_updater(::Boris{true}) = update_velocity
@inline get_velocity_updater(alg::MultistepBoris{N, true}) where {N} =
    MultistepUpdater{N}(alg.n)

@inline function _adaptive_boris_single(
        prob, i, savestepinterval, isoutside, save_start, save_end,
        save_everystep, ::Val{SaveFields}, ::Val{SaveWork}, alg::AbstractBoris
    ) where {SaveFields, SaveWork}
    (; tspan, p, u0) = prob
    q2m, _, _, Bfunc, _ = p
    T = eltype(u0)
    vars_dim = 6
    if SaveFields
        vars_dim += 6
    end
    if SaveWork
        vars_dim += 4
    end

    velocity_updater = get_velocity_updater(alg)
    alg_name = alg isa Boris ? :adaptive_boris : :adaptive_multistep_boris

    # dt = safety * 2π / (abs(q2m) * Bmag)
    C = (2π * alg.safety * sign(tspan[2] - tspan[1])) / abs(q2m)

    initial_capacity = 1000
    traj = Vector{SVector{vars_dim, T}}(undef, 0)
    tsave = Vector{typeof(tspan[1])}(undef, 0)
    sizehint!(traj, initial_capacity)
    sizehint!(tsave, initial_capacity)

    new_prob = prob.prob_func(prob, (i = i, repeat = false))
    u0_i = SVector{6, T}(new_prob.u0)
    r = u0_i[SVector(1, 2, 3)]
    v = u0_i[SVector(4, 5, 6)]
    t = tspan[1]
    ttotal = tspan[2] - tspan[1]

    if save_start
        data = _prepare_saved_data(
            u0_i, p, t,
            Val(SaveFields), Val(SaveWork)
        )
        push!(traj, data)
        push!(tsave, t)
    end

    Bmag = norm(Bfunc(r, t))
    dt = C / Bmag

    # Backstep velocity: v(0) -> v(-1/2)
    v = velocity_updater(v, r, -0.5 * dt, t, p)
    it = 1
    should_save_final = save_end
    retcode = ReturnCode.Success

    @fastmath while abs(t - tspan[1]) < abs(ttotal)
        if abs(t + dt - tspan[1]) > abs(ttotal)
            dt_step = tspan[2] - t
            v = velocity_updater(v, r, 0.5 * dt, t, p)
            v = velocity_updater(v, r, -0.5 * dt_step, t, p)
            dt = dt_step
        end

        if save_everystep &&
                (it - 1) > 0 &&
                (it - 1) % savestepinterval == 0
            v_save = velocity_updater(
                v, r, 0.5 * dt, t, p
            )
            xv_s = vcat(r, v_save)
            data = _prepare_saved_data(
                xv_s, p, t,
                Val(SaveFields), Val(SaveWork)
            )
            push!(traj, data)
            push!(tsave, t)
        end

        t_mid = t + 0.5 * dt
        v_new = velocity_updater(v, r, dt, t_mid, p)

        r_next = r + v_new * dt
        t_next = t + dt

        xv_new = vcat(r_next, v_new)
        if isoutside(xv_new, p, t_next)
            should_save_final = true
            retcode = ReturnCode.Terminated
            break
        end

        r = r_next
        t = t_next
        v = v_new

        Bmag = norm(Bfunc(r, t))
        dt_new = C / Bmag

        # Resync v_{n+1/2}(dt) to v_{n+1/2}(dt_new)
        v = velocity_updater(v, r, 0.5 * dt, t, p)
        v = velocity_updater(v, r, -0.5 * dt_new, t, p)
        dt = dt_new
        it += 1
        if save_everystep && (it - 1) % savestepinterval == 0
            should_save_final = true
        end
    end

    if should_save_final && (isempty(tsave) || tsave[end] != t)
        v_final = velocity_updater(v, r, 0.5 * dt, t, p)
        xv_f = vcat(r, v_final)
        data = _prepare_saved_data(
            xv_f, p, t,
            Val(SaveFields), Val(SaveWork)
        )
        push!(traj, data)
        push!(tsave, t)
    end

    sol_alg = alg_name
    interp = LinearInterpolation(tsave, traj)
    stats = nothing

    return build_solution(prob, sol_alg, tsave, traj; interp, retcode, stats)
end

@muladd function _adaptive_boris!(
        sols, prob, irange, savestepinterval, isoutside,
        save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}, alg::AbstractBoris
    ) where {SaveFields, SaveWork}

    @inbounds for i in irange
        sols[i] = _adaptive_boris_single(
            prob, i, savestepinterval, isoutside, save_start, save_end, save_everystep,
            Val(SaveFields), Val(SaveWork), alg
        )
    end

    return
end

# Adaptive Boris method

struct AdaptiveBoris{T}
    safety::T
end

"""
    AdaptiveBoris(; safety=0.1)

Adaptive Boris method with adaptive time stepping based on local gyroperiod.
The time step is determined by `dt = safety * T_gyro = safety * 2π / |qB/m|`.
"""
function AdaptiveBoris(; safety = 0.1)
    T = typeof(safety)
    return AdaptiveBoris{T}(T(safety))
end
"""
    solve(prob::TraceProblem, alg::AdaptiveBoris,
        ensemblealg::BasicEnsembleAlgorithm=EnsembleSerial();
        trajectories::Int=1, savestepinterval::Int=1,
        isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool=true, save_end::Bool=true, save_everystep::Bool=true,
        save_fields::Bool=false, save_work::Bool=false,
        batch_size::Int = max(1, trajectories ÷ Threads.nthreads()))

Trace particles using the Adaptive Boris method with specified `prob` and `alg`.

# keywords

  - `trajectories::Int`: number of trajectories to trace.
  - `savestepinterval::Int`: saving output interval.
  - `isoutofdomain::Function`: a function with input of position and velocity vector `xv` that determines whether to stop tracing.
  - `save_start::Bool`: save the initial condition. Default is `true`.
  - `save_end::Bool`: save the final condition. Default is `true`.
  - `save_everystep::Bool`: save the state at every `savestepinterval`. Default is `true`.
  - `save_fields::Bool`: save the electric and magnetic fields. Default is `false`.
  - `save_work::Bool`: save the work done by the electric field. Default is `false`.
  - `batch_size::Int`: the number of trajectories to process per worker in `EnsembleDistributed` and `EnsembleSplitThreads`. Default is `max(1, trajectories ÷ nworkers())` for distributed and 1 for others.


"""
function solve(
        prob::TraceProblem, alg::AdaptiveBoris,
        ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1,
        isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        save_fields::Bool = false, save_work::Bool = false,
        batch_size::Int = (ensemblealg isa EnsembleDistributed || ensemblealg isa EnsembleSplitThreads) ?
            max(1, trajectories ÷ nworkers()) : 1
    )

    return if save_fields
        if save_work
            return _solve(
                ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, Val(true), Val(true), batch_size
            )
        else
            return _solve(
                ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, Val(true), Val(false), batch_size
            )
        end
    else
        if save_work
            return _solve(
                ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, Val(false), Val(true), batch_size
            )
        else
            return _solve(
                ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, Val(false), Val(false), batch_size
            )
        end
    end
end

function _solve(
        ::EnsembleSerial, prob, trajectories, alg::AdaptiveBoris, savestepinterval,
        isoutofdomain, save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}, batch_size
    ) where {SaveFields, SaveWork}
    # We cannot precalculate nt for adaptive steps
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)), Val(SaveFields), Val(SaveWork))
    sols = Vector{sol_type}(undef, trajectories)
    irange = 1:trajectories

    _adaptive_boris!(
        sols, prob, irange, alg, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
    )

    return sols
end

function _solve(
        ::EnsembleThreads, prob, trajectories, alg::AdaptiveBoris, savestepinterval,
        isoutofdomain, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        batch_size
    ) where {SaveFields, SaveWork}
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)), Val(SaveFields), Val(SaveWork))
    sols = Vector{sol_type}(undef, trajectories)

    nchunks = Threads.nthreads()
    Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        _adaptive_boris!(
            sols, prob, irange, alg, savestepinterval, isoutofdomain,
            save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
        )
    end

    return sols
end

"See `_solve_single_boris` for the rationale behind the `single_prob` construction."
function _solve_single_adaptive_boris(
        prob, i, alg::AdaptiveBoris, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    new_prob = prob.prob_func(prob, i, false)
    single_prob = TraceProblem(new_prob.u0, new_prob.tspan, new_prob.p)
    sol_type = _get_sol_type(
        single_prob, zero(eltype(single_prob.tspan)), Val(SaveFields), Val(SaveWork)
    )
    local_sols = Vector{sol_type}(undef, 1)
    _adaptive_boris!(
        local_sols, single_prob, 1:1, alg, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
    )
    return local_sols[1]
end

function _solve(
        ::EnsembleDistributed, prob, trajectories, alg::AdaptiveBoris, savestepinterval,
        isoutofdomain, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        batch_size
    ) where {SaveFields, SaveWork}
    return pmap(1:trajectories; batch_size = batch_size) do i
        _solve_single_adaptive_boris(
            prob, i, alg, savestepinterval, isoutofdomain,
            save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
        )
    end
end

function _solve(
        ::EnsembleSplitThreads, prob, trajectories, alg::AdaptiveBoris, savestepinterval,
        isoutofdomain, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        batch_size
    ) where {SaveFields, SaveWork}
    # _solve_single_adaptive_boris wraps each trajectory in a fresh TraceProblem(u0, tspan, p)
    # with DEFAULT_PROB_FUNC. We get a sample problem from prob_func to ensure the
    # u0 type (uType) matches.
    sample_prob = prob.prob_func(prob, 1, false)
    dummy_prob = TraceProblem(sample_prob.u0, sample_prob.tspan, sample_prob.p)
    sol_type = _get_sol_type(
        dummy_prob, zero(eltype(dummy_prob.tspan)), Val(SaveFields), Val(SaveWork)
    )
    ichunks = index_chunks(1:trajectories; size = batch_size)
    results = pmap(ichunks) do irange
        local_sols = Vector{sol_type}(undef, length(irange))
        Threads.@threads for k in eachindex(irange)
            i = irange[k]
            local_sols[k] = _solve_single_adaptive_boris(
                prob, i, alg, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
            )
        end
        local_sols
    end
    return reduce(vcat, results)
end


@muladd function _adaptive_boris!(
        sols, prob, irange, alg, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
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

    # Pre-calculate common factors for time step calculation
    # dt = safety * 2π / (abs(q2m) * Bmag)
    C = (2π * alg.safety * sign(tspan[2] - tspan[1])) / abs(q2m)

    @fastmath @inbounds for i in irange
        # Initialize solution containers
        initial_capacity = 1000
        traj = Vector{SVector{vars_dim, T}}(undef, 0)
        tsave = Vector{typeof(tspan[1])}(undef, 0)
        sizehint!(traj, initial_capacity)
        sizehint!(tsave, initial_capacity)

        new_prob = prob.prob_func(prob, i, false)
        # Load independent r and v SVector from u0
        u0_i = SVector{6, T}(new_prob.u0)
        r = u0_i[SVector(1, 2, 3)]
        v = u0_i[SVector(4, 5, 6)]
        t = tspan[1]
        ttotal = tspan[2] - tspan[1]

        if save_start
            data = _prepare_saved_data(u0_i, p, t, Val(SaveFields), Val(SaveWork))
            push!(traj, data)
            push!(tsave, t)
        end

        # Initial dt calculation
        Bmag = norm(Bfunc(r, t))
        dt = C / Bmag

        # Backstep velocity: v(0) -> v(-1/2) using dt
        v = update_velocity(v, r, -0.5 * dt, t, p)
        it = 1
        should_save_final = save_end
        while abs(t - tspan[1]) < abs(ttotal)
            # Check if next step exceeds tspan[2]
            if abs(t + dt - tspan[1]) > abs(ttotal)
                dt_step = tspan[2] - t
                # Resync v from `t - 0.5*dt` to `t - 0.5*dt_step`
                v = update_velocity(v, r, 0.5 * dt, t, p)
                v = update_velocity(v, r, -0.5 * dt_step, t, p)
                dt = dt_step
            end

            # Saving logic (start of step)
            if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
                # Advance to t to get v_n
                v_save = update_velocity(v, r, 0.5 * dt, t, p)

                xv_s = vcat(r, v_save)
                data = _prepare_saved_data(xv_s, p, t, Val(SaveFields), Val(SaveWork))
                push!(traj, data)
                push!(tsave, t)
            end

            # Update velocity to v_{n+1/2}
            t_mid = t + 0.5 * dt
            v_new = update_velocity(v, r, dt, t_mid, p)

            # Update location x_{n} -> x_{n+1}
            r_next = r + v_new * dt
            t_next = t + dt

            xv_new = vcat(r_next, v_new)
            if isoutofdomain(xv_new, p, t_next)
                should_save_final = true
                break
            end

            r = r_next
            t = t_next
            v = v_new

            # New dt
            Bmag = norm(Bfunc(r, t))
            dt_new = C / Bmag

            # Resync v_{n+1/2}(dt) to v_{n+1/2}(dt_new)
            # v is at t_{new} - 0.5 * dt_old (relative to t_{new})
            # i.e. it is v_{n+1/2} from step we just took.
            v = update_velocity(v, r, 0.5 * dt, t, p)
            v = update_velocity(v, r, -0.5 * dt_new, t, p)
            dt = dt_new
            it += 1
            if save_everystep && (it - 1) % savestepinterval == 0
                should_save_final = true
            end
        end

        if should_save_final && (isempty(tsave) || tsave[end] != t)
            # v is at t - 0.5*dt. To get v at t, advance by 0.5*dt
            v_final = update_velocity(v, r, 0.5 * dt, t, p)

            xv_f = vcat(r, v_final)
            data = _prepare_saved_data(xv_f, p, t, Val(SaveFields), Val(SaveWork))
            push!(traj, data)
            push!(tsave, t)
        end

        # Construct solution
        sol_alg = :adaptive_boris
        interp = LinearInterpolation(tsave, traj)
        retcode = ReturnCode.Default
        stats = nothing

        sols[i] = build_solution(
            prob, sol_alg, tsave, traj; interp = interp, retcode = retcode, stats = stats
        )
    end

    return
end

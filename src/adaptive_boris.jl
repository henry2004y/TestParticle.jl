# Adaptive Boris method

struct AdaptiveBoris{T}
    dtmin::T
    dtmax::T
    safety::T
end

"""
    AdaptiveBoris(; dtmin, dtmax, safety=0.1)

Adaptive Boris method with adaptive time stepping based on local gyrofrequency.
The time step is determined by `dt = safety / |qB/m|`, clamped by `dtmin` and `dtmax`.
"""
function AdaptiveBoris(; dtmax, dtmin = 1.0e-2 * dtmax, safety = 0.1)
    T = promote_type(typeof(dtmin), typeof(dtmax), typeof(safety))
    return AdaptiveBoris{T}(T(dtmin), T(dtmax), T(safety))
end

function solve(
        prob::TraceProblem, alg::AdaptiveBoris,
        ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1,
        isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        save_fields::Bool = false, save_work::Bool = false
    )
    return if save_fields
        if save_work
            return _solve(
                ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, Val(true), Val(true)
            )
        else
            return _solve(
                ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, Val(true), Val(false)
            )
        end
    else
        if save_work
            return _solve(
                ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, Val(false), Val(true)
            )
        else
            return _solve(
                ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, Val(false), Val(false)
            )
        end
    end
end

function _solve(
        ::EnsembleSerial, prob, trajectories, alg::AdaptiveBoris, savestepinterval,
        isoutofdomain, save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}
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
        isoutofdomain, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
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
        isoutofdomain, save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    return pmap(1:trajectories) do i
        _solve_single_adaptive_boris(
            prob, i, alg, savestepinterval, isoutofdomain,
            save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork)
        )
    end
end

function _adaptive_boris!(
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

        if save_start
            data = _prepare_saved_data(u0_i, p, t, Val(SaveFields), Val(SaveWork))
            push!(traj, data)
            push!(tsave, t)
        end

        # Initial dt calculation
        B = Bfunc(r, t)
        B_mag = norm(B)
        omega = abs(q2m) * B_mag
        dt = alg.safety / omega
        dt = clamp(dt, alg.dtmin, alg.dtmax)

        # Backstep velocity: v(0) -> v(-1/2) using dt
        v = update_velocity(v, r, p, -0.5 * dt, t)

        it = 1
        while t < tspan[2]
            # Check if next step exceeds tspan[2]
            if t + dt > tspan[2]
                dt_step = tspan[2] - t
                # Resync v from `t - 0.5*dt` to `t - 0.5*dt_step`
                v = update_velocity(v, r, p, 0.5 * dt, t)
                v = update_velocity(v, r, p, -0.5 * dt_step, t)
                dt = dt_step
            end

            v_prev = v # v_{n-1/2} relative to current dt

            # Saving logic (start of step)
            if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
                # Advance to t to get v_n
                v_save = update_velocity(v_prev, r, p, 0.5 * dt, t)

                xv_s = vcat(r, v_save)
                data = _prepare_saved_data(xv_s, p, t, Val(SaveFields), Val(SaveWork))
                push!(traj, data)
                push!(tsave, t)
            end

            # Update velocity to v_{n+1/2}
            t_mid = t + 0.5 * dt
            v = update_velocity(v, r, p, dt, t_mid)

            # Update location x_{n} -> x_{n+1}
            r += v * dt
            t += dt

            isoutofdomain(vcat(r, v), p, t) && break
            it += 1

            # Calculate new dt for next step
            B = Bfunc(r, t)
            B_mag = norm(B)
            omega = abs(q2m) * B_mag
            dt_new = if omega > 0
                clamp(alg.safety / omega, alg.dtmin, alg.dtmax)
            else
                alg.dtmax
            end

            # Resync for next step
            # v is at t_{new} - 0.5 * dt_old (relative to t_{new})
            # i.e. it is v_{n+1/2} from step we just took.
            v = update_velocity(v, r, p, 0.5 * dt, t)
            v = update_velocity(v, r, p, -0.5 * dt_new, t)

            dt = dt_new
        end

        # Final save
        should_save_final = false
        if save_end
            should_save_final = true
        elseif save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
            should_save_final = true
        end

        if should_save_final
            # xv currently has x at t (which is >= tspan[2] or boundary)
            # xv[4:6] has v at t - 0.5*dt (start of next step)
            # We want v at t, so we just need to advance by 0.5 * dt
            v_final = update_velocity(v, r, p, 0.5 * dt, t)

            xv_s = vcat(r, v_final)
            data = _prepare_saved_data(xv_s, p, t, Val(SaveFields), Val(SaveWork))
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

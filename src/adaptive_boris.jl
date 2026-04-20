# Adaptive Boris inner loop

const AdaptiveBoris = Boris{true}

"""
    AdaptiveBoris(; safety=0.1)

Adaptive Boris method with adaptive time stepping based on
local gyroperiod.
The time step is determined by
`dt = safety * T_gyro = safety * 2π / |qB/m|`.
"""
AdaptiveBoris(; safety = 0.1) = Boris{true}(Float64(safety))

@inline function _adaptive_boris_single(
        prob, i, alg::Boris{true}, savestepinterval, isoutside, save_start, save_end,
        save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
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

    # dt = safety * 2π / (abs(q2m) * Bmag)
    C = (2π * alg.safety * sign(tspan[2] - tspan[1])) / abs(q2m)

    initial_capacity = 1000
    traj = Vector{SVector{vars_dim, T}}(undef, 0)
    tsave = Vector{typeof(tspan[1])}(undef, 0)
    sizehint!(traj, initial_capacity)
    sizehint!(tsave, initial_capacity)

    new_prob = prob.prob_func(prob, i, false)
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
    v = update_velocity(v, r, -0.5 * dt, t, p)
    it = 1
    should_save_final = save_end
    retcode = ReturnCode.Success

    @fastmath while abs(t - tspan[1]) < abs(ttotal)
        if abs(t + dt - tspan[1]) > abs(ttotal)
            dt_step = tspan[2] - t
            v = update_velocity(v, r, 0.5 * dt, t, p)
            v = update_velocity(v, r, -0.5 * dt_step, t, p)
            dt = dt_step
        end

        if save_everystep &&
                (it - 1) > 0 &&
                (it - 1) % savestepinterval == 0
            v_save = update_velocity(
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
        v_new = update_velocity(v, r, dt, t_mid, p)

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
        v = update_velocity(v, r, 0.5 * dt, t, p)
        v = update_velocity(v, r, -0.5 * dt_new, t, p)
        dt = dt_new
        it += 1
        if save_everystep && (it - 1) % savestepinterval == 0
            should_save_final = true
        end
    end

    if should_save_final && (isempty(tsave) || tsave[end] != t)
        v_final = update_velocity(v, r, 0.5 * dt, t, p)
        xv_f = vcat(r, v_final)
        data = _prepare_saved_data(
            xv_f, p, t,
            Val(SaveFields), Val(SaveWork)
        )
        push!(traj, data)
        push!(tsave, t)
    end

    sol_alg = :adaptive_boris
    interp = LinearInterpolation(tsave, traj)
    stats = nothing

    return build_solution(prob, sol_alg, tsave, traj; interp, retcode, stats)
end

@muladd function _adaptive_boris!(
        sols, prob, irange, alg::Boris{true},
        savestepinterval, isoutside,
        save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}

    @inbounds for i in irange
        sols[i] = _adaptive_boris_single(
            prob, i, alg, savestepinterval, isoutside, save_start, save_end, save_everystep,
            Val(SaveFields), Val(SaveWork)
        )
    end

    return
end

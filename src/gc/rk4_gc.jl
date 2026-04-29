"""
    update_rk4(y, param, dt, t)

Update state using the RK4 method.
Returns dy as the update SArray.
"""
@muladd function update_rk4(y, param, dt, t)
    k1 = trace_gc(y, param, t)

    y_tmp = y + 0.5 * dt * k1
    k2 = trace_gc(y_tmp, param, t + 0.5 * dt)

    y_tmp = y + 0.5 * dt * k2
    k3 = trace_gc(y_tmp, param, t + 0.5 * dt)

    y_tmp = y + dt * k3
    k4 = trace_gc(y_tmp, param, t + dt)

    dy = (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return dy
end

"""
Apply RK4 method for particles with index in `irange`.
"""
function _rk4!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutside,
        save_start, save_end, save_everystep, maxiters,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    (; tspan, p, u0) = prob
    ttotal = tspan[2] - tspan[1]
    T = eltype(u0)

    vars_dim = 4
    if SaveFields
        vars_dim += 6
    end
    if SaveWork
        vars_dim += 4
    end

    @fastmath @inbounds for i in irange
        traj = Vector{SVector{vars_dim, T}}(undef, nout)
        tsave = Vector{typeof(tspan[1] + dt)}(undef, nout)

        # set initial conditions for each trajectory i
        iout = 0
        new_prob = prob.prob_func(prob, (sim_id = i, repeat = false))
        xv = new_prob.u0

        if save_start
            iout += 1
            traj[iout] = _prepare_saved_data_gc(
                xv, p, tspan[1], Val(SaveFields), Val(SaveWork)
            )
            tsave[iout] = tspan[1]
        end

        it = 1

        while it <= nt && it <= maxiters
            t = tspan[1] + (it - 1) * dt

            dx = update_rk4(xv, p, dt, t)
            xv_next = xv + dx

            if isoutside(xv_next, p, t + dt)
                break
            end
            xv = xv_next

            if save_everystep && (it % savestepinterval == 0)
                iout += 1
                if iout <= nout
                    traj[iout] = _prepare_saved_data_gc(
                        xv, p, t + dt, Val(SaveFields), Val(SaveWork)
                    )
                    tsave[iout] = t + dt
                end
            end

            it += 1
        end

        # Handle save_end logic
        final_step = it - 1
        should_save_final = false
        if save_end
            should_save_final = true
        elseif save_everystep && (final_step > 0) && (final_step % savestepinterval == 0)
            should_save_final = true
        end

        if iout < nout && should_save_final
            t_final = (final_step == nt) ? tspan[2] : (tspan[1] + final_step * dt)
            if iout == 0 || tsave[iout] < t_final
                iout += 1
                traj[iout] = _prepare_saved_data_gc(
                    xv, p, t_final, Val(SaveFields), Val(SaveWork)
                )
                tsave[iout] = t_final
            end
        end

        retcode = if it <= nt && it <= maxiters
            ReturnCode.Terminated
        elseif it > maxiters
            ReturnCode.MaxIters
        else
            ReturnCode.Success
        end

        if iout < nout
            resize!(traj, iout)
            resize!(tsave, iout)
        end

        alg = :rk4
        t = tsave
        interp = LinearInterpolation(t, traj)
        stats = nothing

        sols[i] = build_solution(prob, alg, t, traj; interp, retcode, stats)
    end

    return
end

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
        sols, prob, irange, plan, dt, nt, nout, isoutside,
        save_start, save_end, save_everystep, maxiters,
        ::Val{SaveFields}, ::Val{SaveWork}, seed = nothing
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
        rng = isnothing(seed) ? default_rng() : Xoshiro(_rng_seed(seed, i))
        new_prob = prob.prob_func(prob, EnsembleContext(i, 1, 0, nothing, rng, seed))
        xv = new_prob.u0

        if save_start
            iout += 1
            traj[iout] = _prepare_saved_data_gc(
                xv, p, tspan[1], Val(SaveFields), Val(SaveWork)
            )
            tsave[iout] = tspan[1]
        end

        it = 1
        nsave = length(plan.times)
        isave = 1

        while it <= nt && it <= maxiters
            t = tspan[1] + (it - 1) * dt
            t_next = t + dt

            xv_prev = xv
            dx = update_rk4(xv, p, dt, t)
            xv_next = xv + dx

            if isoutside(xv_next, p, t_next)
                break
            end
            xv = xv_next

            if use_saveat(plan)
                # Report every requested time this step has passed, interpolating
                # inside the step so the integration itself is untouched.
                while isave <= nsave && _saveat_reached(plan.times[isave], t_next, plan.dir)
                    t_target = plan.times[isave]
                    iout += 1
                    if iout <= nout
                        traj[iout] = _prepare_saved_data_gc(
                            _saveat_interpolate(t, xv_prev, t_next, xv, t_target),
                            p, t_target, Val(SaveFields), Val(SaveWork)
                        )
                        tsave[iout] = t_target
                    end
                    isave += 1
                end
            elseif save_everystep && (it % plan.interval == 0)
                iout += 1
                if iout <= nout
                    traj[iout] = _prepare_saved_data_gc(
                        xv, p, t_next, Val(SaveFields), Val(SaveWork)
                    )
                    tsave[iout] = t_next
                end
            end

            it += 1
        end

        # Handle save_end logic
        final_step = it - 1
        should_save_final = false
        if save_end
            should_save_final = true
        elseif !use_saveat(plan) && save_everystep && (final_step > 0) &&
                (final_step % plan.interval == 0)
            should_save_final = true
        end

        if iout < nout && should_save_final
            t_final = (final_step == nt) ? tspan[2] : (tspan[1] + final_step * dt)
            if iout == 0 || plan.dir * (t_final - tsave[iout]) > 0
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

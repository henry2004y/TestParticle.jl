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

# Steps to take, capped by `maxiters`, steps the span needs, and the saved state count.
function _rk4_steps(prob, dt, plan, save_start, save_end, save_everystep, maxiters)
    isnothing(dt) && throw(ArgumentError("GCRK4 needs a time step `dt`."))
    ttotal = prob.tspan[2] - prob.tspan[1]
    sign(dt) != sign(ttotal) && ttotal != 0 && throw(
        ArgumentError("dt must have the same sign as tspan[2] - tspan[1] for GCRK4.")
    )
    nsteps = round(Int, ttotal / dt)
    nt = trunc(Int, min(nsteps, maxiters))

    nout = save_start ? 1 : 0
    if use_saveat(plan)
        # One slot per requested time, plus the end of the run. A `maxiters` cut
        # short can leave some of them unreached, which the caller resizes away.
        nout += length(plan.times) + (save_end ? 1 : 0)
    elseif save_everystep
        nout += nt
        nt > 0 && !save_end && (nout -= 1)
        nt == 0 && save_end && (nout += 1)
    elseif save_end
        nout += 1
    end

    return nt, nsteps, nout
end

function _trace_gc_native(
        alg::GCRK4, prob, dt, plan, isoutside::F, save_start, save_end, save_everystep,
        maxiters, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {F, SaveFields, SaveWork}
    (; tspan, p, u0) = prob
    nt, nsteps, nout = _rk4_steps(
        prob, dt, plan, save_start, save_end, save_everystep, maxiters
    )
    T = eltype(u0)
    vars_dim = 4 + (SaveFields ? 6 : 0) + (SaveWork ? 4 : 0)

    traj = Vector{SVector{vars_dim, T}}(undef, nout)
    tsave = Vector{typeof(tspan[1] + dt)}(undef, nout)

    iout = 0
    xv = SVector{4, T}(u0)

    if save_start
        iout += 1
        traj[iout] = _prepare_saved_data_gc(xv, p, tspan[1], Val(SaveFields), Val(SaveWork))
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
            while isave <= nsave && saveat_reached(plan.times[isave], t_next, plan.dir)
                t_target = plan.times[isave]
                iout += 1
                if iout <= nout
                    traj[iout] = _prepare_saved_data_gc(
                        saveat_interpolate(t, xv_prev, t_next, xv, t_target),
                        p, t_target, Val(SaveFields), Val(SaveWork)
                    )
                    tsave[iout] = t_target
                end
                isave += 1
            end
        elseif save_everystep
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

    final_step = it - 1
    should_save_final = save_end ||
        (!use_saveat(plan) && save_everystep && final_step > 0)

    if iout < nout && should_save_final
        t_final = (final_step == nsteps) ? tspan[2] : (tspan[1] + final_step * dt)
        if iout == 0 || plan.dir * (t_final - tsave[iout]) > 0
            iout += 1
            traj[iout] = _prepare_saved_data_gc(
                xv, p, t_final, Val(SaveFields), Val(SaveWork)
            )
            tsave[iout] = t_final
        end
    end

    retcode = if it <= nt
        ReturnCode.Terminated
    elseif nt < nsteps
        ReturnCode.MaxIters
    else
        ReturnCode.Success
    end

    if iout < nout
        resize!(traj, iout)
        resize!(tsave, iout)
    end

    interp = LinearInterpolation(tsave, traj)

    return build_solution(prob, alg, tsave, traj; interp, retcode, stats = nothing)
end

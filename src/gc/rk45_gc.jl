"""
    update_dp5(y, param, dt, t, mp = nothing)::(dy, E)

Update state using the Dormand-Prince 5(4) method.
Returns dy as the update and E as the error estimate, both as SArrays.

`mp` is the [`get_magnetic_properties`](@ref) output at `y` and `t`, when the
caller already has it, so the first stage does not recompute the derivatives.
"""
@muladd function update_dp5(y, param, dt, t, mp = nothing)
    # Step in the state precision so Float32 states stay Float32
    T = eltype(y)
    dt = T(dt)
    # Coefficients for DP5
    c2, c3, c4, c5, c6 = T.((1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0))
    a21 = T(1 / 5)
    a31, a32 = T.((3 / 40, 9 / 40))
    a41, a42, a43 = T.((44 / 45, -56 / 15, 32 / 9))
    a51, a52, a53, a54 = T.((19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729))
    a61, a62, a63, a64, a65 = T.((9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656))
    a71, a72, a73, a74, a75, a76 = T.((35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84))

    # Coefficients for 5th order solution and error estimate
    # b2 = b7 = e2 = 0.0
    b1, b3, b4, b5, b6 = T.((35 / 384, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84))
    e1, e3, e4, e5, e6, e7 = T.((71 / 57600, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40))

    k1 = mp === nothing ? trace_gc(y, param, t) : _trace_gc(y, param, t, mp)

    y_tmp = y + dt * (a21 * k1)
    k2 = trace_gc(y_tmp, param, t + c2 * dt)

    y_tmp = y + dt * (a31 * k1 + a32 * k2)
    k3 = trace_gc(y_tmp, param, t + c3 * dt)

    y_tmp = y + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = trace_gc(y_tmp, param, t + c4 * dt)

    y_tmp = y + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    k5 = trace_gc(y_tmp, param, t + c5 * dt)

    y_tmp = y + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = trace_gc(y_tmp, param, t + c6 * dt)

    y_tmp = y + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    k7 = trace_gc(y_tmp, param, t + dt)

    # y_{n+1} update
    dy = dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)
    E = dt * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7)

    return dy, E
end

"""
Apply RK45 method for particles with index in `irange`.
"""
function _rk45!(
        sols, prob, irange, alg, plan, dt_initial, isoutside,
        save_start, save_end, save_everystep, maxiters,
        ::Val{SaveFields}, ::Val{SaveWork}, seed = nothing
    ) where {SaveFields, SaveWork}
    (; tspan, p, u0) = prob
    (; abstol, reltol) = alg
    T = eltype(u0)
    Tt = typeof(tspan[1] + one(T))

    safety = 0.9
    max_growth = 5.0
    min_growth = 0.2

    vars_dim = 4
    if SaveFields
        vars_dim += 6
    end
    if SaveWork
        vars_dim += 4
    end

    @inbounds for i in irange
        traj = SVector{vars_dim, T}[]
        tsave = Tt[]

        rng = isnothing(seed) ? default_rng() : Xoshiro(_rng_seed(seed, i))
        new_prob = prob.prob_func(prob, EnsembleContext(i, 1, 0, nothing, rng, seed))
        xv = SVector{4, T}(new_prob.u0)

        t = Tt(tspan[1])

        if isnothing(dt_initial)
            # Estimate initial step size based on gyroperiod
            q2m = p.q2m
            B_vec = p.B(get_x(xv), t)
            Bmag = norm(B_vec)

            if Bmag == 0 || q2m == 0
                dt = Tt(1.0e-6)
            else
                omega = abs(q2m * Bmag)
                dt = Tt(0.5 * 2π / omega) # 0.5 gyroperiod
            end
        else
            dt = Tt(dt_initial)
        end

        if save_start
            push!(traj, _prepare_saved_data_gc(xv, p, t, Val(SaveFields), Val(SaveWork)))
            push!(tsave, t)
        end

        steps = 0
        nsave = length(plan.times)
        isave = 1
        retcode = ReturnCode.Success
        while t < tspan[2] && steps < maxiters
            if t + dt > tspan[2]
                dt = tspan[2] - t
            end

            dx, E = update_dp5(xv, p, dt, t)

            error_ratio = 0.0

            y_next = xv + dx
            sum_sq_error = 0.0
            for k in 1:4
                sc = abstol + max(abs(xv[k]), abs(y_next[k])) * reltol
                sum_sq_error += (E[k] / sc)^2
            end
            error_ratio = 0.5 * sqrt(sum_sq_error)

            # A non-finite error estimate means the trial step hit invalid fields;
            # either stop on the domain check or reject the step and shrink `dt`.
            if !isfinite(error_ratio)
                if isoutside(y_next, p, t + dt)
                    retcode = ReturnCode.Terminated
                    break
                end
                error_ratio = Inf
            end

            if error_ratio <= 1.0
                y_next = xv + dx

                if isoutside(y_next, p, t + dt)
                    retcode = ReturnCode.Terminated
                    break
                end

                t_prev = t
                xv_prev = xv
                t += dt
                xv = y_next

                if use_saveat(plan)
                    while isave <= nsave && saveat_reached(plan.times[isave], t, plan.dir)
                        t_target = plan.times[isave]
                        push!(
                            traj, _prepare_saved_data_gc(
                                saveat_interpolate(t_prev, xv_prev, t, xv, t_target),
                                p, t_target, Val(SaveFields), Val(SaveWork)
                            )
                        )
                        push!(tsave, t_target)
                        isave += 1
                    end
                elseif save_everystep
                    push!(
                        traj, _prepare_saved_data_gc(
                            xv, p, t, Val(SaveFields), Val(SaveWork)
                        )
                    )
                    push!(tsave, t)
                end

                steps += 1
            end

            scale = if error_ratio == 0.0
                max_growth
            else
                safety * (1.0 / error_ratio)^0.2
            end
            scale = max(min_growth, min(scale, max_growth))
            dt *= Tt(scale)

            dt < 1.0e-14 && break
        end

        if save_end && (isempty(tsave) || tsave[end] != t)
            push!(traj, _prepare_saved_data_gc(xv, p, t, Val(SaveFields), Val(SaveWork)))
            push!(tsave, t)
        end

        t_final = tsave
        u_final = traj
        interp = LinearInterpolation(t_final, u_final)
        stats = nothing

        sols[i] = build_solution(prob, alg, t_final, u_final; interp, retcode, stats)
    end
    return
end

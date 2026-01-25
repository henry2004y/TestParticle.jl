"""
    update_velocity_multistep(v, r, param, dt, t, n)

Update velocity using the Multistep Boris method, returning the new velocity as an SVector.
Reference: [Zenitani & Kato 2025](https://arxiv.org/abs/2505.02270)
"""
@muladd function update_velocity_multistep(v, r, param, dt, t, n::Int)
    q2m, _, Efunc, Bfunc, _ = param
    E = Efunc(r, t)
    B = Bfunc(r, t)

    # t_n and e_n vectors
    factor = q2m * dt / (2 * n)

    t_n = factor * B # (q/m * dt/(2n)) * B
    e_n = factor * E # (q/m * dt/(2n)) * E

    t_n_mag2 = sum(abs2, t_n)
    t_n_mag = sqrt(t_n_mag2)

    # Calculate coefficients
    # Check for small t_n to avoid division by zero or precision loss
    if t_n_mag < 1.0e-4
        # Taylor expansion limits as t_n -> 0
        c_n1 = 1 - 2 * n * n * t_n_mag2
        c_n2 = 2 * n - (4 / 3) * n * n * n * t_n_mag2
        c_n3 = 2 * n * n
        c_n6 = (4 / 3) * n * n * n
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

function _multistep_boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n_steps::Int,
        save_start, save_end, save_everystep, ::Val{SaveFields}
    ) where {SaveFields}
    (; tspan, p, u0) = prob
    q2m, _, Efunc, Bfunc, _ = p
    T = eltype(u0)
    vars_dim = SaveFields ? 12 : 6

    @fastmath @inbounds for i in irange
        traj = Vector{SVector{vars_dim, T}}(undef, nout)
        tsave = Vector{typeof(tspan[1] + dt)}(undef, nout)

        # set initial conditions for each trajectory i
        iout = 0
        new_prob = prob.prob_func(prob, i, false)
        # Load independent r and v SVector from u0
        u0_i = SVector{6, T}(new_prob.u0)
        r = u0_i[SVector(1, 2, 3)]
        v = u0_i[SVector(4, 5, 6)]

        if save_start
            iout += 1
            if SaveFields
                E = SVector{3, T}(Efunc(u0_i, tspan[1]))
                B = SVector{3, T}(Bfunc(u0_i, tspan[1]))
                traj[iout] = vcat(u0_i, E, B)
            else
                traj[iout] = u0_i
            end
            tsave[iout] = tspan[1]
        end

        # push velocity back in time by 1/2 dt
        v = update_velocity_multistep(v, r, p, -0.5 * dt, tspan[1], n_steps)

        it = 1
        while it <= nt
            v_prev = v
            t = (it - 0.5) * dt
            v = update_velocity_multistep(v, r, p, dt, t, n_steps)

            if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
                iout += 1
                if iout <= nout
                    t_current = tspan[1] + (it - 1) * dt
                    v_save = update_velocity_multistep(
                        v_prev, r, p, 0.5 * dt,
                        t_current, n_steps
                    ) # v at t_current

                    xv_s = vcat(r, v_save)

                    if SaveFields
                        E = SVector{3, T}(Efunc(xv_s, t_current))
                        B = SVector{3, T}(Bfunc(xv_s, t_current))
                        traj[iout] = vcat(xv_s, E, B)
                    else
                        traj[iout] = xv_s
                    end
                    tsave[iout] = t_current
                end
            end

            r += v * dt
            isoutofdomain(vcat(r, v), p, it * dt) && break
            it += 1
        end

        final_step = min(it, nt)
        should_save_final = false
        if save_end
            should_save_final = true
        elseif save_everystep && (final_step > 0) && (final_step % savestepinterval == 0)
            should_save_final = true
        end

        if iout < nout && should_save_final
            iout += 1
            t_final = final_step == nt ? tspan[2] : tspan[1] + final_step * dt
            dt_final = t_final - (tspan[1] + (final_step - 0.5) * dt)
            v_final = update_velocity_multistep(
                v, r, p, dt_final,
                t_final, n_steps
            )

            xv_s = vcat(r, v_final)
            if SaveFields
                E = SVector{3, T}(Efunc(xv_s, t_final))
                B = SVector{3, T}(Bfunc(xv_s, t_final))
                traj[iout] = vcat(xv_s, E, B)
            else
                traj[iout] = xv_s
            end
            tsave[iout] = t_final
        end

        if iout < nout
            resize!(traj, iout)
            resize!(tsave, iout)
        end
        traj_save = traj
        t = tsave

        alg = :multistep_boris
        interp = LinearInterpolation(t, traj_save)
        retcode = ReturnCode.Default
        stats = nothing

        sols[i] = build_solution(
            prob, alg, t, traj_save; interp = interp, retcode = retcode, stats = stats
        )
    end

    return
end

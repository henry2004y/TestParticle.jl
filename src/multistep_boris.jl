using LinearAlgebra: dot

struct MultistepBorisMethod{TV}
    t_n::TV
    e_n::TV
    v_cross_t::TV
    e_cross_t::TV
end

function MultistepBorisMethod(T::Type{<:AbstractFloat} = Float64)
    t_n = MVector{3, T}(undef)
    e_n = MVector{3, T}(undef)
    v_cross_t = MVector{3, T}(undef)
    e_cross_t = MVector{3, T}(undef)

    return MultistepBorisMethod{typeof(t_n)}(t_n, e_n, v_cross_t, e_cross_t)
end

"""
    update_velocity_multistep!(xv, paramBoris, param, dt, t, n)

Update velocity using the Multistep Boris method.
Reference: [Zenitani & Kato 2025](https://arxiv.org/abs/2505.02270)
"""
@muladd function update_velocity_multistep!(xv, paramBoris, param, dt, t, n::Int)
    (; t_n, e_n, v_cross_t, e_cross_t) = paramBoris
    q2m, _, Efunc, Bfunc, _ = param
    E = Efunc(xv, t)
    B = Bfunc(xv, t)

    # t_n and e_n vectors
    # Note: t_n in paper is (q/m * dt/(2n)) * B
    # e_n in paper is (q/m * dt/(2n)) * E
    # q2m is q/m
    factor = q2m * dt / (2 * n)

    @inbounds for dim in 1:3
        t_n[dim] = factor * B[dim]
        e_n[dim] = factor * E[dim]
    end

    t_n_mag2 = sum(abs2, t_n)
    t_n_mag = sqrt(t_n_mag2)

    # Calculate coefficients
    # Check for small t_n to avoid division by zero or precision loss
    if t_n_mag < 1.0e-4
        # Taylor expansion limits as t_n -> 0
        c_n1 = 1.0 - 2 * n * n * t_n_mag2
        c_n2 = 2 * n - (4.0 / 3.0) * n * n * n * t_n_mag2
        c_n3 = 2.0 * n * n
        c_n6 = (4.0 / 3.0) * n * n * n
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

    # Extract velocity
    v = @view xv[4:6]

    # Dot products
    v_dot_t = dot(v, t_n)
    e_dot_t = dot(e_n, t_n)

    # Cross products
    # v x t_n
    cross!(v, t_n, v_cross_t)

    # e_n x t_n
    cross!(e_n, t_n, e_cross_t)

    # Update velocity
    # Equation 39:
    # v_new = c_n1*v + c_n2*(v x t_n) + c_n3*(v . t_n)t_n + c_n4*e_n + c_n5*(e_n x t_n) + c_n6*(e_n . t_n)t_n
    @inbounds for i in 1:3
        xv[i + 3] = c_n1 * xv[i + 3] +
            c_n2 * v_cross_t[i] +
            c_n3 * v_dot_t * t_n[i] +
            c_n4 * e_n[i] +
            c_n5 * e_cross_t[i] +
            c_n6 * e_dot_t * t_n[i]
    end

    return
end

function _multistep_boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n_steps::Int,
        save_start, save_end, save_everystep, ::Val{ITD}
    ) where {ITD}
    (; tspan, p, u0) = prob
    T = eltype(u0)
    paramBoris = MultistepBorisMethod(T)
    xv = MVector{6, T}(undef)
    xv_save = MVector{6, T}(undef)
    v_old = MVector{3, T}(undef)

    @fastmath @inbounds for i in irange
        traj = Vector{SVector{6, T}}(undef, nout)
        tsave = Vector{typeof(tspan[1] + dt)}(undef, nout)

        # set initial conditions for each trajectory i
        iout = 0
        new_prob = prob.prob_func(prob, i, false)
        xv .= new_prob.u0

        if save_start
            iout += 1
            traj[iout] = SVector{6, T}(xv)
            tsave[iout] = tspan[1]
        end

        # push velocity back in time by 1/2 dt
        update_velocity_multistep!(xv, paramBoris, p, -0.5 * dt, tspan[1], n_steps)

        it = 1
        while it <= nt
            v_old .= @view xv[4:6]
            t = ITD ? (it - 0.5) * dt : zero(dt)
            update_velocity_multistep!(xv, paramBoris, p, dt, t, n_steps)

            if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
                iout += 1
                if iout <= nout
                    xv_save .= xv
                    xv_save[4:6] .= v_old
                    t_current = tspan[1] + (it - 1) * dt
                    update_velocity_multistep!(
                        xv_save, paramBoris, p, 0.5 * dt,
                        ITD ? t_current : zero(dt), n_steps
                    )
                    traj[iout] = SVector{6, T}(xv_save)
                    tsave[iout] = t_current
                end
            end

            update_location!(xv, dt)
            if isoutofdomain(xv, p, it * dt)
                break
            end
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
            update_velocity_multistep!(
                xv, paramBoris, p, dt_final,
                ITD ? t_final : zero(dt), n_steps
            )
            traj[iout] = SVector{6, T}(xv)
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

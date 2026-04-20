# Multistep Boris solver

struct MultistepUpdater{N_order}
    n::Int
end

@inline function (mu::MultistepUpdater{N_order})(v, r, dt, t, p) where {N_order}
    return update_velocity_multistep(v, r, dt, t, mu.n, Val(N_order), p)
end

function _dispatch_boris!(
        sols, prob::TraceProblem, irange,
        savestepinterval, dt, nt, nout, isoutside::F,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        alg::MultistepBoris{N, false}
    ) where {N, SaveFields, SaveWork, F}
    return _multistep_boris!(
        sols, prob, irange, savestepinterval, dt, nt,
        nout, isoutside,
        save_start, save_end, save_everystep,
        Val(SaveFields), Val(SaveWork), alg.n, Val(N)
    )
end

@inline function _multistep_boris_single(
        prob::TraceProblem, i, savestepinterval, dt, nt, nout, isoutside::F,
        save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}, alg::MultistepBoris{N_order, false}
    ) where {N_order, SaveFields, SaveWork, F}

    return _boris_single(
        prob, i, savestepinterval, dt, nt, nout, isoutside,
        save_start, save_end, save_everystep, Val(SaveFields), Val(SaveWork),
        MultistepUpdater{N_order}(alg.n), :multistep_boris
    )
end

@inline @muladd function _multistep_boris!(
        sols, prob::TraceProblem, irange, savestepinterval, dt, nt, nout, isoutside::F,
        save_start, save_end, save_everystep, ::Val{SaveFields}, ::Val{SaveWork},
        n::Int, ::Val{N_order}
    ) where {N_order, SaveFields, SaveWork, F}

    @inbounds for i in irange
        sols[i] = _multistep_boris_single(
            prob, i, savestepinterval, dt, nt, nout, isoutside,
            save_start, save_end, save_everystep,
            Val(SaveFields), Val(SaveWork), MultistepBoris{N_order, false}(n, 0.0)
        )
    end

    return
end

"""
    update_velocity_multistep(v, r, dt, t, n, ::Val{N}, param)

Update velocity using the Multistep/Hyper Boris method, returning the new velocity as an SVector.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order. When N=2, it corresponds to the Multicycle solver. When N=4 or N=6, it is the Hyper Boris solver.
Reference: [Zenitani \u0026 Kato 2025](https://arxiv.org/abs/2505.02270)
"""
@muladd function update_velocity_multistep(
        v, r, dt, t, n::Int, ::Val{N}, param
    ) where {N}
    q2m, _, Efunc, Bfunc, _ = param
    E = Efunc(r, t)
    B = Bfunc(r, t)

    # t_n and e_n vectors
    factor = q2m * dt / (2 * n)

    t_n = factor * B # (q/m * dt/(2n)) * B
    e_n = factor * E # (q/m * dt/(2n)) * E

    # Hyper Boris N-th order gyrophase correction
    if N != 2
        t_mag2 = sum(abs2, t_n)
        if N == 4
            f_N = 1 + t_mag2 / 3
            e_corr_factor = -1 / 3
        else # N == 6
            f_N = 1 + t_mag2 / 3 + 2 * t_mag2 * t_mag2 / 15
            e_corr_factor = -1 / 3 - 2 * t_mag2 / 15
        end

        e_dot_t = e_n ⋅ t_n
        e_n = f_N * e_n + (e_corr_factor * e_dot_t) * t_n
        t_n = f_N * t_n
    end

    t_n_mag2 = sum(abs2, t_n)
    t_n_mag = sqrt(t_n_mag2)

    # Calculate coefficients
    # Check for small t_n to avoid division by zero or precision loss
    if t_n_mag < TN_MAG_THRESHOLD
        # Taylor expansion limits as t_n -> 0
        c_n1 = 1 - 2 * n * n * t_n_mag2

        n_term1 = 2 * n
        n_term3 = 4 * n * n * n

        c_n2 = n_term1 - (n_term1 + n_term3) / 3 * t_n_mag2
        c_n3 = 2 * n * n - (4 * n * n + 2 * n * n * n * n) / 3 * t_n_mag2
        c_n6 = (n_term1 + n_term3) / 3
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

const TN_MAG_THRESHOLD = 1.0e-4

@inline @muladd function boris_velocity_update(v, E, B, qdt_2m)
    t_rotate = qdt_2m * B
    t_mag2 = sum(abs2, t_rotate)
    s_rotate = 2 * t_rotate / (1 + t_mag2)

    v⁻ = v + qdt_2m * E
    v′ = v⁻ + (v⁻ × t_rotate)
    v⁺ = v⁻ + (v′ × s_rotate)

    v_new = v⁺ + qdt_2m * E

    return v_new
end

@muladd function update_velocity_multistep(v, r, dt, t, n::Int, ::Val{N}, param) where {N}
    q2m, Efunc, Bfunc = get_q2m(param), get_EField(param), get_BField(param)

    E = Efunc(r, t)
    B = Bfunc(r, t)

    factor = q2m * dt / (2 * n)

    t_n = factor * B
    e_n = factor * E

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

    if t_n_mag < TN_MAG_THRESHOLD
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

    v_new = c_n1 * v +
        c_n2 * v_cross_t +
        c_n3 * v_dot_t * t_n +
        c_n4 * e_n +
        c_n5 * e_cross_t +
        c_n6 * e_dot_t * t_n

    return v_new
end

"""
    velocity_update(v, r, dt, t, p, alg)

Advance the velocity `v` by `dt`, evaluating the fields at `(r, t)`.
"""
@inline @muladd function velocity_update(v, r, dt, t, p, ::Union{Boris, AdaptiveBoris})
    qdt_2m = get_q2m(p) * 0.5 * dt
    return boris_velocity_update(v, get_EField(p)(r, t), get_BField(p)(r, t), qdt_2m)
end

@inline @muladd function velocity_update(
        v, r, dt, t, p, alg::Union{MultistepBoris{N}, AdaptiveMultistepBoris{N}}
    ) where {N}
    return update_velocity_multistep(v, r, dt, t, alg.n, Val{N}(), p)
end

# The velocity is carried at the half step, as in a leapfrog scheme: the cache
# holds `v(t - dt/2)` and the node velocity is reconstructed only for output.
# Changing `dt` re-centres the stored velocity onto the new half step, which is
# what keeps the scheme time-reversible under adaptive stepping.
@inline @muladd function boris_initialize!(integrator, cache)
    t = integrator.t
    dt = integrator.dt
    p = integrator.p
    uprev = integrator.uprev
    r = SVector(uprev[1], uprev[2], uprev[3])
    v = SVector(uprev[4], uprev[5], uprev[6])

    cache.v_half = velocity_update(v, r, -0.5 * dt, t, p, integrator.alg)
    cache.dt_prev = dt

    integrator.kshortsize = 0
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    return
end

@inline @muladd function boris_advance!(integrator, cache)
    t = integrator.t
    dt = integrator.dt
    p = integrator.p
    alg = integrator.alg
    uprev = integrator.uprev
    r = SVector(uprev[1], uprev[2], uprev[3])

    if cache.dt_prev != dt
        v_node = velocity_update(cache.v_half, r, 0.5 * cache.dt_prev, t, p, alg)
        cache.v_half = velocity_update(v_node, r, -0.5 * dt, t, p, alg)
    end

    v_half = velocity_update(cache.v_half, r, dt, t + 0.5 * dt, p, alg)
    r_new = r + v_half * dt
    v_new = velocity_update(v_half, r_new, 0.5 * dt, t + dt, p, alg)

    cache.v_half = v_half
    cache.dt_prev = dt

    return r_new, v_new
end

function initialize!(integrator, cache::BorisConstantCache)
    return boris_initialize!(integrator, cache)
end

function initialize!(integrator, cache::BorisCache)
    return boris_initialize!(integrator, cache)
end

function initialize!(integrator, cache::MultistepBorisConstantCache)
    return boris_initialize!(integrator, cache)
end

function initialize!(integrator, cache::MultistepBorisCache)
    return boris_initialize!(integrator, cache)
end

@muladd function perform_step!(integrator, cache::BorisConstantCache, repeat_step = false)
    r_new, v_new = boris_advance!(integrator, cache)
    integrator.u = vcat(r_new, v_new)
    return integrator.u
end

@muladd function perform_step!(integrator, cache::BorisCache, repeat_step = false)
    r_new, v_new = boris_advance!(integrator, cache)
    integrator.u[1] = r_new[1]
    integrator.u[2] = r_new[2]
    integrator.u[3] = r_new[3]
    integrator.u[4] = v_new[1]
    integrator.u[5] = v_new[2]
    integrator.u[6] = v_new[3]
    return
end

@muladd function perform_step!(integrator, cache::MultistepBorisConstantCache, repeat_step = false)
    r_new, v_new = boris_advance!(integrator, cache)
    integrator.u = vcat(r_new, v_new)
    return integrator.u
end

@muladd function perform_step!(integrator, cache::MultistepBorisCache, repeat_step = false)
    r_new, v_new = boris_advance!(integrator, cache)
    integrator.u[1] = r_new[1]
    integrator.u[2] = r_new[2]
    integrator.u[3] = r_new[3]
    integrator.u[4] = v_new[1]
    integrator.u[5] = v_new[2]
    integrator.u[6] = v_new[3]
    return
end

const TN_MAG_THRESHOLD = 1.0e-4

@inline @muladd function _boris_rotate(u⁻, t_rotate)
    T = eltype(u⁻)
    t_mag2 = sum(abs2, t_rotate)
    s_rotate = (2 / (one(T) + t_mag2)) * t_rotate
    u′ = u⁻ + (u⁻ × t_rotate)
    return u⁻ + (u′ × s_rotate)
end

@inline @muladd function update_velocity_boris(v, E, B, qdt2m)
    qE2m = qdt2m * E
    return _boris_rotate(v + qE2m, qdt2m * B) + qE2m
end

@inline @muladd function _quadratic_gamma(σ::T, τ2_plus_ustar2::T) where {T}
    d = sqrt(σ * σ + 4 * τ2_plus_ustar2)
    # Rationalized branch avoids catastrophic cancellation when σ < 0 (large B * dt).
    γ2 = ifelse(σ >= 0, T(0.5) * (σ + d), 2 * τ2_plus_ustar2 / (d - σ))
    return sqrt(γ2)
end

@inline @muladd function update_velocity_vay(u, E, B, qdt2m, inv_c2)
    T = eltype(u)
    γ = sqrt(one(T) + sum(abs2, u) * inv_c2)
    v = u / γ
    τ = qdt2m * B
    u′ = u + (2 * qdt2m) * E + (v × τ)

    τ2 = sum(abs2, τ)
    u_star2 = abs2(u′ ⋅ τ) * inv_c2
    σ = one(T) + sum(abs2, u′) * inv_c2 - τ2
    γ_new = _quadratic_gamma(σ, τ2 + u_star2)

    t_vec = τ / γ_new
    s = inv(one(T) + sum(abs2, t_vec))
    return s * (u′ + (u′ ⋅ t_vec) * t_vec + (u′ × t_vec))
end

@inline @muladd function update_velocity_higuera_cary(u, E, B, qdt2m, inv_c2)
    T = eltype(u)
    qE2m = qdt2m * E
    u⁻ = u + qE2m
    τ = qdt2m * B

    τ2 = sum(abs2, τ)
    u_star2 = abs2(u⁻ ⋅ τ) * inv_c2
    σ = one(T) + sum(abs2, u⁻) * inv_c2 - τ2
    γ_hc = _quadratic_gamma(σ, τ2 + u_star2)

    return _boris_rotate(u⁻, τ / γ_hc) + qE2m
end

@inline @muladd function _hyper_boris_factors(::Type{T}, t_mag2, ::Val{N}) where {T, N}
    if N == 4
        return one(T) + t_mag2 / T(3), -one(T) / T(3)
    elseif N == 6
        return one(T) + t_mag2 / T(3) + T(2) * t_mag2 * t_mag2 / T(15),
            -one(T) / T(3) - T(2) * t_mag2 / T(15)
    else
        return one(T), zero(T)
    end
end

@inline @muladd function _multistep_rotation_coeffs(
        t_n::SVector{3, T}, n::Int
    ) where {T}
    t_n_mag2 = sum(abs2, t_n)
    t_n_mag = sqrt(t_n_mag2)

    if t_n_mag < T(TN_MAG_THRESHOLD)
        c_n1 = one(T) - T(2 * n * n) * t_n_mag2

        n_term1 = T(2 * n)
        n_term3 = T(4 * n * n * n)

        c_n6 = (n_term1 + n_term3) / T(3)
        c_n2 = n_term1 - c_n6 * t_n_mag2
        c_n3 = T(2 * n * n) - T(4 * n * n + 2 * n * n * n * n) / T(3) * t_n_mag2
    else
        alpha_n = atan(t_n_mag)
        n_alpha_n = T(n) * alpha_n
        sin_n_alpha, cos_n_alpha = sincos(n_alpha_n)
        sin_2n_alpha = 2 * sin_n_alpha * cos_n_alpha
        cos_2n_alpha = cos_n_alpha * cos_n_alpha - sin_n_alpha * sin_n_alpha

        c_n1 = cos_2n_alpha
        c_n2 = sin_2n_alpha / t_n_mag
        c_n3 = 2 * sin_n_alpha * sin_n_alpha / t_n_mag2
        c_n6 = (T(2 * n) - c_n2) / t_n_mag2
    end

    return c_n1, c_n2, c_n3, c_n6
end

@muladd function update_velocity_multistep(v, r, dt, t, n::Int, ::Val{N}, param) where {N}
    T = eltype(v)
    q2m, Efunc, Bfunc = param.q2m, param.E, param.B

    E = SVector{3, T}(Efunc(r, t))
    B = SVector{3, T}(Bfunc(r, t))
    c = _speed_of_light(param)

    if c === nothing
        factor = T(q2m) * T(dt) / (2 * n)
        t_n = factor * B
        e_n = factor * E

        if N != 2
            f_N, e_corr_factor = _hyper_boris_factors(T, sum(abs2, t_n), Val{N}())
            e_dot_t = e_n ⋅ t_n
            e_n = f_N * e_n + (e_corr_factor * e_dot_t) * t_n
            t_n = f_N * t_n
        end

        c_n1, c_n2, c_n3, c_n6 = _multistep_rotation_coeffs(t_n, n)
        v_dot_t = v ⋅ t_n
        e_dot_t = e_n ⋅ t_n
        v_cross_t = v × t_n
        e_cross_t = e_n × t_n
        t_coeff = muladd(c_n3, v_dot_t, c_n6 * e_dot_t)

        return c_n1 * v + c_n2 * (v_cross_t + e_n) + c_n3 * e_cross_t + t_coeff * t_n
    end

    half_qdt2m = T(q2m) * (T(0.5) * T(dt))
    qE2m = half_qdt2m * E
    u⁻ = v + qE2m
    γ⁻ = sqrt(one(T) + sum(abs2, u⁻) / T(c)^2)
    t_n = (half_qdt2m / (n * γ⁻)) * B
    if N != 2
        f_N, _ = _hyper_boris_factors(T, sum(abs2, t_n), Val{N}())
        t_n = f_N * t_n
    end
    c_n1, c_n2, c_n3, _ = _multistep_rotation_coeffs(t_n, n)
    u⁺ = c_n1 * u⁻ + c_n2 * (u⁻ × t_n) + (c_n3 * (u⁻ ⋅ t_n)) * t_n
    return u⁺ + qE2m
end

"""
    update_velocity(v, r, dt, t, p, alg)

Advance the velocity `v` (or relativistic momentum-per-mass `γv`) by `dt`, evaluating the
fields at `(r, t)`.
"""
@inline @muladd function update_velocity(v, r, dt, t, p, ::Union{Boris, AdaptiveBoris})
    T = eltype(v)
    qdt2m = T(p.q2m) * (T(0.5) * T(dt))
    E = SVector{3, T}(p.E(r, t))
    B = SVector{3, T}(p.B(r, t))
    c = _speed_of_light(p)
    c === nothing && return update_velocity_boris(v, E, B, qdt2m)

    qE2m = qdt2m * E
    u⁻ = v + qE2m
    γ⁻ = sqrt(one(T) + sum(abs2, u⁻) / T(c)^2)
    return _boris_rotate(u⁻, (qdt2m / γ⁻) * B) + qE2m
end

@inline @muladd function update_velocity(
        v, r, dt, t, p, alg::Union{MultistepBoris{N}, AdaptiveMultistepBoris{N}}
    ) where {N}
    return update_velocity_multistep(v, r, dt, t, alg.n, Val{N}(), p)
end

@inline @muladd function update_velocity(v, r, dt, t, p, alg::Union{Vay, AdaptiveVay})
    c = _speed_of_light(p)
    c === nothing && throw(
        ArgumentError("$(nameof(typeof(alg))) requires speed of light `c` in `p`.")
    )
    T = eltype(v)
    qdt2m = T(p.q2m) * (T(0.5) * T(dt))
    E = SVector{3, T}(p.E(r, t))
    B = SVector{3, T}(p.B(r, t))
    return update_velocity_vay(v, E, B, qdt2m, inv(T(c)^2))
end

@inline @muladd function update_velocity(
        v, r, dt, t, p, alg::Union{HigueraCary, AdaptiveHigueraCary}
    )
    c = _speed_of_light(p)
    c === nothing && throw(
        ArgumentError("$(nameof(typeof(alg))) requires speed of light `c` in `p`.")
    )
    T = eltype(v)
    qdt2m = T(p.q2m) * (T(0.5) * T(dt))
    E = SVector{3, T}(p.E(r, t))
    B = SVector{3, T}(p.B(r, t))
    return update_velocity_higuera_cary(v, E, B, qdt2m, inv(T(c)^2))
end

# The velocity is carried at the half step, as in a leapfrog scheme: the cache
# holds `v(t - dt/2)` and the node velocity is reconstructed only for output.
# Changing `dt` re-centers the stored velocity onto the new half step, which is
# what keeps the scheme time-reversible under adaptive stepping.
#
# The fields are carried the same way. A step advances with the fields at the
# node it starts from and synchronises with the fields at the node it lands on,
# and the second node is where the next step starts, so the pair is handed over
# rather than evaluated again. That is what keeps a step at one field evaluation.
@inline @muladd function boris_initialize!(integrator, cache)
    t = integrator.t
    dt = integrator.dt
    p = integrator.p
    uprev = integrator.uprev
    r = SVector(uprev[1], uprev[2], uprev[3])
    v = SVector(uprev[4], uprev[5], uprev[6])

    fields = _fields_at(cache.fields, p, r, t)
    cache.v_half = update_velocity_half(v, r, dt, t, fields, integrator.alg)
    cache.fields = fields
    cache.dt_prev = dt

    integrator.kshortsize = 0
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    return
end

@inline function _needs_node_velocity(integrator)
    opts = integrator.opts
    opts.save_everystep && return true
    !isempty(opts.saveat) && return true
    !isempty(opts.callback.continuous_callbacks) && return true
    !isempty(opts.callback.discrete_callbacks) && return true
    return false
end

@inline @muladd function advance_boris!(integrator, cache)
    t = integrator.t
    dt = integrator.dt
    p = integrator.p
    alg = integrator.alg
    uprev = integrator.uprev
    r = SVector(uprev[1], uprev[2], uprev[3])

    fields = _fields_at(cache.fields, p, r, t)

    v_half_prev = if cache.dt_prev == dt
        cache.v_half
    else
        update_velocity_resync(cache.v_half, r, cache.dt_prev, dt, t, fields, alg)
    end

    r_new, v_half = advance_boris(v_half_prev, r, dt, t, fields, alg)

    if _needs_node_velocity(integrator)
        fields_new = _node_fields(p, r_new, t + dt)
        v_new = update_velocity_node(v_half, r_new, dt, t + dt, fields_new, alg)
        cache.fields = fields_new
    else
        v_new = v_half
        cache.fields = fields
    end

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

@muladd function perform_step!(integrator, cache::BorisConstantCache, repeat_step = false)
    r_new, v_new = advance_boris!(integrator, cache)
    integrator.u = vcat(r_new, v_new)
    return integrator.u
end

@muladd function perform_step!(integrator, cache::BorisCache, repeat_step = false)
    r_new, v_new = advance_boris!(integrator, cache)
    integrator.u[1] = r_new[1]
    integrator.u[2] = r_new[2]
    integrator.u[3] = r_new[3]
    integrator.u[4] = v_new[1]
    integrator.u[5] = v_new[2]
    integrator.u[6] = v_new[3]
    return
end

function postamble!(integrator::ODEIntegrator{<:AbstractBoris})
    if !_needs_node_velocity(integrator)
        cache = integrator.cache
        t = integrator.t
        r = SVector(integrator.u[1], integrator.u[2], integrator.u[3])
        fields = _fields_at(cache.fields, integrator.p, r, t)
        v_node = update_velocity_node(
            cache.v_half, r, integrator.dt, t, fields, integrator.alg
        )
        if integrator.u isa SVector
            integrator.u = vcat(r, v_node)
        else
            integrator.u[4] = v_node[1]
            integrator.u[5] = v_node[2]
            integrator.u[6] = v_node[3]
        end
        cache.fields = fields
    end
    _postamble!(integrator)
    return nothing
end

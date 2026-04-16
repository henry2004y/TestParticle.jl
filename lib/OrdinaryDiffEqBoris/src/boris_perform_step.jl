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

function OrdinaryDiffEqCore.initialize!(integrator, cache::BorisConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.kshortsize = 0
    return integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

function OrdinaryDiffEqCore.perform_step!(integrator, cache::BorisConstantCache, repeat_step = false)
    t = integrator.t
    dt = integrator.dt
    uprev = integrator.uprev
    f = integrator.f
    p = integrator.p

    # Extract particles from uprev
    # TestParticle.jl standard: u[1:3] is r, u[4:6] is v
    r = uprev[SVector(1, 2, 3)]
    v = uprev[SVector(4, 5, 6)]

    # Expected p structure in TestParticle.jl: (q2m, m, Efunc, Bfunc, Ffunc)
    q2m = p[1]
    Efunc = p[3]
    Bfunc = p[4]

    # Evaluate fields at time t
    # For a full time step, first we push velocity half step?
    # Wait, the traditional Boris pushes v by half step, r by full step, v by half step.
    # In TestParticle.jl, they do:
    # v_new = velocity_updater(v, r, dt, t + dt, p) # actually it uses t + 0.5*dt somewhere

    # We need to closely match what the Boris algorithm conceptually does.
    # Actually, a standard SciML step goes from t to t+dt.
    # If the user expects standard Boris, it is:
    r_half = r + v * (dt / 2)
    # Evaluate fields at t + dt/2, r_half
    t_half = t + dt / 2
    E = Efunc(r_half, t_half)
    B = Bfunc(r_half, t_half)

    qdt_2m = q2m * 0.5 * dt
    v_new = boris_velocity_update(v, E, B, qdt_2m)

    r_new = r_half + v_new * (dt / 2)

    u_new = vcat(r_new, v_new)
    return integrator.u = u_new
end

function OrdinaryDiffEqCore.initialize!(integrator, cache::BorisCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.kshortsize = 0
    return integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

function OrdinaryDiffEqCore.perform_step!(integrator, cache::BorisCache, repeat_step = false)
    # In-place version.
    # Since TestParticle.jl mostly focuses on SVector for small state vectors,
    # we can do a naive copy for the in-place cache.
    t = integrator.t
    dt = integrator.dt
    uprev = integrator.uprev
    f = integrator.f
    p = integrator.p

    r = SVector(uprev[1], uprev[2], uprev[3])
    v = SVector(uprev[4], uprev[5], uprev[6])

    q2m = p[1]
    Efunc = p[3]
    Bfunc = p[4]

    r_half = r + v * (dt / 2)
    t_half = t + dt / 2
    E = Efunc(r_half, t_half)
    B = Bfunc(r_half, t_half)

    qdt_2m = q2m * 0.5 * dt
    v_new = boris_velocity_update(v, E, B, qdt_2m)

    r_new = r_half + v_new * (dt / 2)

    integrator.u[1] = r_new[1]
    integrator.u[2] = r_new[2]
    integrator.u[3] = r_new[3]
    integrator.u[4] = v_new[1]
    integrator.u[5] = v_new[2]
    return integrator.u[6] = v_new[3]
end

@muladd function update_velocity_multistep(v, r, dt, t, n::Int, N::Int, param)
    q2m = param[1]
    Efunc = param[3]
    Bfunc = param[4]

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
    v_new = c_n1 * v +
        c_n2 * v_cross_t +
        c_n3 * v_dot_t * t_n +
        c_n4 * e_n +
        c_n5 * e_cross_t +
        c_n6 * e_dot_t * t_n

    return v_new
end

function OrdinaryDiffEqCore.initialize!(integrator, cache::MultistepBorisConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.kshortsize = 0
    return integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

function OrdinaryDiffEqCore.perform_step!(integrator, cache::MultistepBorisConstantCache, repeat_step = false)
    t = integrator.t
    dt = integrator.dt
    uprev = integrator.uprev
    p = integrator.p
    alg = integrator.alg

    r = uprev[SVector(1, 2, 3)]
    v = uprev[SVector(4, 5, 6)]

    # Half step update
    r_half = r + v * (dt / 2)
    t_half = t + dt / 2

    # Update velocity using multistep
    v_new = update_velocity_multistep(v, r_half, dt, t_half, alg.n, alg.N, p)

    r_new = r_half + v_new * (dt / 2)

    return integrator.u = vcat(r_new, v_new)
end

function OrdinaryDiffEqCore.initialize!(integrator, cache::MultistepBorisCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.kshortsize = 0
    return integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

function OrdinaryDiffEqCore.perform_step!(integrator, cache::MultistepBorisCache, repeat_step = false)
    t = integrator.t
    dt = integrator.dt
    uprev = integrator.uprev
    p = integrator.p
    alg = integrator.alg

    r = SVector(uprev[1], uprev[2], uprev[3])
    v = SVector(uprev[4], uprev[5], uprev[6])

    r_half = r + v * (dt / 2)
    t_half = t + dt / 2

    v_new = update_velocity_multistep(v, r_half, dt, t_half, alg.n, alg.N, p)

    r_new = r_half + v_new * (dt / 2)

    integrator.u[1] = r_new[1]
    integrator.u[2] = r_new[2]
    integrator.u[3] = r_new[3]
    integrator.u[4] = v_new[1]
    integrator.u[5] = v_new[2]
    return integrator.u[6] = v_new[3]
end

function OrdinaryDiffEqCore.initialize!(integrator, cache::AdaptiveBorisConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.kshortsize = 0
    return integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

function OrdinaryDiffEqCore.perform_step!(integrator, cache::AdaptiveBorisConstantCache, repeat_step = false)
    t = integrator.t
    dt = integrator.dt
    uprev = integrator.uprev
    p = integrator.p
    alg = integrator.alg

    r = uprev[SVector(1, 2, 3)]
    v = uprev[SVector(4, 5, 6)]

    q2m = p[1]
    Efunc = p[3]
    Bfunc = p[4]

    r_half = r + v * (dt / 2)
    t_half = t + dt / 2
    E = Efunc(r_half, t_half)
    B = Bfunc(r_half, t_half)

    qdt_2m = q2m * 0.5 * dt
    v_new = boris_velocity_update(v, E, B, qdt_2m)

    r_new = r_half + v_new * (dt / 2)
    integrator.u = vcat(r_new, v_new)

    # Adaptive Step proposition based on local gyroperiod
    Bmag = norm(Bfunc(r_new, t + dt))
    dt_new = (2π * alg.safety) / (abs(q2m) * Bmag)
    return set_proposed_dt!(integrator, dt_new)
end

function OrdinaryDiffEqCore.initialize!(integrator, cache::AdaptiveBorisCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.kshortsize = 0
    return integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

function OrdinaryDiffEqCore.perform_step!(integrator, cache::AdaptiveBorisCache, repeat_step = false)
    t = integrator.t
    dt = integrator.dt
    uprev = integrator.uprev
    p = integrator.p
    alg = integrator.alg

    r = SVector(uprev[1], uprev[2], uprev[3])
    v = SVector(uprev[4], uprev[5], uprev[6])

    q2m = p[1]
    Efunc = p[3]
    Bfunc = p[4]

    r_half = r + v * (dt / 2)
    t_half = t + dt / 2
    E = Efunc(r_half, t_half)
    B = Bfunc(r_half, t_half)

    qdt_2m = q2m * 0.5 * dt
    v_new = boris_velocity_update(v, E, B, qdt_2m)

    r_new = r_half + v_new * (dt / 2)

    integrator.u[1] = r_new[1]
    integrator.u[2] = r_new[2]
    integrator.u[3] = r_new[3]
    integrator.u[4] = v_new[1]
    integrator.u[5] = v_new[2]
    integrator.u[6] = v_new[3]

    Bmag = norm(Bfunc(r_new, t + dt))
    dt_new = (2π * alg.safety) / (abs(q2m) * Bmag)
    return set_proposed_dt!(integrator, dt_new)
end

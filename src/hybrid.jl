# Adaptive hybrid GC-FO solver

"""
    AdaptiveHybrid{T}

Hybrid algorithm that switches between adaptive Boris (FO) and RK45 (GC).
"""
struct AdaptiveHybrid{T}
    threshold::T
    dtmin::T
    dtmax::T
    safety_fo::T
    abstol::T
    reltol::T
    maxiters::Int
end

function AdaptiveHybrid(; threshold = 0.1, dtmax, dtmin = 1.0e-2 * dtmax, safety_fo = 0.1, abstol = 1.0e-6, reltol = 1.0e-6, maxiters = 10000)
    T = promote_type(typeof(threshold), typeof(dtmin), typeof(dtmax), typeof(safety_fo), typeof(abstol), typeof(reltol))
    return AdaptiveHybrid{T}(T(threshold), T(dtmin), T(dtmax), T(safety_fo), T(abstol), T(reltol), maxiters)
end

"""
    TraceHybridProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF}

Problem type for hybrid GC-FO tracing. Initial condition `u0` should be the full orbit state (6-element).
"""
struct TraceHybridProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <:
    AbstractODEProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P # (q2m, m, E, B, F)
    prob_func::PF
end

function TraceHybridProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    _f = ODEFunction{true, DEFAULT_SPECIALIZATION}(x -> nothing) # dummy func
    return TraceHybridProblem{
        typeof(u0), typeof(tspan), true, typeof(p), typeof(_f),
        typeof(prob_func),
    }(
        _f,
        u0,
        tspan,
        p,
        prob_func
    )
end

function solve(
        prob::TraceHybridProblem, alg::AdaptiveHybrid,
        ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1,
        isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true
    )
    return _solve(
        ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep
    )
end

function _solve(
        ::EnsembleSerial, prob, trajectories, alg::AdaptiveHybrid, savestepinterval,
        isoutofdomain, save_start, save_end, save_everystep
    )
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)))
    sols = Vector{sol_type}(undef, trajectories)
    irange = 1:trajectories

    _hybrid_adaptive!(
        sols, prob, irange, alg, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep
    )

    return sols
end

function _solve(
        ::EnsembleThreads, prob, trajectories, alg::AdaptiveHybrid, savestepinterval,
        isoutofdomain, save_start, save_end, save_everystep
    )
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)))
    sols = Vector{sol_type}(undef, trajectories)

    nchunks = Threads.nthreads()
    Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        _hybrid_adaptive!(
            sols, prob, irange, alg, savestepinterval, isoutofdomain,
            save_start, save_end, save_everystep
        )
    end

    return sols
end

function _get_sol_type(prob::TraceHybridProblem, dt)
    u0 = prob.u0
    tspan = prob.tspan
    T_t = typeof(tspan[1] + dt)
    t = Vector{T_t}(undef, 0)
    T = eltype(u0)
    u = Vector{SVector{6, T}}(undef, 0)
    interp = LinearInterpolation(t, u)
    alg = :hybrid

    sol = build_solution(prob, alg, t, u; interp = interp)
    return typeof(sol)
end

# Internal helpers to handle field calls with time
function _get_gc_parameters_at_t(xv, E, B, q, m, t)
    x, v = xv[SA[1:3...]], xv[SA[4:6...]]

    bparticle = B(x, t)
    Bmag_particle = norm(bparticle)
    b̂particle = bparticle / Bmag_particle
    # vector of Larmor radius
    ρ = (b̂particle × v) / (q / m * Bmag_particle)
    # Get the guiding center location
    X = x - ρ
    # Get EM field at guiding center
    b = B(X, t)
    Bmag = norm(b)
    b̂ = b / Bmag
    vpar = b̂ ⋅ v

    vperp = v - vpar * b̂
    e = E(X, t)
    vE = (e × b̂) / Bmag
    w = vperp - vE
    μ = m * (w ⋅ w) / (2 * Bmag)

    return X, vpar, μ
end

function _gc_to_full_at_t(state_gc, E_field, B_field, q, m, μ, t, phase = 0)
    R = get_x(state_gc)
    vpar = state_gc[4]

    E = E_field(R, t)
    B = B_field(R, t)
    Bmag = norm(B)
    b̂ = B / Bmag

    v_E = (E × b̂) / Bmag

    # perp speed, μ = m * w^2 / (2B)
    w_mag = sqrt(2 * μ * Bmag / m)

    # perp vector
    e1, e2 = get_perp_vector(b̂)
    v_gyr = w_mag * (cos(phase) * e1 + sin(phase) * e2)
    v_perp = v_gyr + v_E

    v = vpar * b̂ + v_perp

    # gyroradius
    Ω = q * Bmag / m
    ρ_vec = (b̂ × v_perp) / Ω

    x = R + ρ_vec

    return vcat(x, v)
end

function _hybrid_adaptive!(
        sols, prob, irange, alg, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep
    )
    (; tspan, p, u0) = prob
    # p = (q2m, m, Efunc, Bfunc, Ffunc)
    q2m, m, Efunc, Bfunc = p[1], p[2], p[3], p[4]
    q = q2m * m
    T = eltype(u0)
    paramBoris = BorisMethod(T)
    xv_fo = MVector{6, T}(undef)
    v_old = MVector{3, T}(undef)
    xv_save = MVector{6, T}(undef)

    # Determine if fields are time dependent
    is_td = is_time_dependent(Efunc) || is_time_dependent(Bfunc)

    # Safety parameters for RK45 (from gc_solver.jl)
    safety_gc = 0.9
    max_growth = 5.0
    min_growth = 0.2

    @fastmath @inbounds for i in irange
        traj = SVector{6, T}[]
        tsave = typeof(tspan[1])[]

        new_prob = prob.prob_func(prob, i, false)
        xv_fo .= new_prob.u0
        t = tspan[1]

        # Initial Mode Determination
        X_gc, vpar, μ = _get_gc_parameters_at_t(SVector{6, T}(xv_fo), Efunc, Bfunc, q, m, t)

        ϵ = get_adiabaticity(xv_fo[SA[1:3...]], Bfunc, q, m, μ, t)
        is_adiabatic = ϵ < alg.threshold

        if is_adiabatic
            mode = :GC
            xv_gc = SVector{4, T}(X_gc[1], X_gc[2], X_gc[3], vpar)
            p_gc = (q, q2m, μ, Efunc, Bfunc)

            B_vec = Bfunc(X_gc, t)
            Bmag = norm(B_vec)
            omega = abs(q2m * Bmag)
            dt = 0.5 * 2π / omega
        else
            mode = :FO
            # Initial dt for FO
            B_vec = Bfunc(xv_fo[SA[1:3...]], t)
            Bmag = norm(B_vec)
            omega = abs(q2m) * Bmag
            dt = alg.safety_fo / omega
            dt = clamp(dt, alg.dtmin, alg.dtmax)
            update_velocity!(xv_fo, paramBoris, (q2m, m, Efunc, Bfunc), -0.5 * dt, t)
        end

        if save_start
            push!(traj, SVector{6, T}(new_prob.u0))
            push!(tsave, t)
        end

        steps = 0
        it = 1
        while t < tspan[2] && steps < alg.maxiters
            if mode == :GC
                # Adiabaticity check
                ϵ = get_adiabaticity(xv_gc[SA[1:3...]], Bfunc, q, m, μ, t)
                if ϵ >= alg.threshold
                    # Switch to FO (GC -> FO)
                    # Switch to FO (GC -> FO)
                    @info "Switching from GC to FO at t = $t, ϵ = $ϵ"
                    mode = :FO
                    xv_fo_vec = _gc_to_full_at_t(xv_gc, Efunc, Bfunc, q, m, μ, t)
                    xv_fo .= xv_fo_vec

                    B_mag = norm(Bfunc(xv_fo[SA[1:3...]], t))
                    omega = abs(q2m * B_mag)
                    dt = alg.safety_fo / omega
                    dt = clamp(dt, alg.dtmin, alg.dtmax)
                    update_velocity!(xv_fo, paramBoris, (q2m, m, Efunc, Bfunc), -0.5 * dt, t)
                    continue
                end

                if t + dt > tspan[2]
                    dt = tspan[2] - t
                end

                dx, E_err = update_dp5(xv_gc, p_gc, dt, t)

                sum_sq_error = 0.0
                y_next = xv_gc + dx
                for k in 1:4
                    sc = alg.abstol + max(abs(xv_gc[k]), abs(y_next[k])) * alg.reltol
                    sum_sq_error += (E_err[k] / sc)^2
                end
                error_ratio = 0.5 * sqrt(sum_sq_error)

                if error_ratio <= 1.0
                    t += dt
                    xv_gc = y_next

                    if isoutofdomain(xv_gc, p_gc, t) && break end

                    if save_everystep && (it % savestepinterval == 0)
                        push!(traj, _gc_to_full_at_t(xv_gc, Efunc, Bfunc, q, m, μ, t))
                        push!(tsave, t)
                    end
                    steps += 1
                    it += 1
                end

                scale = error_ratio == 0.0 ? max_growth : safety_gc * (1.0 / error_ratio)^0.2
                scale = max(min_growth, min(scale, max_growth))
                dt *= scale
                dt < 1.0e-14 && break

            else # Mode == :FO
                t_sync = is_td ? t : zero(T)
                xv_save .= xv_fo
                update_velocity!(xv_save, paramBoris, (q2m, m, Efunc, Bfunc), 0.5 * dt, t_sync)

                _, _, μ = _get_gc_parameters_at_t(SVector{6, T}(xv_save), Efunc, Bfunc, q, m, t)

                ϵ = get_adiabaticity(xv_fo[SA[1:3...]], Bfunc, q, m, μ, t)
                if ϵ < alg.threshold
                    # Switch to GC
                    # Switch to GC
                    @info "Switching from FO to GC at t = $t, ϵ = $ϵ"
                    mode = :GC
                    xv_fo_sync = SVector{6, T}(xv_save)
                    X_gc, vpar, _ = _get_gc_parameters_at_t(xv_fo_sync, Efunc, Bfunc, q, m, t)
                    xv_gc = SVector{4, T}(X_gc[1], X_gc[2], X_gc[3], vpar)
                    p_gc = (q, q2m, μ, Efunc, Bfunc)

                    Bmag = norm(Bfunc(X_gc, t))
                    omega = abs(q2m * Bmag)
                    dt = 0.5 * 2π / omega
                    continue
                end

                if t + dt > tspan[2]
                    dt_step = tspan[2] - t
                    update_velocity!(xv_fo, paramBoris, (q2m, m, Efunc, Bfunc), 0.5 * dt, t_sync)
                    update_velocity!(xv_fo, paramBoris, (q2m, m, Efunc, Bfunc), -0.5 * dt_step, t_sync)
                    dt = dt_step
                end

                v_old .= @view xv_fo[4:6]

                if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
                    xv_save .= xv_fo
                    xv_save[4:6] .= v_old
                    update_velocity!(xv_save, paramBoris, (q2m, m, Efunc, Bfunc), 0.5 * dt, t_sync)
                    push!(traj, SVector{6, T}(xv_save))
                    push!(tsave, t)
                end

                t_mid = is_td ? t + 0.5 * dt : zero(T)
                update_velocity!(xv_fo, paramBoris, (q2m, m, Efunc, Bfunc), dt, t_mid)
                update_location!(xv_fo, dt)
                t += dt

                if isoutofdomain(xv_fo, (q2m, m, Efunc, Bfunc), t) && break end
                it += 1
                steps += 1

                # New dt for FO
                Bmag = norm(Bfunc(xv_fo[SA[1:3...]], t))
                if Bmag > 0
                    dt_new = alg.safety_fo / (abs(q2m) * Bmag)
                    dt_new = clamp(dt_new, alg.dtmin, alg.dtmax)
                else
                    dt_new = alg.dtmax
                end

                t_sync = is_td ? t : zero(T)
                update_velocity!(xv_fo, paramBoris, (q2m, m, Efunc, Bfunc), 0.5 * dt, t_sync)
                update_velocity!(xv_fo, paramBoris, (q2m, m, Efunc, Bfunc), -0.5 * dt_new, t_sync)
                dt = dt_new
            end
        end

        # Final Save
        if save_end && (isempty(tsave) || tsave[end] != t)
            if mode == :GC
                push!(traj, _gc_to_full_at_t(xv_gc, Efunc, Bfunc, q, m, μ, t))
            else
                t_final = is_td ? t : zero(T)
                update_velocity!(xv_fo, paramBoris, (q2m, m, Efunc, Bfunc), 0.5 * dt, t_final)
                push!(traj, SVector{6, T}(xv_fo))
            end
            push!(tsave, t)
        end

        sol_alg = :hybrid
        interp = LinearInterpolation(tsave, traj)
        retcode = steps >= alg.maxiters ? ReturnCode.MaxIters : ReturnCode.Success
        stats = nothing

        sols[i] = build_solution(
            prob, sol_alg, tsave, traj; interp = interp, retcode = retcode, stats = stats
        )
    end
    return
end

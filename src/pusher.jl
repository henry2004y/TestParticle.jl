# Native particle pusher

struct TraceProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <:
    AbstractODEProblem{uType, tType, isinplace}
    f::F
    "initial condition"
    u0::uType
    "time span"
    tspan::tType
    "(q2m, E, B)"
    p::P
    "function for setting initial conditions"
    prob_func::PF
end

get_EField(p::AbstractODEProblem) = get_EField(p.p)
get_BField(p::AbstractODEProblem) = get_BField(p.p)

get_BField(sol::AbstractODESolution) = get_BField(sol.prob)
get_EField(sol::AbstractODESolution) = get_EField(sol.prob)

DEFAULT_PROB_FUNC(prob, i, repeat) = prob

function TraceProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    _f = ODEFunction{true, DEFAULT_SPECIALIZATION}(x -> nothing) # dummy func
    return TraceProblem{
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
# For remake
function TraceProblem{iip}(; f, u0, tspan, p, prob_func) where {iip}
    return TraceProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f),
        typeof(prob_func),
    }(
        f,
        u0,
        tspan,
        p,
        prob_func
    )
end

struct BorisMethod{TV}
    # intermediate variables used in the solver
    v⁻::TV
    v′::TV
    v⁺::TV
    t_rotate::TV
    s_rotate::TV
    v⁻_cross_t::TV
    v′_cross_s::TV
end

function BorisMethod(T::Type{<:AbstractFloat} = Float64)
    v⁻ = MVector{3, T}(undef)
    v′ = MVector{3, T}(undef)
    v⁺ = MVector{3, T}(undef)
    t_rotate = MVector{3, T}(undef)
    s_rotate = MVector{3, T}(undef)
    v⁻_cross_t = MVector{3, T}(undef)
    v′_cross_s = MVector{3, T}(undef)

    return BorisMethod{typeof(v⁻)}(v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s)
end

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false

"""
    update_velocity(v, r, param, dt, t)

Update velocity using the Boris method, returning the new velocity as an SVector.
"""

function update_velocity(v, r, param, dt, t)
    q2m, _, Efunc, Bfunc, _ = param

    E = Efunc(r, t)
    B = Bfunc(r, t)

    t_rotate = q2m * B * 0.5 * dt
    t_mag2 = sum(abs2, t_rotate)
    s_rotate = 2 * t_rotate / (1 + t_mag2)

    v_minus = v + q2m * E * 0.5 * dt
    v_prime = v_minus + (v_minus × t_rotate)
    v_plus = v_minus + (v_prime × s_rotate)

    v_new = v_plus + q2m * E * 0.5 * dt

    return v_new
end


"""
    update_velocity!(xv, paramBoris, param, dt, t)

Update velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation.
Reference: [DTIC](https://apps.dtic.mil/sti/citations/ADA023511)
"""
@muladd function update_velocity!(xv, paramBoris, param, dt, t)
    (; v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s) = paramBoris
    q2m, _, Efunc, Bfunc, _ = param

    E = Efunc(xv, t)
    B = Bfunc(xv, t)
    # t vector

    for dim in 1:3
        t_rotate[dim] = q2m * B[dim] * 0.5 * dt
    end
    t_mag2 = sum(abs2, t_rotate)
    # s vector
    for dim in 1:3
        s_rotate[dim] = 2 * t_rotate[dim] / (1 + t_mag2)
    end
    # v-
    for dim in 1:3
        v⁻[dim] = xv[dim + 3] + q2m * E[dim] * 0.5 * dt
    end
    # v′
    cross!(v⁻, t_rotate, v⁻_cross_t)
    for dim in 1:3
        v′[dim] = v⁻[dim] + v⁻_cross_t[dim]
    end
    # v+
    cross!(v′, s_rotate, v′_cross_s)
    for dim in 1:3
        v⁺[dim] = v⁻[dim] + v′_cross_s[dim]
    end
    # v[n+1/2]
    for dim in 1:3
        xv[dim + 3] = v⁺[dim] + q2m * E[dim] * 0.5 * dt
    end

    return
end

"""
Update location in one timestep `dt`.
"""
@muladd function update_location!(xv, dt)
    xv[1] = xv[1] + xv[4] * dt
    xv[2] = xv[2] + xv[5] * dt
    xv[3] = xv[3] + xv[6] * dt

    return
end

"""
In-place cross product.
"""
@muladd function cross!(v1, v2, vout)
    vout[1] = v1[2] * v2[3] - v1[3] * v2[2]
    vout[2] = v1[3] * v2[1] - v1[1] * v2[3]
    vout[3] = v1[1] * v2[2] - v1[2] * v2[1]

    return
end

"""
    solve(prob::TraceProblem; trajectories::Int=1, dt::AbstractFloat,
    savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN,
        n::Int=1)

Trace particles using the Boris method with specified `prob`.

# keywords

  - `trajectories::Int`: number of trajectories to trace.
  - `dt::AbstractFloat`: time step.
  - `savestepinterval::Int`: saving output interval.
  - `isoutofdomain::Function`: a function with input of position and velocity vector `xv` that determines whether to stop tracing.
  - `n::Int`: number of substeps for the Multistep Boris method. Default is 1 (standard Boris).
  - `save_start::Bool`: save the initial condition. Default is `true`.
  - `save_end::Bool`: save the final condition. Default is `true`.
  - `save_everystep::Bool`: save the state at every `savestepinterval`. Default is `true`.
"""
function solve(
        prob::TraceProblem, ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1, dt::AbstractFloat,
        isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN, n::Int = 1,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        save_fields::Bool = false, save_work::Bool = false
    )
    return sols = _solve(
        ensemblealg, prob, trajectories, dt, savestepinterval, isoutofdomain, n,
        save_start, save_end, save_everystep, save_fields, save_work
    )
end


function _dispatch_boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
        save_start, save_end, save_everystep, save_fields, save_work
    )
    is_td = is_time_dependent(get_EField(prob)) || is_time_dependent(get_BField(prob))
    return if n == 1
        _boris!(
            sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
            save_start, save_end, save_everystep, Val(is_td), Val(save_fields), Val(save_work)
        )
    else
        # Multi-step boris does not support extra savings yet
        if save_fields || save_work
            @warn "Multi-step Boris does not support save_fields or save_work yet. Ignoring."
        end
        _multistep_boris!(
            sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
            save_start, save_end, save_everystep, Val(is_td)
        )
    end
end


function _solve(
        ::EnsembleSerial, prob, trajectories, dt, savestepinterval, isoutofdomain, n,
        save_start, save_end, save_everystep, save_fields, save_work
    )
    sols, nt,
        nout = _prepare(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, save_fields, save_work
    )
    irange = 1:trajectories
    _dispatch_boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
        save_start, save_end, save_everystep, save_fields, save_work
    )

    return sols
end


function _solve(
        ::EnsembleThreads, prob, trajectories, dt, savestepinterval, isoutofdomain, n,
        save_start, save_end, save_everystep, save_fields, save_work
    )
    sols, nt,
        nout = _prepare(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, save_fields, save_work
    )

    nchunks = Threads.nthreads()
    Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        _dispatch_boris!(
            sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
            save_start, save_end, save_everystep, save_fields, save_work
        )
    end

    return sols
end


function _get_sol_type(prob, dt, save_fields, save_work)
    u0 = prob.u0
    tspan = prob.tspan
    T_t = typeof(tspan[1] + dt)
    t = Vector{T_t}(undef, 0)
    T = eltype(u0)

    dim = 6
    if save_fields
        dim += 6
    end
    if save_work
        dim += 4
    end

    u = Vector{SVector{dim, T}}(undef, 0)
    interp = LinearInterpolation(t, u)
    alg = :boris

    sol = build_solution(prob, alg, t, u; interp = interp)
    return typeof(sol)
end


"""
Prepare for advancing.
"""
function _prepare(
        prob::TraceProblem, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, save_fields, save_work
    )

    ttotal = prob.tspan[2] - prob.tspan[1]
    nt = round(Int, ttotal / dt) |> abs

    nout = 0
    if save_start
        nout += 1
    end

    if save_everystep
        steps = nt ÷ savestepinterval
        last_is_step = (nt > 0) && (nt % savestepinterval == 0)
        nout += steps
        if !save_end && last_is_step
            nout -= 1
        end
        if save_end && !last_is_step
            nout += 1
        end
    elseif save_end
        nout += 1
    end

    sol_type = _get_sol_type(prob, dt, save_fields, save_work)
    sols = Vector{sol_type}(undef, trajectories)

    return sols, nt, nout
end

"""
Apply Boris method for particles with index in `irange`.
"""
function _boris!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
        save_start, save_end, save_everystep, ::Val{ITD}, ::Val{SAVE_FIELDS}, ::Val{SAVE_WORK}
    ) where {ITD, SAVE_FIELDS, SAVE_WORK}
    (; tspan, p, u0) = prob
    T = eltype(u0)

    q2m, m, Efunc, Bfunc, _ = p

    dim = 6
    if SAVE_FIELDS
        dim += 6
    end
    if SAVE_WORK
        dim += 4
    end

    @fastmath @inbounds for i in irange
        traj = Vector{SVector{dim, T}}(undef, nout)
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
            if !SAVE_FIELDS && !SAVE_WORK
                traj[iout] = u0_i
            else
                t0 = tspan[1]

                E = Efunc(r, t0)
                if SAVE_WORK
                    B_vec, Bmag, b_hat, ∇B, JB = TestParticle.get_B_parameters(r, t0, Bfunc)
                    # For work at t=0, everything is instantaneous.
                    v_par_vec = (v ⋅ b_hat) .* b_hat
                    Ω = q2m * Bmag
                    # P_par = q v_par . E
                    # P_gradB = q v_gradB . E
                    # P_curv = q v_curv . E
                    # P_ind = mu dB/dt

                    # We need q. q = q2m * m.
                    q = q2m * m
                    P_par = q * (v_par_vec ⋅ E)

                    # v_gradB = (m v_perp^2 / (2 q B^2)) * (b x ∇B)
                    # v_perp = v - v_par_vec
                    # Instead of manual calc, use existing drifts formula implicitly or re-derive
                    v_perp = v - v_par_vec
                    # m/q = 1/q2m
                    mu = m * norm(v_perp)^2 / (2 * Bmag) # magnetic moment
                    # gradB drift term: (mu/q) * (b x ∇B)
                    # (mu/q) = (m v_perp^2 / 2B) / q = (v_perp^2 / 2 * q2m * B) ?
                    # q/m * B = Omega.
                    # mu/q = (m/q) * v_perp^2 / 2B = v_perp^2 / (2 * q2m * B)
                    term_gradB = norm(v_perp)^2 / (2 * Ω)
                    v_gradB_vec = term_gradB * (b_hat × ∇B) / Bmag # b x ∇B / B
                    # Actually eq 238: norm(w)^2 * (b × ∇B) / (2 * Ω * Bmag)
                    # Matches.
                    P_gradB = q * (v_gradB_vec ⋅ E)

                    # v_curv
                    # eq 239: norm(v_par)^2 * (b × κ) / Ω
                    # κ = (JB * b + b * (-∇B ⋅ b)) / Bmag
                    κ = (JB * b_hat + b_hat * (-∇B ⋅ b_hat)) / Bmag
                    v_curv_vec = norm(v_par_vec)^2 * (b_hat × κ) / Ω
                    P_curv = q * (v_curv_vec ⋅ E)

                    # P_ind = mu * dB/dt
                    # dB/dt needs ForwardDiff on magnitude of B
                    dB_dt = ForwardDiff.derivative(t -> norm(Bfunc(r, t)), t0)
                    P_ind = mu * dB_dt

                    state_extra = SVector{4, T}(P_par, P_gradB, P_curv, P_ind)
                end

                if SAVE_FIELDS
                    B = Bfunc(r, t0)
                    state_fields = vcat(E, B)
                end

                if SAVE_FIELDS && SAVE_WORK

                    traj[iout] = vcat(r, v, state_fields, state_extra)
                elseif SAVE_FIELDS
                    traj[iout] = vcat(r, v, state_fields)
                elseif SAVE_WORK
                    traj[iout] = vcat(r, v, state_extra)
                end
            end
            tsave[iout] = tspan[1]
        end


        # push velocity back in time by 1/2 dt
        v = update_velocity(v, r, p, -0.5 * dt, tspan[1])

        it = 1
        while it <= nt
            v_prev = v
            t = ITD ? (it - 0.5) * dt : zero(dt)
            v = update_velocity(v, r, p, dt, t)

            if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
                iout += 1
                if iout <= nout
                    t_current = tspan[1] + (it - 1) * dt
                    # Approximate v_n from v_{n-1/2} (v_prev)
                    # Use v_prev because update_velocity logic is:
                    # v_{n+1/2} = update(v_{n-1/2}, ...)
                    # We want v_n = update(v_{n-1/2}, dt/2, ...)
                    v_save = update_velocity(
                        v_prev, r, p, 0.5 * dt,
                        ITD ? t_current : zero(dt)
                    )

                    if !SAVE_FIELDS && !SAVE_WORK
                        traj[iout] = vcat(r, v_save)
                    else
                        if SAVE_WORK

                            B_vec, Bmag_val, b_hat, ∇B, JB = TestParticle.get_B_parameters(r, t_current, Bfunc)
                            Ω = q2m * Bmag_val

                            E_val = Efunc(r, t_current)
                            v_par_vec = (v_save ⋅ b_hat) .* b_hat
                            q = q2m * m
                            P_par = q * (v_par_vec ⋅ E_val)

                            v_perp = v_save - v_par_vec
                            mu = m * norm(v_perp)^2 / (2 * Bmag_val)

                            term_gradB = norm(v_perp)^2 / (2 * Ω)
                            v_gradB_vec = term_gradB * (b_hat × ∇B) / Bmag_val
                            P_gradB = q * (v_gradB_vec ⋅ E_val)

                            κ = (JB * b_hat + b_hat * (-∇B ⋅ b_hat)) / Bmag_val
                            v_curv_vec = norm(v_par_vec)^2 * (b_hat × κ) / Ω
                            P_curv = q * (v_curv_vec ⋅ E_val)

                            dB_dt = ForwardDiff.derivative(t -> norm(Bfunc(r, t)), t_current)
                            P_ind = mu * dB_dt

                            state_extra = SVector{4, T}(P_par, P_gradB, P_curv, P_ind)
                        end
                        if SAVE_FIELDS
                            E_val = Efunc(r, t_current)
                            B_val = Bfunc(r, t_current)
                            state_fields = vcat(E_val, B_val)
                        end

                        if SAVE_FIELDS && SAVE_WORK
                            traj[iout] = vcat(r, v_save, state_fields, state_extra)
                        elseif SAVE_FIELDS
                            traj[iout] = vcat(r, v_save, state_fields)
                        elseif SAVE_WORK
                            traj[iout] = vcat(r, v_save, state_extra)
                        end
                    end
                    tsave[iout] = t_current
                end
            end


            r += v * dt
            if isoutofdomain(vcat(r, v), p, it * dt)
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
            v_final = update_velocity(v, r, p, dt_final, ITD ? t_final : zero(dt))

            if !SAVE_FIELDS && !SAVE_WORK
                traj[iout] = vcat(r, v_final)
            else
                # Calculate extra (reusing logic would be better but simple copy-paste for now to ensure capturing local vars)
                if SAVE_WORK
                    B_vec, Bmag_val, b_hat, ∇B, JB = TestParticle.get_B_parameters(r, t_final, Bfunc)
                    Ω = q2m * Bmag_val

                    E_val = Efunc(r, t_final)
                    v_par_vec = (v_final ⋅ b_hat) .* b_hat
                    q = q2m * m
                    P_par = q * (v_par_vec ⋅ E_val)

                    v_perp = v_final - v_par_vec
                    mu = m * norm(v_perp)^2 / (2 * Bmag_val)

                    term_gradB = norm(v_perp)^2 / (2 * Ω)
                    v_gradB_vec = term_gradB * (b_hat × ∇B) / Bmag_val
                    P_gradB = q * (v_gradB_vec ⋅ E_val)

                    κ = (JB * b_hat + b_hat * (-∇B ⋅ b_hat)) / Bmag_val
                    v_curv_vec = norm(v_par_vec)^2 * (b_hat × κ) / Ω
                    P_curv = q * (v_curv_vec ⋅ E_val)

                    dB_dt = ForwardDiff.derivative(t -> norm(Bfunc(r, t)), t_final)
                    P_ind = mu * dB_dt

                    state_extra = SVector{4, T}(P_par, P_gradB, P_curv, P_ind)
                end
                if SAVE_FIELDS
                    E_val = Efunc(r, t_final)
                    B_val = Bfunc(r, t_final)
                    state_fields = vcat(E_val, B_val)
                end

                if SAVE_FIELDS && SAVE_WORK
                    traj[iout] = vcat(r, v_final, state_fields, state_extra)
                elseif SAVE_FIELDS
                    traj[iout] = vcat(r, v_final, state_fields)
                elseif SAVE_WORK
                    traj[iout] = vcat(r, v_final, state_extra)
                end
            end
            tsave[iout] = t_final

        end

        if iout < nout
            resize!(traj, iout)
            resize!(tsave, iout)
        end
        traj_save = traj
        t = tsave

        alg = :boris
        interp = LinearInterpolation(t, traj_save)
        retcode = ReturnCode.Default
        stats = nothing

        sols[i] = build_solution(
            prob, alg, t, traj_save; interp = interp, retcode = retcode, stats = stats
        )
    end

    return
end

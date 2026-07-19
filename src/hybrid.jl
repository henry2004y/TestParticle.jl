# Adaptive hybrid GC-FO solver

const MIN_DT = 1.0e-14

"""
    AdaptiveHybrid{T}

Hybrid algorithm that switches between adaptive Boris (FO) and RK45 (GC).

Mode switching uses hysteresis with two adiabaticity thresholds to avoid
rapid GC ↔ FO chattering near a single bound:
- `threshold_gc_to_fo` (α): switch GC → FO when `ε ≥ α`.
- `threshold_fo_to_gc` (β): switch FO → GC when `ε < β`.
With `α ≥ β` the band `β ≤ ε < α` forms a buffer where the current mode is
retained, so a particle hovering near the boundary does not oscillate.

The adiabaticity parameter `ε` compared against the thresholds is selected by
the `adiabaticity` field, which may be:
- `:curvature` (default) → `ε_curv = ρ_L / R_c` (curvature drift; reproduces
  the legacy `get_adiabaticity` behavior exactly),
- `:gradB` → `ε_gradB = ρ_L / L_B`, with `L_B = |B| / |∇B|`,
- `:both` → OR of the two criteria: the solver switches to the full orbit
  whenever *either* `ε_curv ≥ α` *or* `ε_gradB ≥ α` (equivalently
  `max(ε_curv, ε_gradB) ≥ α`),
- `:jacobian` → `ε_jac = ρ_L · ‖JB‖_F / |B|`, CHIMP's original all-in-one
  criterion (Bates & Fuentes, 2020), where `JB` is the Jacobian of `B` and
  `‖JB‖_F` its Frobenius norm. It captures curvature, grad-B, torsion and
  shear at once, so it switches to full orbit in non-adiabatic regions that
  `:curvature`/`:gradB` miss (e.g. a field whose direction rotates in space
  while keeping constant magnitude).

Each check point records only the adiabaticity value used by the selected
`adiabaticity` mode (a single scalar per check) together with the active mode
in `sol.stats.adiabaticity`. This diagnostic is only saved when
`save_adiabaticity = true` (the default); set it to `false` to skip the
buffers and reduce memory/allocation overhead.
"""
struct AdaptiveHybrid{T}
    threshold_gc_to_fo::T
    threshold_fo_to_gc::T
    dtmin::T
    dtmax::T
    safety_fo::T
    abstol::T
    reltol::T
    maxiters::Int
    check_interval::Int
    adiabaticity::Symbol
    save_adiabaticity::Bool
end

function AdaptiveHybrid(;
        threshold = 0.1,
        threshold_gc_to_fo = threshold,
        threshold_fo_to_gc = threshold,
        dtmax, dtmin = 1.0e-2 * dtmax, safety_fo = 0.1,
        abstol = 1.0e-6, reltol = 1.0e-6, maxiters = 10000, check_interval = 10,
        adiabaticity = :curvature, save_adiabaticity = true
    )
    check_interval > 0 || throw(ArgumentError("check_interval must be positive."))
    threshold_gc_to_fo >= threshold_fo_to_gc ||
        throw(ArgumentError("threshold_gc_to_fo must be >= threshold_fo_to_gc for hysteresis."))
    adiabaticity in (:curvature, :gradB, :both, :jacobian) ||
        throw(
        ArgumentError(
            "adiabaticity must be one of :curvature, :gradB, :both, :jacobian."
        )
    )
    T = promote_type(
        typeof(threshold_gc_to_fo), typeof(threshold_fo_to_gc),
        typeof(dtmin), typeof(dtmax), typeof(safety_fo), typeof(abstol),
        typeof(reltol)
    )
    return AdaptiveHybrid{T}(
        T(threshold_gc_to_fo), T(threshold_fo_to_gc), T(dtmin), T(dtmax), T(safety_fo),
        T(abstol), T(reltol), maxiters, check_interval, adiabaticity, save_adiabaticity
    )
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

# For remake
function TraceHybridProblem{iip}(; f, u0, tspan, p, prob_func) where {iip}
    return TraceHybridProblem{
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


@inline function solve(
        prob::TraceHybridProblem, alg::AdaptiveHybrid,
        ensemblealg::EA = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true,
        save_everystep::Bool = true, verbose::Bool = false,
        seed::Union{Nothing, Integer} = nothing
    ) where {EA <: BasicEnsembleAlgorithm, F}
    return _solve(
        ensemblealg, prob, trajectories, alg,
        savestepinterval, isoutside,
        save_start, save_end, save_everystep, verbose, seed
    )
end

@inline function _solve(
        ::EnsembleSerial, prob::TraceHybridProblem,
        trajectories, alg::AdaptiveHybrid,
        savestepinterval, isoutside::F,
        save_start, save_end, save_everystep, verbose, seed
    ) where {F}
    sample_stats = _adia_sample_stats(prob, alg.save_adiabaticity)
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)), sample_stats)
    sols = Vector{sol_type}(undef, trajectories)
    irange = 1:trajectories

    elapsed_time = @elapsed _hybrid_adaptive!(
        sols, prob, irange, alg, savestepinterval,
        isoutside, save_start, save_end, save_everystep, verbose, seed
    )

    return EnsembleSolution(sols, elapsed_time, true)
end

@inline function _solve(
        ::EnsembleThreads, prob::TraceHybridProblem,
        trajectories, alg::AdaptiveHybrid,
        savestepinterval, isoutside::F,
        save_start, save_end, save_everystep, verbose, seed
    ) where {F}
    sample_stats = _adia_sample_stats(prob, alg.save_adiabaticity)
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)), sample_stats)
    sols = Vector{sol_type}(undef, trajectories)

    nchunks = Threads.nthreads()
    elapsed_time = @elapsed Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        _hybrid_adaptive!(
            sols, prob, irange, alg, savestepinterval,
            isoutside, save_start, save_end,
            save_everystep, verbose, seed
        )
    end

    return EnsembleSolution(sols, elapsed_time, true)
end

function _get_sol_type(prob::TraceHybridProblem, dt, stats = nothing)
    u0 = prob.u0
    tspan = prob.tspan
    T_t = typeof(tspan[1] + dt)
    t = Vector{T_t}(undef, 0)
    T = eltype(u0)
    u = Vector{SVector{6, T}}(undef, 0)
    interp = LinearInterpolation(t, u)
    alg = :hybrid

    sol = build_solution(prob, alg, t, u; interp = interp, stats = stats)
    return typeof(sol)
end

# Select the adiabaticity criterion used for the switching decision, compared
# against the shared thresholds. For `:both` we use OR logic: switch whenever
# *either* component exceeds the threshold.
@inline function _adia_select(comps, alg::AdaptiveHybrid)
    if alg.adiabaticity === :curvature
        return comps[1]
    elseif alg.adiabaticity === :gradB
        return comps[2]
    elseif alg.adiabaticity === :jacobian
        return comps[3]  # Frobenius-norm criterion
    else
        return max(comps[1], comps[2])  # :both → OR of the two criteria
    end
end

# Representative diagnostics structure used only to infer the solution type
# (the actual per-trajectory vectors are built inside the solver loop).
# Returns `nothing` when `save_adiabaticity` is false, so the affected
# `build_solution` is constructed without the `adiabaticity` stats field.
function _adia_sample_stats(prob::TraceHybridProblem, save_adiabaticity::Bool)
    save_adiabaticity || return nothing
    T = eltype(prob.u0)
    return (
        adiabaticity = (
            t = typeof(prob.tspan[1])[],
            components = T[],
            mode = Symbol[],
        ),
    )
end

# Fixed FO (Boris) time step for the hybrid solver, derived from the local field
# magnitude and the same `safety_fo` convention used by `AdaptiveBoris`.
@inline function _fo_dt(alg, q2m, Bfunc, r, t = 0.0)
    Bmag = norm(Bfunc(r, t))
    ω = abs(q2m) * Bmag
    dt = 2π * alg.safety_fo / ω
    return clamp(dt, alg.dtmin, alg.dtmax)
end

@inline function _hybrid_adaptive!(
        sols, prob::TraceHybridProblem, irange, alg,
        savestepinterval, isoutside::F,
        save_start, save_end, save_everystep, verbose, seed
    ) where {F}
    (; tspan, p, u0) = prob
    q2m, m, Efunc, Bfunc, _ = p
    q = q2m * m
    T = eltype(u0)

    is_td = is_time_dependent(Efunc) || is_time_dependent(Bfunc)

    # Safety parameters for RK45 (from gc_solver.jl)
    safety_gc = 0.9
    max_growth = 5.0
    min_growth = 0.2

    @inbounds for i in irange
        rng = isnothing(seed) ? default_rng() : Xoshiro(_rng_seed(seed, i))
        new_prob = prob.prob_func(prob, EnsembleContext(i, 1, 0, nothing, rng, seed))
        xv_fo = SVector{6, T}(new_prob.u0)
        r = xv_fo[SVector(1, 2, 3)]
        v = xv_fo[SVector(4, 5, 6)]
        t = tspan[1]

        # Initial Mode Determination
        X_gc, vpar, _, _, μ, phase = _get_gc_parameters(xv_fo, Efunc, Bfunc, q, m, t)

        comps0 = adiabaticity_components(r, Bfunc, q, m, μ, t)
        ϵ = _adia_select(comps0, alg)
        # Start in GC only when clearly adiabatic (ε < α); otherwise FO.
        is_adiabatic = ϵ < alg.threshold_gc_to_fo

        # Per-check diagnostic buffers, aligned with the decision points.
        # Allocated only when diagnostics are requested.
        if alg.save_adiabaticity
            adia_t = typeof(tspan[1])[]
            adia_vals = T[]
            adia_mode = Symbol[]
        end

        if is_adiabatic
            mode = :GC
            xv_gc = SVector{4, T}(X_gc[1], X_gc[2], X_gc[3], vpar)
            p_gc = (q, q2m, μ, Efunc, Bfunc)

            B_vec = Bfunc(X_gc, t)
            Bmag = norm(B_vec)
            ω = abs(q2m * Bmag)
            dt = 0.5 * 2π / ω
            verbose && @info "Initial mode: GC" ϵ t
        else
            mode = :FO
            dt = _fo_dt(alg, q2m, Bfunc, r, t)
            v = update_velocity(v, r, -0.5 * dt, t, p)
            verbose && @info "Initial mode: FO" ϵ t
        end
        if alg.save_adiabaticity
            push!(adia_t, t)
            push!(adia_vals, _adia_select(comps0, alg))
            push!(adia_mode, mode)
        end

        # Initialize solution containers with sizehint!
        # Estimate initial capacity based on timespan and initial dt
        estimated_steps = ceil(Int, (tspan[2] - tspan[1]) / dt)
        initial_capacity = min(max(10, estimated_steps + div(estimated_steps, 10)), 10000)
        traj = Vector{SVector{6, T}}(undef, 0)
        tsave = Vector{typeof(tspan[1])}(undef, 0)
        sizehint!(traj, initial_capacity)
        sizehint!(tsave, initial_capacity)

        if save_start
            push!(traj, SVector{6, T}(new_prob.u0))
            push!(tsave, t)
        end

        steps = 0
        it = 1
        while t < tspan[2] && steps < alg.maxiters
            if mode == :GC
                # Adiabaticity check
                if it % alg.check_interval == 0
                    comps = adiabaticity_components(
                        xv_gc[SVector(1, 2, 3)], Bfunc, q, m, μ, t
                    )
                    if alg.save_adiabaticity
                        push!(adia_t, t)
                        push!(adia_vals, _adia_select(comps, alg))
                        push!(adia_mode, :GC)
                    end
                    ϵ = _adia_select(comps, alg)
                    if ϵ >= alg.threshold_gc_to_fo
                        # Switch to FO (GC -> FO)
                        mode = :FO
                        verbose && @info "Switch GC → FO" ϵ t r = xv_gc[SVector(1, 2, 3)]
                        xv_fo_vec = _gc_to_full(
                            xv_gc, Efunc, Bfunc, q, m, μ, t, phase
                        )
                        xv_fo = xv_fo_vec
                        r = xv_fo[SVector(1, 2, 3)]
                        v = xv_fo[SVector(4, 5, 6)]

                        dt = _fo_dt(alg, q2m, Bfunc, r, t)
                        v = update_velocity(v, r, -0.5 * dt, t, p)
                        continue
                    end
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
                    y_next = xv_gc + dx

                    if isoutside(y_next, p_gc, t + dt)
                        break
                    end

                    t += dt
                    xv_gc = y_next

                    # Evolve phase using the (signed) gyrofrequency at the
                    # start-of-step guiding-center position.
                    Bmag_step = norm(Bfunc(get_x(xv_gc), t))
                    phase = mod2pi(phase - dt * (q2m * Bmag_step))

                    if save_everystep && (it % savestepinterval == 0)
                        push!(
                            traj,
                            _gc_to_full(xv_gc, Efunc, Bfunc, q, m, μ, t, phase)
                        )
                        push!(tsave, t)
                    end
                    steps += 1
                    it += 1
                end

                scale = error_ratio == 0.0 ?
                    max_growth : safety_gc * (1.0 / error_ratio)^0.2
                scale = max(min_growth, min(scale, max_growth))
                dt *= scale
                dt < MIN_DT && break

            else # Mode == :FO
                # Adiabaticity check only at the chosen cadence, using the
                # synchronized (position, velocity) pair at the integer time `t`.
                if it % alg.check_interval == 0
                    t_sync = is_td ? t : zero(T)
                    v_sync = update_velocity(v, r, 0.5 * dt, t_sync, p)
                    xv_sync = vcat(r, v_sync)
                    X_gc, vpar, _, _, μ_fo, phase_fo = _get_gc_parameters(
                        xv_sync, Efunc, Bfunc, q, m, t
                    )
                    comps = adiabaticity_components(X_gc, Bfunc, q, m, μ_fo, t)
                    if alg.save_adiabaticity
                        push!(adia_t, t)
                        push!(adia_vals, _adia_select(comps, alg))
                        push!(adia_mode, :FO)
                    end
                    ϵ = _adia_select(comps, alg)
                    if ϵ < alg.threshold_fo_to_gc
                        # Switch to GC (FO -> GC)
                        mode = :GC
                        verbose && @info "Switch FO → GC" ϵ t r = r
                        xv_gc = vcat(X_gc, vpar)
                        p_gc = (q, q2m, μ_fo, Efunc, Bfunc)
                        μ = μ_fo
                        phase = phase_fo

                        Bmag = norm(Bfunc(r, t))
                        omega = abs(q2m * Bmag)
                        dt = 0.5 * 2π / omega
                        continue
                    end
                end

                # Truncate the final step so we never overshoot the time span.
                if t + dt > tspan[2]
                    dt = tspan[2] - t
                end

                # One Boris leapfrog step with a fixed dt, identical in
                # convention to `Boris(dt)` so FO segments match it exactly.
                v_prev = v
                v = update_velocity(v, r, dt, t + 0.5 * dt, p)
                r_next = r + v * dt
                t_next = t + dt
                if isoutside(vcat(r_next, v), p, t_next)
                    break
                end

                if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
                    v_save = update_velocity(v_prev, r, 0.5 * dt, t, p)
                    push!(traj, vcat(r, v_save))
                    push!(tsave, t)
                end

                r = r_next
                t = t_next
                it += 1
                steps += 1
            end
        end

        # Final Save
        if save_end && (isempty(tsave) || tsave[end] != t)
            if mode == :GC
                push!(traj, _gc_to_full(xv_gc, Efunc, Bfunc, q, m, μ, t, phase))
            else
                t_final = is_td ? t : zero(T)
                v_final = update_velocity(v, r, 0.5 * dt, t_final, p)
                push!(traj, vcat(r, v_final))
            end
            push!(tsave, t)
        end

        sol_alg = :hybrid
        interp = LinearInterpolation(tsave, traj)
        retcode = steps >= alg.maxiters ? ReturnCode.MaxIters : ReturnCode.Success
        stats = alg.save_adiabaticity ?
            (adiabaticity = (t = adia_t, components = adia_vals, mode = adia_mode),) :
            nothing

        sols[i] = build_solution(
            prob, sol_alg, tsave, traj; interp = interp, retcode = retcode, stats = stats
        )
    end
    return
end

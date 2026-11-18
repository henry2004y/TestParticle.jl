module TestParticleBorisPushersExt

using TestParticle
import TestParticle: TraceProblem, FullOrbit, RelativisticOrbit, AdaptiveHybrid,
    _equations_model, _speed_of_light, _lorentz_factor, _model,
    SavingPlan, use_saveat, saveat_reached, saveat_interpolate, _span_direction,
    Field, ZeroField, FieldInterpolator, FieldInterpolator1D,
    FieldInterpolator2D, SphericalFieldInterpolator,
    _to_gpu_grid,
    get_x, get_magnetic_properties, _promote_p,
    _gc_to_full, _gc_state_at_switch, _gc_switch_state_checked,
    init_gc_integrator, reinit_gc_integrator!, gc_step!,
    is_time_dependent, _foreach_chunk, ODE_DEFAULT_ISOUTOFDOMAIN, DEFAULT_PROB_FUNC

using BorisPushers
import BorisPushers: AbstractBoris, Boris, AdaptiveBoris, MultistepBoris,
    AdaptiveMultistepBoris, Vay, AdaptiveVay, HigueraCary, AdaptiveHigueraCary,
    CachedFields, update_velocity, adapt_field_to_gpu, _node_fields

using SciMLBase: AbstractODEProblem, BasicEnsembleAlgorithm,
    DiscreteCallback, EnsembleProblem, EnsembleSolution, EnsembleContext,
    FullSpecialize, LinearInterpolation, ODEProblem, ReturnCode,
    build_solution, isadaptive, terminate!, EnsembleSerial, EnsembleThreads
import SciMLBase: solve

using LinearAlgebra: norm, ⋅, ×
using MuladdMacro: @muladd
using StaticArrays: SVector, SA
import Adapt
using BorisPushers: Backend, CPU
using Random: default_rng, Xoshiro
using PrecompileTools: @setup_workload, @compile_workload

# --- Boris ODE Problem Adaptation ---

const FOProblem = TraceProblem{FullOrbit}
const RelFOProblem = TraceProblem{RelativisticOrbit}
const BorisProblem = Union{FOProblem, RelFOProblem}
const RelOnlyBoris = Union{Vay, AdaptiveVay, HigueraCary, AdaptiveHigueraCary}

_boris_rhs(u, p, t) = nothing
_boris_rel_rhs(u, p, t) = nothing
_equations_model(::typeof(_boris_rhs)) = FullOrbit
_equations_model(::typeof(_boris_rel_rhs)) = RelativisticOrbit

# Boris solvers build their own update from `p.q2m`, `p.E` and `p.B` (plus `c` for
# relativistic orbits); `AdaptiveHybrid` only supports `FullOrbit`.
function _check_boris_problem(prob, alg::AbstractBoris = Boris())
    M = _equations_model(prob.f)
    M in (FullOrbit, RelativisticOrbit) && prob.p.F isa ZeroField || throw(
        ArgumentError(
            "Boris solvers need the default `FullOrbit` or `RelativisticOrbit` " *
                "equations without `F`."
        )
    )
    if alg isa RelOnlyBoris && M !== RelativisticOrbit
        throw(
            ArgumentError(
                "$(nameof(typeof(alg))) requires `TraceProblem{RelativisticOrbit}`."
            )
        )
    end
    return
end

function _check_hybrid_problem(prob)
    _equations_model(prob.f) === FullOrbit && prob.p.F isa ZeroField || throw(
        ArgumentError(
            "AdaptiveHybrid needs the default `FullOrbit` equations without `F`."
        )
    )
    return
end

@inline function _boris_promote_p(::FOProblem, ::Type{T}, p) where {T}
    p_T = _promote_p(T, p)
    p_T isa NamedTuple || return (; q2m = T(p.q2m), E = p.E, B = p.B)
    return hasproperty(p_T, :c) ? (; p_T..., c = nothing) : p_T
end

@inline function _boris_promote_p(::RelFOProblem, ::Type{T}, p) where {T}
    p_T = _promote_p(T, p)
    c_val = T(_speed_of_light(p))
    return p_T isa NamedTuple ? (; p_T..., c = c_val) :
        (; q2m = T(p.q2m), E = p.E, B = p.B, c = c_val)
end

function SciMLBase.__solve(
        eprob::EnsembleProblem{<:TraceProblem}, alg::AbstractBoris, ens::EnsembleKernel;
        kwargs...
    )
    base_prob = eprob.prob
    _check_boris_problem(base_prob, alg)
    base_model = _equations_model(base_prob.f)
    T = eltype(base_prob.u0)
    adapted_base = TraceProblem{base_model}(
        base_prob.f, base_prob.u0, base_prob.tspan,
        _boris_promote_p(base_prob, T, base_prob.p);
        prob_func = base_prob.prob_func
    )
    prob_func = eprob.prob_func
    wrapped_prob_func = if prob_func === DEFAULT_PROB_FUNC
        DEFAULT_PROB_FUNC
    else
        (prob, ctx) -> (
            new_prob = prob_func(prob, ctx);
            _check_boris_problem(new_prob, alg);
            _equations_model(new_prob.f) === base_model || throw(
                ArgumentError("prob_func cannot change the orbit model in EnsembleKernel.")
            );
            new_prob
        )
    end
    eprob_adapted = EnsembleProblem(
        adapted_base;
        prob_func = wrapped_prob_func,
        output_func = eprob.output_func,
        reduction = eprob.reduction,
        u_init = eprob.u_init,
        safetycopy = eprob.safetycopy
    )
    return invoke(
        SciMLBase.__solve,
        Tuple{SciMLBase.AbstractEnsembleProblem, AbstractBoris, EnsembleKernel},
        eprob_adapted, alg, ens; kwargs...
    )
end

function _boris_problem(prob::BorisProblem)
    T = eltype(prob.u0)
    tspan = (T(prob.tspan[1]), T(prob.tspan[2]))
    p = _boris_promote_p(prob, T, prob.p)
    rhs = prob isa RelFOProblem ? _boris_rel_rhs : _boris_rhs
    return ODEProblem{false, FullSpecialize}(
        rhs, SVector{6, T}(prob.u0), tspan, p
    )
end

function _boris_initial_dt(prob::BorisProblem, alg)
    hasproperty(alg, :safety) ||
        throw(ArgumentError("$(nameof(typeof(alg))) needs a time step `dt`."))
    q2m = prob.p.q2m
    u0 = prob.u0
    r0 = SVector{3}(u0[1], u0[2], u0[3])
    v0 = SVector{3}(u0[4], u0[5], u0[6])
    Bmag = norm(prob.p.B(r0, prob.tspan[1]))
    tdir = sign(prob.tspan[2] - prob.tspan[1])
    γ = _lorentz_factor(v0, _model(prob), prob.p)

    return tdir * (2π * alg.safety * γ) / (abs(q2m) * Bmag)
end

function _boris_check_limits(prob::BorisProblem, dt, alg, maxiters)
    T = typeof(dt)
    tspan = prob.tspan
    timescale = max(abs(tspan[1]), abs(tspan[2]), abs(tspan[2] - tspan[1]))
    min_dt = 10 * eps(T) * timescale
    if abs(dt) < min_dt
        throw(
            ArgumentError(
                "time step dt is too small, violating min_dt = 10 * eps(typeof(dt))"
            )
        )
    end

    if !isadaptive(alg)
        ttotal = prob.tspan[2] - prob.tspan[1]
        nt = abs(round(Int, ttotal / dt))
        if nt > maxiters
            throw(
                ArgumentError("number of iterations nt ($nt) exceeds maxiters ($maxiters)")
            )
        end
    end

    return
end

function _boris_callback(isoutside::F) where {F}
    isoutside === ODE_DEFAULT_ISOUTOFDOMAIN && return nothing

    condition = (u, t, integrator) -> isoutside(u, integrator.p, t)
    affect! = function (integrator)
        integrator.u = integrator.uprev
        return terminate!(integrator)
    end

    return DiscreteCallback(condition, affect!)
end

function SciMLBase.solve(
        prob::TraceProblem, alg::AbstractBoris;
        saveat = (),
        dt::Union{Nothing, Real} = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true,
        save_end::Bool = true,
        save_everystep::Bool = isempty(saveat),
        maxiters::Int = 1_000_000,
        rng = nothing,
        seed = nothing,
        kwargs...
    ) where {F}
    _check_boris_problem(prob, alg)
    T = eltype(prob.u0)
    step = dt === nothing ? T(_boris_initial_dt(prob, alg)) : T(dt)
    _boris_check_limits(prob, step, alg, maxiters)

    ode_prob = _boris_problem(prob)
    callback = _boris_callback(isoutside)

    return if isempty(saveat)
        solve(
            ode_prob, alg; dt = step, save_start, save_end, save_everystep,
            maxiters, callback, dense = false, kwargs...
        )
    else
        solve(
            ode_prob, alg; dt = step, saveat, save_start, save_end,
            maxiters, callback, dense = false, kwargs...
        )
    end
end

# --- AdaptiveHybrid Solver Implementation ---

@inline function SciMLBase.solve(
        prob::TraceProblem, alg::AdaptiveHybrid,
        ensemblealg::EA = EnsembleSerial();
        trajectories::Int = 1,
        saveat = (),
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true,
        save_everystep::Bool = true, verbose::Bool = false,
        seed::Union{Nothing, Integer} = nothing
    ) where {EA <: BasicEnsembleAlgorithm, F}
    _check_hybrid_problem(prob)
    length(prob.u0) == 6 ||
        throw(ArgumentError("u0 must have 6 elements (full orbit state)."))
    prob.tspan[2] != prob.tspan[1] ||
        throw(ArgumentError("tspan[2] must not equal tspan[1]."))
    plan = SavingPlan(
        saveat, prob.tspan, _span_direction(prob.tspan), typeof(prob.tspan[1])
    )
    args = (
        ensemblealg, prob, trajectories, alg,
        plan, isoutside, save_start, save_end, save_everystep, verbose, seed,
    )
    return alg.save_adiabaticity ? _solve(args..., Val(true)) : _solve(args..., Val(false))
end

@inline function _solve(
        ensemblealg::Union{EnsembleSerial, EnsembleThreads}, prob::FOProblem,
        trajectories, alg::AdaptiveHybrid,
        plan, isoutside::F,
        save_start, save_end, save_everystep, verbose, seed, v_save::Val{SaveAdia}
    ) where {F, SaveAdia}
    # The seeds of SciMLBase's ensemble solve, so every solver draws the same states.
    sim_seeds = (isnothing(seed) || prob.prob_func === DEFAULT_PROB_FUNC) ? nothing :
        SciMLBase.generate_sim_seeds(nothing, seed, trajectories)
    if trajectories == 1
        t0 = time_ns()
        sol_1 = _hybrid_one(
            prob, 1, alg, plan, isoutside, save_start, save_end, save_everystep,
            verbose, sim_seeds, nothing, v_save, Val(false)
        )
        elapsed_time = (time_ns() - t0) * 1.0e-9
        return EnsembleSolution([sol_1], elapsed_time, true)
    end
    sample_stats = _adia_sample_stats(prob, SaveAdia)
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)), sample_stats)
    sols = Vector{sol_type}(undef, trajectories)

    elapsed_time = @elapsed _foreach_chunk(ensemblealg, trajectories) do irange
        _hybrid_adaptive!(
            sols, prob, irange, alg, plan, isoutside,
            save_start, save_end, save_everystep, verbose, sim_seeds, v_save
        )
    end

    return EnsembleSolution(sols, elapsed_time, true)
end

function _get_sol_type(prob::FOProblem, dt, stats = nothing)
    u0 = prob.u0
    tspan = prob.tspan
    T_t = typeof(tspan[1] + dt)
    t = Vector{T_t}(undef, 0)
    T = eltype(u0)
    u = Vector{SVector{6, T}}(undef, 0)
    interp = LinearInterpolation(t, u)
    alg = :hybrid

    sol = build_solution(prob, alg, t, u; interp, stats)
    return typeof(sol)
end

@inline function _adia_value(mp, q, m, vpar, μ, alg::AdaptiveHybrid)
    (; ∇B, κ, Bmag, JB) = mp
    iszero(Bmag) && return oftype(Bmag, Inf)
    inv_Bmag = inv(Bmag)
    q2m = q / m
    ρ = hypot(vpar, sqrt(2 * (μ / m) * Bmag)) * inv_Bmag / abs(q2m)
    if alg.adiabaticity === :curvature
        return ρ * norm(κ)
    elseif alg.adiabaticity === :gradB
        return ρ * norm(∇B) * inv_Bmag
    elseif alg.adiabaticity === :jacobian
        return ρ * sqrt(sum(abs2, JB)) * inv_Bmag
    else
        return ρ * max(norm(κ), norm(∇B) * inv_Bmag)
    end
end

@muladd @inline function _frame_vpar_μ(v, E, b̂, Bmag, m)
    inv_Bmag = inv(Bmag)
    vpar = b̂ ⋅ v
    w = v - vpar * b̂ - (E × b̂) * inv_Bmag
    return vpar, m * (w ⋅ w) * (inv_Bmag / 2)
end

# Check both the particle and its candidate guiding center: near a field null
# `r - ρ(r)` can land far away in an adiabatic field, while on a symmetry axis
# `κ(r)` can vanish even though `X_gc` sits in strong curvature.
@inline function _fo_adia_from_mp(
        r::SVector{3, T}, v, E_r, mp_r, Efunc, Bfunc, q, m, t, alg::AdaptiveHybrid,
        short_circuit::Bool
    ) where {T}
    iszero(mp_r.Bmag) && return oftype(mp_r.Bmag, Inf)
    vpar_r, μ_r = _frame_vpar_μ(v, E_r, mp_r.b̂, mp_r.Bmag, m)
    ϵ_r = _adia_value(mp_r, q, m, vpar_r, μ_r, alg)
    short_circuit && ϵ_r >= T(alg.threshold_fo_to_gc) && return ϵ_r

    q2m = q / m
    X_gc = r - (mp_r.B × v) / (q2m * (mp_r.B ⋅ mp_r.B))
    mp_gc = get_magnetic_properties(X_gc, t, Bfunc)
    iszero(mp_gc.Bmag) && return oftype(mp_gc.Bmag, Inf)
    E_gc = SVector{3, T}(Efunc(X_gc, t))
    vpar_gc, μ_gc = _frame_vpar_μ(v, E_gc, mp_gc.b̂, mp_gc.Bmag, m)
    ϵ_gc = _adia_value(mp_gc, q, m, vpar_gc, μ_gc, alg)
    return max(ϵ_r, ϵ_gc)
end

@inline function _fo_adia_value(
        r::SVector{3, T}, v, E_r, Efunc, Bfunc, q, m, t, alg::AdaptiveHybrid
    ) where {T}
    mp_r = get_magnetic_properties(r, t, Bfunc)
    ϵ = _fo_adia_from_mp(r, v, E_r, mp_r, Efunc, Bfunc, q, m, t, alg, false)
    return ϵ, mp_r.Bmag
end

function _adia_sample_stats(prob::FOProblem, save_adiabaticity::Bool)
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

@inline function _fo_dt(alg, q2m, Bmag)
    ω = abs(q2m) * Bmag
    dt = 2π * alg.safety_fo / ω # Inf for a vanishing field, hence clipped to dtmax
    return clamp(dt, alg.dtmin, alg.dtmax)
end

# Compare against the absolute deadline also used to clip steps. Subtracting
# `t - t_last_check` can round just below `check_dt` right after a step landed on
# the deadline, which would skip the check and leave the following step unclipped.
@inline function _check_due(it, t, t_last_check, check_interval, deadline, tdir)
    scale = max(abs(t), abs(t_last_check), abs(deadline))
    tol = scale isa AbstractFloat ? 4 * eps(scale) : zero(scale)
    return tdir * (t - t_last_check) > zero(t) &&
        (it % check_interval == 0 || tdir * (t - deadline) + tol >= zero(t))
end

"""
Multiplier on the check interval while the solver is in the full-orbit mode: `1`
next to the exit threshold `β`, growing with `ϵ/β` up to `check_margin_max`. `ϵ`
changes by `O(ϵ)` over one curvature radius, so an `ϵ` that is a factor `k` above
`β` needs of order `k` of those excursions before the mode can flip; checking in
between only costs time. The guiding-center side is deliberately left at the
nominal interval: shortening it there would delay detection of a transition
*into* a non-adiabatic region, which is the case the check exists for.
"""
@inline function _check_scale(ϵ, alg, mode)
    mode === :FO || return 1.0
    β = alg.threshold_fo_to_gc
    iszero(β) && return 1.0
    return clamp(ϵ / β - 1, 1.0, alg.check_margin_max)
end

@inline function _hybrid_save!(
        traj, tsave, plan, isave, t_last, y_last, t_out, y_out
    )
    if use_saveat(plan)
        nsave = length(plan.times)
        while isave <= nsave && saveat_reached(plan.times[isave], t_out, plan.dir)
            t_target = plan.times[isave]
            push!(traj, saveat_interpolate(t_last, y_last, t_out, y_out, t_target))
            push!(tsave, t_target)
            isave += 1
        end
    else
        push!(traj, y_out)
        push!(tsave, t_out)
    end

    return isave, t_out, y_out
end

# Out of line: inlined `@info` expansions in the hot loop slow it down whenever
# `verbose` is not a compile-time constant.
@noinline _log_mode(msg, ϵ, t) = @info msg ϵ t
@noinline _log_mode(msg, ϵ, t, r) = @info msg ϵ t r

@inline _record_adia!(::Val{false}, adia_t, adia_vals, adia_mode, tc, ec, mc) = nothing
@inline function _record_adia!(::Val{true}, adia_t, adia_vals, adia_mode, tc, ec, mc)
    if !isempty(adia_t) && adia_t[end] == tc
        adia_vals[end] = ec
        adia_mode[end] = mc
    else
        push!(adia_t, tc)
        push!(adia_vals, ec)
        push!(adia_mode, mc)
    end
    return nothing
end

# Tighten the guiding-center step check cadence down to every step once ε exceeds
# a quarter of `threshold_gc_to_fo`, leaving deeply adiabatic runs at the nominal
# `check_interval`.
@inline function _gc_check_interval(ϵ::T, alg::AdaptiveHybrid) where {T}
    α = T(alg.threshold_gc_to_fo)
    iszero(α) && return 1
    r = ϵ / α
    r <= T(1 // 4) && return alg.check_interval
    scale = clamp((T(2 // 5) - r) / T(3 // 20), zero(T), one(T))
    return max(1, round(Int, alg.check_interval * scale))
end

@inline function _advance_gc_state(
        t_next, dt_step, y_next::SVector{4, T}, B_next, Bmag_next::T,
        phase::T, Bmag_curr::T, q2m, q, m, μ, Efunc, need_full::Bool
    ) where {T}
    phase_next = T(mod2pi(phase - dt_step * (q2m * (Bmag_curr + Bmag_next)) / 2))
    xv_full = need_full ? _gc_to_full(
            y_next, Efunc, B_next, Bmag_next, q, m, μ, t_next, phase_next
        ) : zero(SVector{6, T})
    return phase_next, xv_full
end

@inline function _enter_fo_from_gc(
        xv_fo, p_fo, q2m, t, dt_gc, tdir, alg, clip_dt::Bool
    )
    r = xv_fo[SVector(1, 2, 3)]
    v = xv_fo[SVector(4, 5, 6)]
    fields = _node_fields(p_fo, r, t)
    Tt = typeof(t)
    dt_lim = Tt(_fo_dt(alg, q2m, norm(fields.B(r, t))))
    dt = clip_dt ? tdir * Tt(clamp(abs(dt_gc), alg.dtmin, dt_lim)) : tdir * dt_lim
    v = update_velocity_half(v, r, dt, t, fields, Boris())
    return r, v, dt
end

struct _GCCheckpoint{Tt, T}
    t::Tt
    xv_gc::SVector{4, T}
    phase::T
    Bmag::T
    ϵ::T
    dt::Tt
    steps::Int
    isave::Int
    t_last::Tt
    y_last::SVector{6, T}
    ntraj::Int
end

@inline function _restore_gc_checkpoint!(gc_chk::_GCCheckpoint, traj, tsave)
    resize!(traj, gc_chk.ntraj)
    resize!(tsave, gc_chk.ntraj)
    return (
        gc_chk.t, gc_chk.xv_gc, gc_chk.phase, gc_chk.Bmag, gc_chk.ϵ,
        gc_chk.dt, gc_chk.steps, gc_chk.isave, gc_chk.t_last, gc_chk.y_last,
    )
end

@inline function _bisect_gc_to_fo!(
        integrator_gc, t::Tt, t_hi::Tt, xv_gc::SVector{4, T}, phase::T, Bmag_curr::T,
        ϵ::T, p_gc, q2m, q, m, μ::T, Efunc, Bfunc, alg::AdaptiveHybrid,
        tspan_end::Tt, isoutside::F, need_save_gc::Bool,
        traj, tsave, plan, isave::Int, t_last::Tt, y_last::SVector{6, T},
        v_save::Val{SaveAdia}, adia_t, adia_vals, adia_mode,
        steps::Int, overshoot_limit::T
    ) where {Tt, T, F, SaveAdia}
    switched_at_boundary = false
    xv_fo = zero(SVector{6, T})
    for _ in 1:16
        dt_fo_loc = Tt(_fo_dt(alg, q2m, Bmag_curr))
        abs(t_hi - t) <= dt_fo_loc / 2 && break
        dt_mid = (t_hi - t) / 2
        reinit_gc_integrator!(integrator_gc, xv_gc, (t, tspan_end), p_gc, dt_mid)
        code_m = gc_step!(integrator_gc, t + dt_mid, isoutside, p_gc, tspan_end)
        if code_m == ReturnCode.Terminated
            return (
                t, t_hi, xv_gc, phase, Bmag_curr, ϵ, xv_fo,
                isave, t_last, y_last, steps, false, code_m,
            )
        end
        code_m != ReturnCode.Success && break
        t_m = Tt(integrator_gc.t)
        y_m = SVector{4, T}(integrator_gc.u)
        mp_m = get_magnetic_properties(get_x(y_m), t_m, Bfunc)
        ϵ_m = T(_adia_value(mp_m, q, m, y_m[4], μ, alg))
        if ϵ_m < alg.threshold_gc_to_fo
            Bmag_m = T(mp_m.Bmag)
            phase, y_out = _advance_gc_state(
                t_m, t_m - t, y_m, mp_m.B, Bmag_m,
                phase, Bmag_curr, q2m, q, m, μ, Efunc, need_save_gc
            )
            Bmag_curr = Bmag_m
            t = t_m
            xv_gc = y_m
            ϵ = ϵ_m
            if need_save_gc
                isave, t_last, y_last = _hybrid_save!(
                    traj, tsave, plan, isave, t_last, y_last, t, y_out
                )
            end
            _record_adia!(v_save, adia_t, adia_vals, adia_mode, t, ϵ, :GC)
            steps += 1
            steps >= alg.maxiters && break
        elseif ϵ_m <= overshoot_limit
            Bmag_m = T(mp_m.Bmag)
            phase, xv_fo = _advance_gc_state(
                t_m, t_m - t, y_m, mp_m.B, Bmag_m,
                phase, Bmag_curr, q2m, q, m, μ, Efunc, true
            )
            Bmag_curr = Bmag_m
            t = t_m
            xv_gc = y_m
            ϵ = ϵ_m
            if need_save_gc
                isave, t_last, y_last = _hybrid_save!(
                    traj, tsave, plan, isave, t_last, y_last, t, xv_fo
                )
            end
            steps += 1
            switched_at_boundary = true
            break
        else
            t_hi = t_m
        end
    end
    if !switched_at_boundary
        xv_fo = _gc_to_full(xv_gc, Efunc, Bfunc, q, m, μ, t, phase)
    end
    return (
        t, t_hi, xv_gc, phase, Bmag_curr, ϵ, xv_fo,
        isave, t_last, y_last, steps, switched_at_boundary, ReturnCode.Success,
    )
end

function _hybrid_one(
        prob::FOProblem, i::Int, alg::AdaptiveHybrid,
        plan, isoutside::F,
        save_start::Bool, save_end::Bool, save_everystep::Bool, verbose::Bool,
        sim_seeds, integrator_gc, v_save::Val{SaveAdia},
        ::Val{KeepInteg} = Val(true)
    ) where {F, SaveAdia, KeepInteg}
    (; tspan, p) = prob
    q2m, m, Efunc, Bfunc = p.q2m, p.m, p.E, p.B
    p_fo = (; q2m, m, E = Efunc, B = Bfunc)
    q = q2m * m
    T = eltype(prob.u0)
    Tt = typeof(tspan[1])
    tdir = Tt(_span_direction(tspan))

    check_interval = alg.check_interval
    check_dt = Tt(alg.check_dt)
    n_correct = alg.switch_correction
    overshoot_limit = T(2) * T(alg.threshold_gc_to_fo)

    new_prob = if prob.prob_func === DEFAULT_PROB_FUNC
        prob
    else
        sim_seed = isnothing(sim_seeds) ? nothing : sim_seeds[i]
        rng = isnothing(sim_seed) ? default_rng() : Xoshiro(sim_seed)
        p_i = prob.prob_func(prob, EnsembleContext(i, 1, 0, sim_seed, rng, nothing))
        _check_hybrid_problem(p_i)
        p_i
    end
    xv_fo = SVector{6, T}(new_prob.u0)
    r = xv_fo[SVector(1, 2, 3)]
    v = xv_fo[SVector(4, 5, 6)]
    t = tspan[1]

    fields = _node_fields(p_fo, r, t)
    ϵ_init, Bmag_r = _fo_adia_value(r, v, fields.E(r, t), Efunc, Bfunc, q, m, t, alg)
    ϵ = T(ϵ_init)
    if ϵ < alg.threshold_gc_to_fo
        X_gc, vpar, μ_init, phase_init, valid_gc = _gc_switch_state_checked(
            xv_fo, Efunc, Bfunc, q, m, t, n_correct
        )
    else
        X_gc = r
        vpar = zero(T)
        μ_init = zero(T)
        phase_init = zero(T)
        valid_gc = false
    end
    is_adiabatic = valid_gc

    if is_adiabatic
        mode = :GC
        xv_gc = SVector{4, T}(X_gc[1], X_gc[2], X_gc[3], vpar)
        μ = T(μ_init)
        phase = T(phase_init)
        p_gc = (; q, q2m, μ, E = Efunc, B = Bfunc)

        Bmag_curr = T(norm(Bfunc(X_gc, t)))
        ω = abs(q2m * Bmag_curr)
        dt = tdir * Tt(0.5 * 2π / ω)
        if integrator_gc === nothing
            integrator_gc = init_gc_integrator(xv_gc, (t, tspan[2]), p_gc, dt, alg)
        else
            reinit_gc_integrator!(integrator_gc, xv_gc, (t, tspan[2]), p_gc, dt)
        end
        verbose && _log_mode("Initial mode: GC", ϵ, t)
    else
        mode = :FO
        xv_gc = zero(SVector{4, T})
        μ = zero(T)
        phase = zero(T)
        p_gc = (; q, q2m, μ, E = Efunc, B = Bfunc)
        dt = tdir * Tt(_fo_dt(alg, q2m, Bmag_r))
        v = update_velocity_half(v, r, dt, t, fields, Boris())
        Bmag_curr = zero(T)
        verbose && _log_mode("Initial mode: FO", ϵ, t)
    end

    initial_capacity = if use_saveat(plan)
        length(plan.times) + 2
    elseif !save_everystep
        2
    elseif is_adiabatic
        64
    else
        estimated_steps = ceil(Int, abs((tspan[2] - tspan[1]) / dt))
        max_cap = iszero(alg.threshold_fo_to_gc) ? 10000 : 512
        clamp(estimated_steps + div(estimated_steps, 10), 10, max_cap)
    end
    tsave = sizehint!(Vector{Tt}(undef, 0), initial_capacity)
    traj = sizehint!(Vector{SVector{6, T}}(undef, 0), initial_capacity)

    if save_start
        push!(traj, xv_fo)
        push!(tsave, t)
    end

    isave = 1
    t_last = t
    y_last = xv_fo
    steps = 0
    it = 1
    retcode = ReturnCode.Success
    t_last_check = t
    t_fo_holdoff = t
    check_scale = Tt(_check_scale(ϵ, alg, mode))
    # A step chosen by a switch must be the one actually taken, otherwise the
    # staggered velocity and the step disagree.
    dt_fixed = false

    adia_cap = SaveAdia ? max(8, div(initial_capacity, check_interval) + 2) : 0
    adia_t = SaveAdia ? sizehint!(Vector{Tt}(undef, 0), adia_cap) : Vector{Tt}(undef, 0)
    adia_vals = SaveAdia ? sizehint!(Vector{T}(undef, 0), adia_cap) : Vector{T}(undef, 0)
    adia_mode = SaveAdia ? sizehint!(Vector{Symbol}(undef, 0), adia_cap) :
        Vector{Symbol}(undef, 0)
    _record_adia!(v_save, adia_t, adia_vals, adia_mode, t, ϵ, mode)

    need_save_gc = use_saveat(plan) || save_everystep
    gc_chk = _GCCheckpoint(
        t, xv_gc, phase, Bmag_curr, ϵ, dt, steps, isave, t_last, y_last, length(traj)
    )

    while tdir * (tspan[2] - t) > zero(Tt) && steps < alg.maxiters
        if mode == :GC
            deadline = t_last_check + tdir * check_dt

            # Never step past the next check deadline: a long guiding-center
            # step could otherwise enter a non-adiabatic region and only be
            # examined after leaving it again.
            code = gc_step!(integrator_gc, deadline, isoutside, p_gc, tspan[2])
            if code != ReturnCode.Success
                retcode = code
                break
            end

            t_next = Tt(integrator_gc.t)
            y_next = SVector{4, T}(integrator_gc.u)
            it_next = it + 1
            it_gc = _gc_check_interval(ϵ, alg)
            # Also check right before clipping if the next natural GC step would
            # exceed `deadline`, avoiding a short remainder step just to check.
            dt_prop = SciMLBase.get_proposed_dt(integrator_gc)
            dt_rem = tspan[2] - t_next
            dt_try = tdir * dt_prop > tdir * dt_rem ? dt_rem : dt_prop
            scale = max(abs(t_next), abs(t_last_check), abs(deadline))
            tol = scale isa AbstractFloat ? 4 * eps(scale) : zero(scale)
            pre_deadline = tdir * (t_next - t_last_check) > zero(Tt) &&
                tdir * (t_next + dt_try - deadline) > tol

            if _check_due(it_next, t_next, t_last_check, it_gc, deadline, tdir) ||
                    pre_deadline
                mp_next = get_magnetic_properties(get_x(y_next), t_next, Bfunc)
                ϵ_next = T(_adia_value(mp_next, q, m, y_next[4], μ, alg))

                if ϵ_next ≥ alg.threshold_gc_to_fo
                    # issue #667: bisect/roll back when GC steps cross deep into
                    # a non-adiabatic region instead of crossing with frozen μ.
                    need_bisect = ϵ_next > overshoot_limit
                    switched_at_prev = false
                    xv_fo = zero(SVector{6, T})
                    t_hi = t_next
                    if t != gc_chk.t
                        mp_t = get_magnetic_properties(get_x(xv_gc), t, Bfunc)
                        ϵ_t = T(_adia_value(mp_t, q, m, xv_gc[4], μ, alg))
                        if ϵ_t > overshoot_limit
                            need_bisect = true
                            t_hi = t
                            t, xv_gc, phase, Bmag_curr, ϵ, dt, steps, isave,
                                t_last, y_last = _restore_gc_checkpoint!(
                                gc_chk, traj, tsave
                            )
                        elseif ϵ_t ≥ alg.threshold_gc_to_fo
                            ϵ = ϵ_t
                            xv_fo = _gc_to_full(
                                xv_gc, Efunc, mp_t.B, mp_t.Bmag, q, m, μ, t, phase
                            )
                            switched_at_prev = true
                        else
                            ϵ = ϵ_t
                        end
                    end

                    if switched_at_prev
                        t_fo_holdoff = t
                        clip_dt = true
                    elseif need_bisect
                        t, t_hi, xv_gc, phase, Bmag_curr, ϵ, xv_fo,
                            isave, t_last, y_last, steps, switched_at_boundary, code_b =
                            _bisect_gc_to_fo!(
                            integrator_gc, t, t_hi, xv_gc, phase, Bmag_curr,
                            ϵ, p_gc, q2m, q, m, μ, Efunc, Bfunc, alg,
                            tspan[2], isoutside, need_save_gc,
                            traj, tsave, plan, isave, t_last, y_last,
                            v_save, adia_t, adia_vals, adia_mode,
                            steps, overshoot_limit
                        )
                        if code_b == ReturnCode.Terminated
                            retcode = code_b
                            break
                        end
                        if steps >= alg.maxiters && !switched_at_boundary
                            break
                        end
                        t_fo_holdoff = switched_at_boundary ? t : t_hi
                        clip_dt = switched_at_boundary
                    else
                        Bmag_next = T(mp_next.Bmag)
                        phase, xv_fo = _advance_gc_state(
                            t_next, t_next - t, y_next, mp_next.B, Bmag_next,
                            phase, Bmag_curr, q2m, q, m, μ, Efunc, true
                        )
                        Bmag_curr = Bmag_next
                        t = t_next
                        xv_gc = y_next
                        dt = Tt(SciMLBase.get_proposed_dt(integrator_gc))
                        if need_save_gc
                            isave, t_last, y_last = _hybrid_save!(
                                traj, tsave, plan, isave, t_last, y_last, t, xv_fo
                            )
                        end
                        steps += 1
                        ϵ = ϵ_next
                        t_fo_holdoff = t
                        clip_dt = true
                    end

                    mode = :FO
                    it = 1
                    t_last_check = t
                    check_scale = Tt(_check_scale(ϵ, alg, mode))
                    verbose && _log_mode("Switch GC → FO", ϵ, t, get_x(xv_gc))
                    _record_adia!(v_save, adia_t, adia_vals, adia_mode, t, ϵ, :FO)
                    r, v, dt = _enter_fo_from_gc(
                        xv_fo, p_fo, q2m, t, dt, tdir, alg, clip_dt
                    )
                    dt_fixed = true
                    continue
                else
                    Bmag_next = T(mp_next.Bmag)
                    phase, y_out = _advance_gc_state(
                        t_next, t_next - t, y_next, mp_next.B, Bmag_next,
                        phase, Bmag_curr, q2m, q, m, μ, Efunc, need_save_gc
                    )
                    Bmag_curr = Bmag_next
                    t = t_next
                    xv_gc = y_next
                    dt = Tt(SciMLBase.get_proposed_dt(integrator_gc))
                    if need_save_gc
                        isave, t_last, y_last = _hybrid_save!(
                            traj, tsave, plan, isave, t_last, y_last, t, y_out
                        )
                    end
                    ϵ = ϵ_next
                    t_last_check = t
                    _record_adia!(v_save, adia_t, adia_vals, adia_mode, t, ϵ, :GC)
                    steps += 1
                    it = it_next
                    gc_chk = _GCCheckpoint(
                        t, xv_gc, phase, Bmag_curr, ϵ, dt,
                        steps, isave, t_last, y_last, length(traj)
                    )
                    continue
                end
            end

            B_next = Bfunc(get_x(y_next), t_next)
            Bmag_next = T(norm(B_next))
            phase, y_out = _advance_gc_state(
                t_next, t_next - t, y_next, B_next, Bmag_next,
                phase, Bmag_curr, q2m, q, m, μ, Efunc, need_save_gc
            )
            Bmag_curr = Bmag_next
            t = t_next
            xv_gc = y_next
            dt = Tt(SciMLBase.get_proposed_dt(integrator_gc))
            if need_save_gc
                isave, t_last, y_last = _hybrid_save!(
                    traj, tsave, plan, isave, t_last, y_last, t, y_out
                )
            end
            steps += 1
            it = it_next

        else # Mode == :FO
            it_eff = max(1, round(Int, check_scale * check_interval))
            deadline = t_last_check + tdir * check_scale * check_dt
            past_holdoff = tdir * (t - t_fo_holdoff) >= zero(Tt)
            check_now = past_holdoff && _check_due(
                it, t, t_last_check, it_eff, deadline, tdir
            )

            if check_now
                t_last_check = t
                E_r = SVector{3, T}(Efunc(r, t))
                mp_r = get_magnetic_properties(r, t, Bfunc)
                Bmag_r = T(mp_r.Bmag)
                fields = CachedFields(
                    q2m, Returns(E_r), Returns(SVector{3, T}(mp_r.B)), r, t, nothing
                )
                v_node = update_velocity_node(v, r, dt, t, fields, Boris())
                have_v_node = true
                short_circuit = !SaveAdia && alg.check_margin_max == 1
                ϵ = T(
                    _fo_adia_from_mp(
                        r, v_node, E_r, mp_r, Efunc, Bfunc, q, m, t, alg, short_circuit
                    )
                )
                check_scale = Tt(_check_scale(ϵ, alg, mode))
                deadline = t + tdir * check_scale * check_dt
                if ϵ < alg.threshold_fo_to_gc
                    xv_sync = vcat(r, v_node)
                    X_gc, vpar, μ_fo, phase_fo, valid_gc = _gc_switch_state_checked(
                        xv_sync, Efunc, Bfunc, q, m, t, n_correct
                    )
                    if valid_gc
                        mode = :GC
                        it = 1
                        check_scale = one(Tt)
                        verbose && _log_mode("Switch FO → GC", ϵ, t, r)
                        _record_adia!(v_save, adia_t, adia_vals, adia_mode, t, ϵ, :GC)
                        if need_save_gc
                            isave, t_last, y_last = _hybrid_save!(
                                traj, tsave, plan, isave, t_last, y_last, t, xv_sync
                            )
                        end
                        xv_gc = SVector{4, T}(X_gc[1], X_gc[2], X_gc[3], vpar)
                        μ = T(μ_fo)
                        phase = T(phase_fo)
                        p_gc = (; q, q2m, μ, E = Efunc, B = Bfunc)

                        Bmag_curr = T(norm(Bfunc(X_gc, t)))
                        # Carry the full-orbit step over instead of resetting
                        # it, capped at half a gyroperiod.
                        omega = abs(q2m * Bmag_curr)
                        dt = tdir * min(abs(dt), Tt(0.5 * 2π / omega))
                        if integrator_gc === nothing
                            integrator_gc = init_gc_integrator(
                                xv_gc, (t, tspan[2]), p_gc, dt, alg
                            )
                        else
                            reinit_gc_integrator!(
                                integrator_gc, xv_gc, (t, tspan[2]), p_gc, dt
                            )
                        end
                        gc_chk = _GCCheckpoint(
                            t, xv_gc, phase, Bmag_curr, ϵ, dt,
                            steps, isave, t_last, y_last, length(traj)
                        )
                        continue
                    end
                    check_scale = one(Tt)
                    deadline = t + tdir * check_dt
                    _record_adia!(v_save, adia_t, adia_vals, adia_mode, t, ϵ, :FO)
                else
                    _record_adia!(v_save, adia_t, adia_vals, adia_mode, t, ϵ, :FO)
                end
            else
                fields = _node_fields(p_fo, r, t)
                Bmag_r = T(norm(fields.B(r, t)))
                v_node = v
                have_v_node = false
            end

            dt_new = dt_fixed ? dt : tdir * Tt(_fo_dt(alg, q2m, Bmag_r))
            dt_fixed = false
            dt_rem = tspan[2] - t
            if tdir * dt_new > tdir * dt_rem
                dt_new = dt_rem
            end
            # Never step past the next check deadline (see the GC branch).
            dt_check = deadline - t
            if zero(dt_check) < tdir * dt_check < tdir * dt_new
                dt_new = dt_check
            end

            need_save = use_saveat(plan) || (save_everystep && it > 1)
            need_resync = dt_new != dt
            if (need_save || need_resync) && !have_v_node
                v_node = update_velocity_node(v, r, dt, t, fields, Boris())
            end
            if need_resync
                v = update_velocity_half(v_node, r, dt_new, t, fields, Boris())
                dt = dt_new
            end

            v_prev = v
            r_next, v = advance_boris(v, r, dt, t, fields, Boris())
            t_next = t + dt
            if isoutside(vcat(r_next, v), p, t_next)
                retcode = ReturnCode.Terminated
                v = v_prev
                break
            end

            if need_save
                isave, t_last, y_last = _hybrid_save!(
                    traj, tsave, plan, isave, t_last, y_last, t, vcat(r, v_node)
                )
            end

            r = r_next
            t = t_next
            it += 1
            steps += 1
        end
    end

    if mode == :FO && use_saveat(plan) && t != t_last
        v_final = update_velocity_node(v, r, dt, t, p_fo, Boris())
        y_final = vcat(r, v_final)
        isave, t_last, y_last = _hybrid_save!(
            traj, tsave, plan, isave, t_last, y_last, t, y_final
        )
        if save_end && (isempty(tsave) || tsave[end] != t)
            push!(traj, y_final)
            push!(tsave, t)
        end
    elseif save_end && (isempty(tsave) || tsave[end] != t)
        if mode == :GC
            push!(traj, _gc_to_full(xv_gc, Efunc, Bfunc, q, m, μ, t, phase))
        else
            v_final = update_velocity_node(v, r, dt, t, p_fo, Boris())
            push!(traj, vcat(r, v_final))
        end
        push!(tsave, t)
    end

    sol_alg = :hybrid
    interp = LinearInterpolation(tsave, traj)
    if steps >= alg.maxiters && retcode == ReturnCode.Success
        retcode = ReturnCode.MaxIters
    end
    stats = SaveAdia ?
        (adiabaticity = (t = adia_t, components = adia_vals, mode = adia_mode),) :
        nothing

    sol = build_solution(prob, sol_alg, tsave, traj; interp, retcode, stats)
    return KeepInteg ? (sol, integrator_gc) : sol
end

@inline function _hybrid_adaptive!(
        sols, prob::FOProblem, irange, alg::AdaptiveHybrid,
        plan, isoutside::F,
        save_start, save_end, save_everystep, verbose, sim_seeds, v_save::Val{SaveAdia}
    ) where {F, SaveAdia}
    integrator_gc = nothing
    @inbounds for i in irange
        sol, integrator_gc = _hybrid_one(
            prob, i, alg, plan, isoutside, save_start, save_end, save_everystep,
            verbose, sim_seeds, integrator_gc, v_save
        )
        sols[i] = sol
    end
    return
end

# --- GPU Field Adaptation ---

BorisPushers.adapt_field_to_gpu(field::Field, ::CPU) = field
BorisPushers.adapt_field_to_gpu(field::ZeroField, ::CPU) = field
BorisPushers.adapt_field_to_gpu(field::ZeroField, ::Backend) = field

function BorisPushers.adapt_field_to_gpu(field::Field, backend::Backend)
    backend isa CPU && return field

    raw = field.field_function
    adapted_func = raw isa Function ? BorisPushers.adapt_field_to_gpu(raw, backend) :
        Adapt.adapt(backend, raw)
    return Field{is_time_dependent(field), typeof(adapted_func)}(adapted_func)
end

BorisPushers.adapt_field_to_gpu(f::Function, backend::Backend) = Adapt.adapt(backend, f)

function BorisPushers.adapt_field_to_gpu(
        fi::Union{
            FieldInterpolator, FieldInterpolator2D, FieldInterpolator1D,
            SphericalFieldInterpolator,
        },
        backend::Backend
    )
    backend isa CPU && return fi
    return _to_gpu_grid(fi, backend)
end

@setup_workload begin
    @compile_workload begin
        x = range(-10, 10, length = 4)
        y = range(-10, 10, length = 6)
        z = range(-10, 10, length = 8)
        B = fill(0.0, 3, length(x), length(y), length(z))
        E = fill(0.0, 3, length(x), length(y), length(z))
        B[3, :, :, :] .= 10.0e-9
        E[3, :, :, :] .= 5.0e-10
        B_analytic(r) = SA[0.0, 0.0, 1.0e-9]
        E_analytic(r, t) = SA[1.0e-10, 0.0, 0.0]

        tspan = (0.0, 1.0)
        dt = 0.5
        # The slow particle stays inside the grid, where the interpolated field is finite.
        for (stateinit, param) in (
                (
                    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0], prepare(;
                        E = build_interpolator(E, (x, y, z); bc = TestParticle.FillExtrap(NaN)),
                        B = build_interpolator(B, (x, y, z); bc = TestParticle.FillExtrap(NaN))
                    ),
                ),
                ([0.0, 0.0, 0.0, 0.0, 1.0e5, 0.0], prepare(; B = B_analytic)),
                (
                    [0.0, 0.0, 0.0, 0.0, 1.0e5, 0.0],
                    prepare(; E = E_analytic, B = B_analytic),
                ),
            )
            prob = TraceProblem(stateinit, tspan, param)

            solve(prob, Boris(); dt, saveat = 100 * dt)
            solve(prob, Boris(), EnsembleThreads(); dt, saveat = 100 * dt)
            solve(prob, Boris(), CPU(); dt, saveat = 100 * dt)

            alg_adaptive = AdaptiveBoris(safety = 0.1)
            sol = solve(prob, alg_adaptive)
            get_fields(sol)
            get_work_rates(sol)

            solve(prob, MultistepBoris2(; n = 2); dt)

            # `threshold = 0` never leaves the full-orbit mode, so the hybrid loop
            # is compiled here without needing the guiding-center integrator that
            # `TestParticleVernerExt` supplies.
            solve(prob, AdaptiveHybrid(threshold = 0.0, dtmax = 1.0))
        end

        let
            tspan32 = (0.0f0, 1.0f0)
            dt32 = 0.5f0
            stateinit32 = Float32[0.0, 0.0, 0.0, 0.0, 1.0e5, 0.0]
            param32_B = prepare(; B = B_analytic, type = Float32)
            param32_EB = prepare(; E = E_analytic, B = B_analytic, type = Float32)

            for param in (param32_B, param32_EB)
                prob32 = TraceProblem(stateinit32, tspan32, param)
                solve(prob32, Boris(); dt = dt32, save_everystep = false)
                solve(
                    prob32, Boris(), EnsembleThreads();
                    dt = dt32, trajectories = 2, save_everystep = false
                )
                solve(prob32, Boris(), CPU(); dt = dt32, save_everystep = false)
            end
        end
    end
end

end # module TestParticleBorisPushersExt

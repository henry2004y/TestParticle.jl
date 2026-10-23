module TestParticleBorisPushersExt

using TestParticle
import TestParticle: TraceProblem, TraceHybridProblem, AdaptiveHybrid,
    SavingPlan, use_saveat, saveat_reached, saveat_interpolate, _span_direction,
    Field, ZeroField, FieldInterpolator, FieldInterpolator1D,
    FieldInterpolator2D, SphericalFieldInterpolator,
    GPUGrid1D, GPUGrid2D, GPUGrid3D, GPUSphericalGrid,
    _to_gpu_grid, _to_gpu_spherical_grid,
    get_x, get_magnetic_properties, get_work_rates, _promote_p,
    _gc_to_full, _get_gc_parameters, update_dp5,
    adiabaticity_components, _rng_seed, is_time_dependent, _foreach_chunk,
    ODE_DEFAULT_ISOUTOFDOMAIN

using BorisPushers
import BorisPushers: AbstractBoris, Boris, AdaptiveBoris, MultistepBoris,
    AdaptiveMultistepBoris, update_velocity, adapt_field_to_gpu

using SciMLBase: AbstractODEProblem, BasicEnsembleAlgorithm,
    DiscreteCallback, EnsembleProblem, EnsembleSolution, EnsembleContext,
    FullSpecialize, LinearInterpolation, ODEProblem, ReturnCode,
    build_solution, isadaptive, terminate!, EnsembleSerial, EnsembleThreads
import SciMLBase: solve

using LinearAlgebra: norm
using StaticArrays: SVector, SA
import Adapt
using KernelAbstractions: Backend, CPU
using Random: default_rng, Xoshiro
using PrecompileTools: @setup_workload, @compile_workload

# --- Boris ODE Problem Adaptation ---

_boris_rhs(u, p, t) = nothing

function _boris_problem(prob::TraceProblem)
    T = eltype(prob.u0)
    tspan = (T(prob.tspan[1]), T(prob.tspan[2]))
    p = _promote_p(T, prob.p)
    return ODEProblem{false, FullSpecialize}(
        _boris_rhs, SVector{6, T}(prob.u0), tspan, p
    )
end

function _boris_initial_dt(prob::TraceProblem, alg)
    q2m = prob.p.q2m
    u0 = prob.u0
    r0 = SVector{3}(u0[1], u0[2], u0[3])
    Bmag = norm(prob.p.B(r0, prob.tspan[1]))
    tdir = sign(prob.tspan[2] - prob.tspan[1])

    return tdir * (2π * alg.safety) / (abs(q2m) * Bmag)
end

function _boris_check_limits(prob::TraceProblem, dt, alg, maxiters)
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

@inline function _prepare_saved_data(
        xv, p, t, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    data = xv
    local E_field, magnetic_props

    if SaveFields
        r = get_x(xv)
        T = eltype(xv)
        E_field = SVector{3, T}(p.E(r, t))

        if SaveWork
            magnetic_props = get_magnetic_properties(r, t, p.B)
            B_vec = SVector{3, T}(magnetic_props[1])
            data = vcat(data, E_field, B_vec)
        else
            B_vec = SVector{3, T}(p.B(r, t))
            data = vcat(data, E_field, B_vec)
        end
    end

    if SaveWork
        if SaveFields
            work = get_work_rates(xv, p, t, magnetic_props, E_field)
        else
            work = get_work_rates(xv, p, t)
        end
        data = vcat(data, work)
    end

    return data
end

_boris_finalize(sol, p, ::Val{false}, ::Val{false}) = sol

function _boris_finalize(
        sol, p, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    u = [
        _prepare_saved_data(sol.u[i], p, sol.t[i], Val(SaveFields), Val(SaveWork))
            for i in eachindex(sol.t)
    ]

    return build_solution(
        sol.prob, sol.alg, sol.t, u;
        interp = LinearInterpolation(sol.t, u), retcode = sol.retcode, stats = nothing
    )
end

function SciMLBase.solve(
        prob::TraceProblem, alg::AbstractBoris;
        saveat = (),
        dt::Union{Nothing, Real} = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true,
        save_end::Bool = true,
        save_everystep::Bool = isempty(saveat),
        save_fields::Bool = false,
        save_work::Bool = false,
        maxiters::Int = 1_000_000,
        rng = nothing,
        seed = nothing,
        kwargs...
    ) where {F}
    T = eltype(prob.u0)
    step = dt === nothing ? T(_boris_initial_dt(prob, alg)) : T(dt)
    _boris_check_limits(prob, step, alg, maxiters)

    ode_prob = _boris_problem(prob)
    callback = _boris_callback(isoutside)

    sol = if isempty(saveat)
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

    return if save_fields
        save_work ? _boris_finalize(sol, prob.p, Val(true), Val(true)) :
            _boris_finalize(sol, prob.p, Val(true), Val(false))
    else
        save_work ? _boris_finalize(sol, prob.p, Val(false), Val(true)) :
            _boris_finalize(sol, prob.p, Val(false), Val(false))
    end
end

function SciMLBase.solve(
        prob::TraceProblem, alg::AbstractBoris,
        ensemblealg::BasicEnsembleAlgorithm;
        trajectories::Int = 1,
        seed::Union{Nothing, Integer} = nothing,
        safetycopy::Bool = false,
        kwargs...
    )
    ensemble_prob = EnsembleProblem(
        prob; prob_func = prob.prob_func, safetycopy
    )

    return solve(ensemble_prob, alg, ensemblealg; trajectories, seed, kwargs...)
end

# --- AdaptiveHybrid Solver Implementation ---

@inline function SciMLBase.solve(
        prob::TraceHybridProblem, alg::AdaptiveHybrid,
        ensemblealg::EA = EnsembleSerial();
        trajectories::Int = 1,
        saveat = (),
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true,
        save_everystep::Bool = true, verbose::Bool = false,
        seed::Union{Nothing, Integer} = nothing
    ) where {EA <: BasicEnsembleAlgorithm, F}
    plan = SavingPlan(
        saveat, prob.tspan, _span_direction(prob.tspan), typeof(prob.tspan[1])
    )
    return _solve(
        ensemblealg, prob, trajectories, alg,
        plan, isoutside,
        save_start, save_end, save_everystep, verbose, seed
    )
end

@inline function _solve(
        ensemblealg::Union{EnsembleSerial, EnsembleThreads}, prob::TraceHybridProblem,
        trajectories, alg::AdaptiveHybrid,
        plan, isoutside::F,
        save_start, save_end, save_everystep, verbose, seed
    ) where {F}
    sample_stats = _adia_sample_stats(prob, alg.save_adiabaticity)
    sol_type = _get_sol_type(prob, zero(eltype(prob.tspan)), sample_stats)
    sols = Vector{sol_type}(undef, trajectories)

    elapsed_time = @elapsed _foreach_chunk(ensemblealg, trajectories) do irange
        _hybrid_adaptive!(
            sols, prob, irange, alg, plan,
            isoutside, save_start, save_end, save_everystep, verbose, seed
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

    sol = build_solution(prob, alg, t, u; interp, stats)
    return typeof(sol)
end

@inline function _adia_select(comps, alg::AdaptiveHybrid)
    if alg.adiabaticity === :curvature
        return comps[1]
    elseif alg.adiabaticity === :gradB
        return comps[2]
    elseif alg.adiabaticity === :jacobian
        return comps[3]
    else
        return max(comps[1], comps[2])
    end
end

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

@inline function _fo_dt(alg, q2m, Bfunc, r, t = 0.0)
    Bmag = norm(Bfunc(r, t))
    ω = abs(q2m) * Bmag
    dt = 2π * alg.safety_fo / ω
    return clamp(dt, alg.dtmin, alg.dtmax)
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

@inline function _hybrid_adaptive!(
        sols, prob::TraceHybridProblem, irange, alg,
        plan, isoutside::F,
        save_start, save_end, save_everystep, verbose, seed
    ) where {F}
    (; tspan, p, u0) = prob
    q2m, m, Efunc, Bfunc = p.q2m, p.m, p.E, p.B
    q = q2m * m
    T = eltype(u0)

    is_td = is_time_dependent(Efunc) || is_time_dependent(Bfunc)

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

        X_gc, vpar, _, _, μ, phase = _get_gc_parameters(xv_fo, Efunc, Bfunc, q, m, t)

        comps0 = adiabaticity_components(X_gc, Bfunc, q, m, μ, t)
        ϵ = _adia_select(comps0, alg)
        is_adiabatic = ϵ < alg.threshold_gc_to_fo

        if is_adiabatic
            mode = :GC
            xv_gc = SVector{4, T}(X_gc[1], X_gc[2], X_gc[3], vpar)
            p_gc = (; q, q2m, μ, E = Efunc, B = Bfunc)

            B_vec = Bfunc(X_gc, t)
            Bmag = norm(B_vec)
            Bmag_curr = Bmag
            ω = abs(q2m * Bmag)
            dt = 0.5 * 2π / ω
            verbose && @info "Initial mode: GC" ϵ t
        else
            mode = :FO
            dt = _fo_dt(alg, q2m, Bfunc, r, t)
            v = update_velocity(v, r, -0.5 * dt, t, p, Boris())
            Bmag_curr = zero(T)
            verbose && @info "Initial mode: FO" ϵ t
        end

        estimated_steps = ceil(Int, (tspan[2] - tspan[1]) / dt)
        initial_capacity = if use_saveat(plan)
            length(plan.times) + 2
        elseif !save_everystep
            2
        else
            min(max(10, estimated_steps + div(estimated_steps, 10)), 10000)
        end
        tsave = Vector{T}(undef, 0)
        sizehint!(tsave, initial_capacity)
        traj = Vector{SVector{6, T}}(undef, 0)
        sizehint!(traj, initial_capacity)

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

        adia_t = Vector{typeof(tspan[1])}(undef, 0)
        adia_vals = Vector{T}(undef, 0)
        adia_mode = Vector{Symbol}(undef, 0)
        if alg.save_adiabaticity
            adia_cap = max(10, div(initial_capacity, alg.check_interval) + 2)
            sizehint!(adia_t, adia_cap)
            sizehint!(adia_vals, adia_cap)
            sizehint!(adia_mode, adia_cap)
            push!(adia_t, t)
            push!(adia_vals, _adia_select(comps0, alg))
            push!(adia_mode, mode)
        end

        while t < tspan[2] && steps < alg.maxiters
            if mode == :GC
                if it % alg.check_interval == 0 && t != t_last_check
                    t_last_check = t
                    comps = adiabaticity_components(
                        xv_gc[SVector(1, 2, 3)], Bfunc, q, m, μ, t
                    )
                    ϵ = _adia_select(comps, alg)
                    if ϵ ≥ alg.threshold_gc_to_fo
                        mode = :FO
                        it = 1
                        verbose && @info "Switch GC → FO" ϵ t r = xv_gc[SVector(1, 2, 3)]
                        if alg.save_adiabaticity
                            push!(adia_t, t)
                            push!(adia_vals, ϵ)
                            push!(adia_mode, :FO)
                        end
                        xv_fo = _gc_to_full(
                            xv_gc, Efunc, Bfunc, q, m, μ, t, phase
                        )
                        r = xv_fo[SVector(1, 2, 3)]
                        v = xv_fo[SVector(4, 5, 6)]

                        dt = _fo_dt(alg, q2m, Bfunc, r, t)
                        v = update_velocity(v, r, -0.5 * dt, t, p, Boris())
                        continue
                    elseif alg.save_adiabaticity
                        push!(adia_t, t)
                        push!(adia_vals, ϵ)
                        push!(adia_mode, :GC)
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
                    if isoutside(y_next, p_gc, t + dt)
                        retcode = ReturnCode.Terminated
                        break
                    end

                    t += dt
                    xv_gc = y_next
                    Bmag_next = norm(Bfunc(get_x(xv_gc), t))
                    phase = mod2pi(phase - 0.5 * dt * (q2m * (Bmag_curr + Bmag_next)))
                    Bmag_curr = Bmag_next

                    if use_saveat(plan) || save_everystep
                        y_out = _gc_to_full(xv_gc, Efunc, Bfunc, q, m, μ, t, phase)
                        isave, t_last, y_last = _hybrid_save!(
                            traj, tsave, plan, isave, t_last, y_last, t, y_out
                        )
                    end
                    steps += 1
                    it += 1
                end

                scale = error_ratio == 0.0 ?
                    max_growth : safety_gc * (1.0 / error_ratio)^0.2
                scale = max(min_growth, min(scale, max_growth))
                dt *= scale
                if t < tspan[2] && dt < 1.0e-14
                    retcode = ReturnCode.DtLessThanMin
                    break
                end

            else # Mode == :FO
                if it % alg.check_interval == 0 && t != t_last_check
                    t_last_check = t
                    t_sync = is_td ? t : zero(T)
                    v_sync = update_velocity(v, r, 0.5 * dt, t_sync, p, Boris())
                    xv_sync = SVector{6, T}(
                        r[1], r[2], r[3], v_sync[1], v_sync[2], v_sync[3]
                    )
                    X_gc, vpar, _, _, μ_fo, phase_fo = _get_gc_parameters(
                        xv_sync, Efunc, Bfunc, q, m, t
                    )
                    comps = adiabaticity_components(X_gc, Bfunc, q, m, μ_fo, t)
                    ϵ = _adia_select(comps, alg)
                    if ϵ < alg.threshold_fo_to_gc
                        mode = :GC
                        it = 1
                        verbose && @info "Switch FO → GC" ϵ t r = r
                        if alg.save_adiabaticity
                            push!(adia_t, t)
                            push!(adia_vals, ϵ)
                            push!(adia_mode, :GC)
                        end
                        if use_saveat(plan) || save_everystep
                            isave, t_last, y_last = _hybrid_save!(
                                traj, tsave, plan, isave, t_last, y_last, t, xv_sync
                            )
                        end
                        xv_gc = SVector{4, T}(X_gc[1], X_gc[2], X_gc[3], vpar)
                        p_gc = (; q, q2m, μ = μ_fo, E = Efunc, B = Bfunc)
                        μ = μ_fo
                        phase = phase_fo

                        Bmag_curr = norm(Bfunc(X_gc, t))
                        omega = abs(q2m * Bmag_curr)
                        dt = 0.5 * 2π / omega
                        continue
                    elseif alg.save_adiabaticity
                        push!(adia_t, t)
                        push!(adia_vals, ϵ)
                        push!(adia_mode, :FO)
                    end
                end

                if t + dt > tspan[2]
                    dt = tspan[2] - t
                end

                v_prev = v
                v = update_velocity(v, r, dt, t + 0.5 * dt, p, Boris())
                r_next = r + v * dt
                t_next = t + dt
                if isoutside(vcat(r_next, v), p, t_next)
                    retcode = ReturnCode.Terminated
                    v = v_prev
                    break
                end

                if use_saveat(plan) || (save_everystep && (it - 1) > 0)
                    v_save = update_velocity(v_prev, r, 0.5 * dt, t, p, Boris())
                    isave, t_last, y_last = _hybrid_save!(
                        traj, tsave, plan, isave, t_last, y_last, t, vcat(r, v_save)
                    )
                end

                r = r_next
                t = t_next
                it += 1
                steps += 1
            end
        end

        if save_end && (isempty(tsave) || tsave[end] != t)
            if mode == :GC
                push!(traj, _gc_to_full(xv_gc, Efunc, Bfunc, q, m, μ, t, phase))
            else
                t_final = is_td ? t : zero(T)
                v_final = update_velocity(v, r, 0.5 * dt, t_final, p, Boris())
                push!(traj, vcat(r, v_final))
            end
            push!(tsave, t)
        end

        sol_alg = :hybrid
        interp = LinearInterpolation(tsave, traj)
        if steps >= alg.maxiters && retcode == ReturnCode.Success
            retcode = ReturnCode.MaxIters
        end
        stats = alg.save_adiabaticity ?
            (adiabaticity = (t = adia_t, components = adia_vals, mode = adia_mode),) :
            nothing

        sols[i] = build_solution(
            prob, sol_alg, tsave, traj; interp, retcode, stats
        )
    end
    return
end

# --- GPU Field Adaptation ---

BorisPushers.adapt_field_to_gpu(field::Field, ::CPU) = field
BorisPushers.adapt_field_to_gpu(field::ZeroField, ::CPU) = field
BorisPushers.adapt_field_to_gpu(field::ZeroField, ::Backend) = field

function BorisPushers.adapt_field_to_gpu(field::Field, backend::Backend)
    backend isa CPU && return field

    adapted_func = BorisPushers.adapt_field_to_gpu(field.field_function, backend)
    return Field{is_time_dependent(field), typeof(adapted_func)}(adapted_func)
end

BorisPushers.adapt_field_to_gpu(f::Function, backend::Backend) = Adapt.adapt(backend, f)

function BorisPushers.adapt_field_to_gpu(fi::FieldInterpolator, backend::Backend)
    backend isa CPU && return fi
    return _to_gpu_grid(fi.itp, backend)
end

function BorisPushers.adapt_field_to_gpu(fi::FieldInterpolator2D, backend::Backend)
    backend isa CPU && return fi
    return _to_gpu_grid(fi.itp, backend)
end

function BorisPushers.adapt_field_to_gpu(fi::FieldInterpolator1D, backend::Backend)
    backend isa CPU && return fi
    return _to_gpu_grid(fi.itp, backend; dir = fi.dir)
end

function BorisPushers.adapt_field_to_gpu(
        fi::SphericalFieldInterpolator, backend::Backend
    )
    backend isa CPU && return fi
    return _to_gpu_spherical_grid(fi.itp, backend)
end

BorisPushers.adapt_field_to_gpu(g::GPUGrid3D, backend::Backend) =
    backend isa CPU ? g : Adapt.adapt(backend, g)
BorisPushers.adapt_field_to_gpu(g::GPUGrid2D, backend::Backend) =
    backend isa CPU ? g : Adapt.adapt(backend, g)
BorisPushers.adapt_field_to_gpu(g::GPUGrid1D, backend::Backend) =
    backend isa CPU ? g : Adapt.adapt(backend, g)
BorisPushers.adapt_field_to_gpu(g::GPUSphericalGrid, backend::Backend) =
    backend isa CPU ? g : Adapt.adapt(backend, g)

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

        tspan = (0.0, 1.0)
        dt = 0.5
        # The slow particle stays inside the grid, where the interpolated field is finite.
        for (stateinit, param) in (
                ([0.0, 0.0, 0.0, 0.0, 1.0, 0.0], prepare(x, y, z, E, B)),
                ([0.0, 0.0, 0.0, 0.0, 1.0e5, 0.0], prepare(B_analytic)),
            )
            prob = TraceProblem(stateinit, tspan, param)

            solve(prob, Boris(); dt, saveat = 100 * dt)
            solve(prob, Boris(), EnsembleThreads(); dt, saveat = 100 * dt)
            solve(prob, Boris(), CPU(); dt, saveat = 100 * dt)

            alg_adaptive = AdaptiveBoris(safety = 0.1)
            sol = solve(prob, alg_adaptive)
            solve(prob, alg_adaptive; save_fields = true, save_work = true)
            get_fields(sol)
            get_work(sol)

            solve(prob, MultistepBoris2(; n = 2); dt)
            solve(prob, MultistepBoris2(; n = 2); dt, save_fields = true)

            prob_hybrid = TraceHybridProblem(stateinit, tspan, param)
            solve(prob_hybrid, AdaptiveHybrid(threshold = 0.1, dtmax = 1.0))
        end
    end
end

end # module TestParticleBorisPushersExt

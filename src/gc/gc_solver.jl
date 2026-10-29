# Native GC solver

"""
    GCRK4()

Fixed-step RK4 native solver of [`TraceGCProblem`](@ref); needs `dt`.
"""
struct GCRK4 end

"""
    GCVern7(; abstol = 1e-6, reltol = 1e-6)

Adaptive native solver of [`TraceGCProblem`](@ref); `dt` is the optional initial step.
Each accepted step is taken by `Vern7` from OrdinaryDiffEqVerner, so it requires
`using OrdinaryDiffEqVerner`. That package is a weak dependency to keep it out of
the import path of users who never trace a guiding center; without it this solver
raises an `ArgumentError` and [`GCRK4`](@ref) remains available.
"""
struct GCVern7{T}
    abstol::T
    reltol::T
end

GCVern7(; abstol = 1.0e-6, reltol = 1.0e-6) = GCVern7(promote(abstol, reltol)...)

const GCAlgorithm = Union{GCRK4, GCVern7}

"""
    solve(prob::TraceGCProblem, alg::Union{GCRK4, GCVern7}, ensemblealg = EnsembleSerial();
        trajectories = 1, dt = nothing, saveat = (), isoutside = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start = true, save_end = true, save_everystep = true, maxiters = 10000,
        save_fields = false, save_work = false, seed = nothing)

Trace guiding centers with a native solver.

`saveat` names the times at which to report the state, as a collection or as an
interval; the state is interpolated linearly between the two steps that bracket
each time, leaving the trajectory unchanged.
"""
function solve(
        prob::TraceGCProblem, alg::GCAlgorithm,
        ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1,
        saveat = (),
        dt::Union{AbstractFloat, Nothing} = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        maxiters = 10000, save_fields::Bool = false, save_work::Bool = false,
        seed::Union{Nothing, Integer} = nothing
    ) where {F}
    plan = SavingPlan(saveat, prob.tspan, _gc_direction(prob), _gc_time_type(prob, dt))

    # Ternaries give `Union{Val{true}, Val{false}}`, which inference splits into static
    # calls; `Val(::Bool)` would dispatch dynamically and allocate on every solve.
    return _solve(
        ensemblealg, prob, trajectories, dt, plan, isoutside,
        save_start, save_end, save_everystep, alg, maxiters,
        save_fields ? Val(true) : Val(false), save_work ? Val(true) : Val(false), seed
    )
end

function solve(
        prob::TraceGCProblem, ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        alg::Symbol, abstol = 1.0e-6, reltol = 1.0e-6, kwargs...
    )
    gcalg = if alg === :rk4
        GCRK4()
    elseif alg === :rk45
        GCVern7(; abstol, reltol)
    else
        throw(ArgumentError("alg must be :rk4 or :rk45 for TraceGCProblem; got $(repr(alg))."))
    end
    Base.depwarn(
        "`solve(prob::TraceGCProblem; alg = $(repr(alg)))` is deprecated, use `solve(prob, $(nameof(typeof(gcalg)))(); ...)` instead.",
        :solve
    )
    return solve(prob, gcalg, ensemblealg; kwargs...)
end

function _solve(
        ensemblealg::Union{EnsembleSerial, EnsembleThreads}, prob::TraceGCProblem,
        trajectories, dt, plan, isoutside::F, save_start, save_end, save_everystep, alg,
        maxiters, ::Val{SaveFields}, ::Val{SaveWork}, seed
    ) where {F, SaveFields, SaveWork}
    sols, nt, nout = _prepare_gc(
        prob, trajectories, dt, plan,
        save_start, save_end, save_everystep, alg, maxiters, Val(SaveFields), Val(SaveWork)
    )

    elapsed_time = @elapsed _foreach_chunk(ensemblealg, trajectories) do irange
        if alg isa GCVern7
            _vern7!(
                sols, prob, irange, alg, plan, dt, isoutside,
                save_start, save_end, save_everystep, maxiters,
                Val(SaveFields), Val(SaveWork), seed
            )
        else
            _rk4!(
                sols, prob, irange, plan, dt, nt, nout, isoutside,
                save_start, save_end, save_everystep, maxiters,
                Val(SaveFields), Val(SaveWork), seed
            )
        end
    end

    return EnsembleSolution(sols, elapsed_time, true)
end

"""
    _gc_direction(prob) -> Int

`+1` for a forward time span and `-1` for a backward one.
"""
_gc_direction(prob::TraceGCProblem) = _span_direction(prob.tspan)

"""
    _gc_time_type(prob, dt) -> Type

The type of the time values a solve will produce.
"""
function _gc_time_type(prob::TraceGCProblem, dt)
    dt_guess = isnothing(dt) ? one(eltype(prob.u0)) : dt
    return typeof(prob.tspan[1] + dt_guess)
end

function _get_sol_type(prob::TraceGCProblem, dt, alg, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    u0 = prob.u0
    tspan = prob.tspan
    T_t = _gc_time_type(prob, dt)
    t = Vector{T_t}(undef, 0)
    # Force u to be Vector{SVector{4, T}}
    T = eltype(u0)

    n_vars = 4
    if SaveFields
        n_vars += 6
    end
    if SaveWork
        n_vars += 4
    end

    u = Vector{SVector{n_vars, T}}(undef, 0)
    interp = LinearInterpolation(t, u)

    sol = build_solution(prob, alg, t, u; interp = interp)
    return typeof(sol)
end

function _prepare_gc(
        prob::TraceGCProblem, trajectories, dt, plan,
        save_start, save_end, save_everystep, alg, maxiters, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    ttotal = prob.tspan[2] - prob.tspan[1]
    if isnothing(dt)
        if alg isa GCRK4
            error("dt must be provided for the fixed step solver GCRK4")
        end
        nt = 0
    else
        if alg isa GCRK4 && sign(dt) != sign(ttotal) && ttotal != 0
            error("dt must have the same sign as tspan[2] - tspan[1] for the fixed step solver GCRK4")
        end
        nt = round(Int, ttotal / dt)
        if alg isa GCRK4 && nt > maxiters
            nt = maxiters
        end
    end

    nout = save_start ? 1 : 0

    if use_saveat(plan)
        # One slot per requested time, plus the end of the run. A `maxiters` cut
        # short can leave some of them unreached, which the caller resizes away.
        nout += length(plan.times) + (save_end ? 1 : 0)
    elseif alg isa GCRK4 && save_everystep
        # For GCVern7, nout stays at the count above since we don't know the exact steps.
        last_is_step = nt > 0
        nout += nt
        if !save_end && last_is_step
            nout -= 1
        end
        if save_end && !last_is_step
            nout += 1
        end
    elseif alg isa GCRK4 && save_end
        nout += 1
    end

    sol_type = _get_sol_type(prob, dt, alg, Val(SaveFields), Val(SaveWork))
    sols = Vector{sol_type}(undef, trajectories)

    return sols, nt, nout
end

@inline function _prepare_saved_data_gc(xv, p, t, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    data = xv
    if SaveFields
        r = get_x(xv)
        E = p.E(r, t)
        B = p.B(r, t)
        data = vcat(data, E, B)
    end
    if SaveWork
        data = vcat(data, get_work_rates(xv, p, t))
    end
    return data
end

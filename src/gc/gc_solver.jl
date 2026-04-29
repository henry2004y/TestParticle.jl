# Native GC solver

"""
    TraceGCProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <: AbstractODEProblem{uType, tType, isinplace}

Problem type for tracing guiding centers.
"""
struct TraceGCProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <:
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

function TraceGCProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    _f = ODEFunction{true, DEFAULT_SPECIALIZATION}(x -> nothing) # dummy func
    return TraceGCProblem{
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
function TraceGCProblem{iip}(; f, u0, tspan, p, prob_func) where {iip}
    return TraceGCProblem{
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


"""
    solve(prob::TraceGCProblem; trajectories::Int=1, dt::AbstractFloat,
        savestepinterval::Int=1, isoutside::Function=ODE_DEFAULT_ISOUTOFDOMAIN,
        alg::Symbol=:rk4, abstol=1e-6, reltol=1e-6, maxiters=10000,
        save_fields::Bool=false, save_work::Bool=false)

Trace guiding centers using the RK4 method with specified `prob`.
If `alg` is `:rk45`, uses adaptive time stepping.
"""
function solve(
        prob::TraceGCProblem, ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1,
        dt::Union{AbstractFloat, Nothing} = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        alg::Symbol = :rk4, abstol = 1.0e-6, reltol = 1.0e-6, maxiters = 10000,
        save_fields::Bool = false, save_work::Bool = false
    ) where {F}
    if alg != :rk4 && alg != :rk45
        @warn "Only :rk4 and :rk45 are supported for native TraceGCProblem currently. Using :rk4."
    end

    return if save_fields
        if save_work
            return _solve(
                ensemblealg, prob, trajectories, dt, savestepinterval, isoutside,
                save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
                Val(true), Val(true)
            )
        else
            return _solve(
                ensemblealg, prob, trajectories, dt, savestepinterval, isoutside,
                save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
                Val(true), Val(false)
            )
        end
    else
        if save_work
            return _solve(
                ensemblealg, prob, trajectories, dt, savestepinterval, isoutside,
                save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
                Val(false), Val(true)
            )
        else
            return _solve(
                ensemblealg, prob, trajectories, dt, savestepinterval, isoutside,
                save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
                Val(false), Val(false)
            )
        end
    end
end

function _solve(
        ::EnsembleSerial, prob::TraceGCProblem, trajectories, dt, savestepinterval,
        isoutside, save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    sols, nt, nout = _prepare_gc(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, alg, maxiters, Val(SaveFields), Val(SaveWork)
    )
    irange = 1:trajectories

    elapsed_time = @elapsed if alg == :rk45
        _rk45!(
            sols, prob, irange, dt, isoutside,
            save_start, save_end, save_everystep, abstol, reltol, maxiters,
            Val(SaveFields), Val(SaveWork)
        )
    else
        _rk4!(
            sols, prob, irange, savestepinterval, dt, nt, nout, isoutside,
            save_start, save_end, save_everystep, maxiters,
            Val(SaveFields), Val(SaveWork)
        )
    end

    return EnsembleSolution(sols, elapsed_time, true)
end

function _solve(
        ::EnsembleThreads, prob::TraceGCProblem, trajectories, dt, savestepinterval,
        isoutside, save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    sols, nt,
        nout = _prepare_gc(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, alg, maxiters, Val(SaveFields), Val(SaveWork)
    )

    nchunks = Threads.nthreads()
    elapsed_time = @elapsed Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        if alg == :rk45
            _rk45!(
                sols, prob, irange, dt, isoutside,
                save_start, save_end, save_everystep, abstol, reltol, maxiters,
                Val(SaveFields), Val(SaveWork)
            )
        else
            _rk4!(
                sols, prob, irange, savestepinterval, dt, nt, nout, isoutside,
                save_start, save_end, save_everystep, maxiters,
                Val(SaveFields), Val(SaveWork)
            )
        end
    end

    return EnsembleSolution(sols, elapsed_time, true)
end

function _get_sol_type(prob::TraceGCProblem, dt, alg, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    u0 = prob.u0
    tspan = prob.tspan
    dt_guess = isnothing(dt) ? one(eltype(u0)) : dt
    T_t = typeof(tspan[1] + dt_guess)
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
        prob::TraceGCProblem, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, alg, maxiters, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    ttotal = prob.tspan[2] - prob.tspan[1]
    if isnothing(dt)
        if alg == :rk4
            error("dt must be provided for fixed step solver :rk4")
        end
        nt = 0
    else
        if alg == :rk4 && sign(dt) != sign(ttotal) && ttotal != 0
            error("dt must have the same sign as tspan[2] - tspan[1] for fixed step solver :rk4")
        end
        nt = round(Int, ttotal / dt)
        if alg == :rk4 && nt > maxiters
            nt = maxiters
        end
    end

    nout = 0
    if save_start
        nout += 1
    end

    # For :rk45, we initialize nout=0 since we don't know the exact steps.
    if alg == :rk4 && save_everystep
        steps = nt ÷ savestepinterval
        last_is_step = (nt > 0) && (nt % savestepinterval == 0)
        nout += steps
        if !save_end && last_is_step
            nout -= 1
        end
        if save_end && !last_is_step
            nout += 1
        end
    elseif alg == :rk4 && save_end
        nout += 1
    end

    sol_type = _get_sol_type(prob, dt, alg, Val(SaveFields), Val(SaveWork))
    sols = Vector{sol_type}(undef, trajectories)

    return sols, nt, nout
end

@inline function _prepare_saved_data_gc(xv, p, t, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    data = xv
    if SaveFields
        # p = (q, q2m, μ, Efunc, Bfunc)
        r = get_x(xv)
        T = eltype(xv)
        E = p[4](r, t)
        B = p[5](r, t)
        data = vcat(data, E, B)
    end
    if SaveWork
        work = get_work_rates_gc(xv, p, t)
        data = vcat(data, work)
    end
    return data
end



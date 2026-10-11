# TraceProblem constructors and SciML interface integration

@inline function _promote_p(::Type{T}, p::Tuple) where {T <: AbstractFloat}
    if length(p) >= 2 && p[1] isa Number && p[2] isa Number
        return (T(p[1]), T(p[2]), Base.tail(Base.tail(p))...)
    end
    return p
end
@inline _promote_p(::Type{T}, p) where {T} = p

@inline function _promote_trace_args(u0, tspan, p)
    T = eltype(u0)
    if T <: AbstractFloat
        tspan_T = (T(tspan[1]), T(tspan[2]))
        p_T = _promote_p(T, p)
        return tspan_T, p_T
    else
        return tspan, p
    end
end

"""
    TraceProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    TraceProblem(f, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)

Problem specification for particle tracing. Supports both Boris solvers and standard ODE
solvers (e.g. Tsit5, Vern7). When `f` is not provided, defaults to `trace!` (in-place)
for mutable arrays or `trace` (out-of-place) for `StaticArray`s.
"""
function TraceProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    tspan, p = _promote_trace_args(u0, tspan, p)
    isinplace = !(u0 isa StaticArray)
    _func = isinplace ? trace! : trace
    _f = ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(_func)
    return TraceProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function TraceProblem(f::Function, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    tspan, p = _promote_trace_args(u0, tspan, p)
    isinplace = !(u0 isa StaticArray)
    _f = ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(f)
    return TraceProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function TraceProblem(
        f::AbstractODEFunction{iip}, u0, tspan, p;
        prob_func = DEFAULT_PROB_FUNC
    ) where {iip}
    tspan, p = _promote_trace_args(u0, tspan, p)
    return TraceProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function TraceProblem{iip}(; f, u0, tspan, p, prob_func = DEFAULT_PROB_FUNC) where {iip}
    tspan, p = _promote_trace_args(u0, tspan, p)
    return TraceProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function SciMLBase.remake(
        prob::TraceProblem;
        f = prob.f,
        u0 = prob.u0,
        tspan = prob.tspan,
        p = prob.p,
        prob_func = prob.prob_func,
        kwargs...
    )
    tspan, p = _promote_trace_args(u0, tspan, p)
    isinplace = !(u0 isa StaticArray)
    _f = if f === prob.f && (prob.f.f === trace! || prob.f.f === trace)
        isinplace ? ODEFunction{true, DEFAULT_SPECIALIZATION}(trace!) :
            ODEFunction{false, DEFAULT_SPECIALIZATION}(trace)
    else
        f
    end
    return TraceProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function solve(
        prob::TraceProblem, alg::SciMLBase.AbstractODEAlgorithm,
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

function TraceCanonicalProblem(
        u0, tspan, p;
        relativistic::Bool = false,
        prob_func = DEFAULT_PROB_FUNC
    )
    tspan, p = _promote_trace_args(u0, tspan, p)
    isinplace = !(u0 isa StaticArray)
    _func = if relativistic
        isinplace ? trace_canonical_relativistic! : trace_canonical_relativistic
    else
        isinplace ? trace_canonical! : trace_canonical
    end
    _f = ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(_func)
    return TraceCanonicalProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function TraceCanonicalProblem(f::Function, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    tspan, p = _promote_trace_args(u0, tspan, p)
    isinplace = !(u0 isa StaticArray)
    _f = ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(f)
    return TraceCanonicalProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function TraceCanonicalProblem(
        f::AbstractODEFunction{iip}, u0, tspan, p;
        prob_func = DEFAULT_PROB_FUNC
    ) where {iip}
    tspan, p = _promote_trace_args(u0, tspan, p)
    return TraceCanonicalProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function TraceCanonicalProblem{iip}(; f, u0, tspan, p, prob_func = DEFAULT_PROB_FUNC) where {iip}
    tspan, p = _promote_trace_args(u0, tspan, p)
    return TraceCanonicalProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function SciMLBase.remake(
        prob::TraceCanonicalProblem;
        f = prob.f,
        u0 = prob.u0,
        tspan = prob.tspan,
        p = prob.p,
        prob_func = prob.prob_func,
        kwargs...
    )
    tspan, p = _promote_trace_args(u0, tspan, p)
    isinplace = !(u0 isa StaticArray)
    _f = if f === prob.f && (
            prob.f.f === trace_canonical! || prob.f.f === trace_canonical ||
                prob.f.f === trace_canonical_relativistic! || prob.f.f === trace_canonical_relativistic
        )
        is_rel = prob.f.f === trace_canonical_relativistic! ||
            prob.f.f === trace_canonical_relativistic
        if is_rel
            isinplace ? ODEFunction{true, DEFAULT_SPECIALIZATION}(trace_canonical_relativistic!) :
                ODEFunction{false, DEFAULT_SPECIALIZATION}(trace_canonical_relativistic)
        else
            isinplace ? ODEFunction{true, DEFAULT_SPECIALIZATION}(trace_canonical!) :
                ODEFunction{false, DEFAULT_SPECIALIZATION}(trace_canonical)
        end
    else
        f
    end
    return TraceCanonicalProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function solve(
        prob::TraceCanonicalProblem, alg::SciMLBase.AbstractODEAlgorithm,
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

# TraceProblem constructors and SciML interface integration

const FO_PARAMS = (:q2m, :m, :E, :B, :F)
const GC_PARAMS = (:q, :q2m, :μ, :E, :B)
const CANONICAL_PARAMS = (:q, :m, :c, :pf)

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false

"""
    _named_params(p, names)

Name the slots of a positional parameter tuple `p` after `names`, so that tuples
written in the documented layout keep working. Other containers pass through.
"""
function _named_params(p::NTuple{N, Any}, names) where {N}
    N <= length(names) || throw(
        ArgumentError("expected at most $(length(names)) parameters $names, got $N")
    )
    return NamedTuple{ntuple(i -> names[i], Val(N))}(p)
end
_named_params(p, names) = p

@inline _promote_p(::Type{T}, p::NamedTuple) where {T <: AbstractFloat} =
    map(x -> x isa Real ? T(x) : x, p)
@inline _promote_p(::Type{T}, p) where {T} = p

@inline function _promote_trace_args(u0, tspan, p, names)
    p = _named_params(p, names)
    T = eltype(u0)
    if T <: AbstractFloat
        tspan_T = (T(tspan[1]), T(tspan[2]))
        p_T = _promote_p(T, p)
        return tspan_T, p_T
    else
        return tspan, p
    end
end

_fo_equations(isinplace, relativistic) = relativistic ?
    (isinplace ? trace_relativistic! : trace_relativistic) :
    (isinplace ? trace! : trace)

"""
    TraceProblem(u0, tspan, p; relativistic = false, prob_func = DEFAULT_PROB_FUNC)
    TraceProblem(f, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)

Problem specification for particle tracing. Supports both Boris solvers and standard ODE
solvers (e.g. Tsit5, Vern7). When `f` is not provided, defaults to `trace!` (in-place)
for mutable arrays or `trace` (out-of-place) for `StaticArray`s, and to
`trace_relativistic!`/`trace_relativistic` with `relativistic = true`, in which case `u0`
holds `γv` in place of `v`. Boris solvers are non-relativistic.

Normalized units need no separate equations: prepare `p` with `q = 1, m = 1`, and
`c = 1` for the relativistic case.
"""
function TraceProblem(
        u0, tspan, p; relativistic::Bool = false, prob_func = DEFAULT_PROB_FUNC
    )
    tspan, p = _promote_trace_args(u0, tspan, p, FO_PARAMS)
    isinplace = !(u0 isa StaticArray)
    _func = _fo_equations(isinplace, relativistic)
    _f = ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(_func)
    return TraceProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function TraceProblem(f::Function, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    tspan, p = _promote_trace_args(u0, tspan, p, FO_PARAMS)
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
    tspan, p = _promote_trace_args(u0, tspan, p, FO_PARAMS)
    return TraceProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function TraceProblem{iip}(; f, u0, tspan, p, prob_func = DEFAULT_PROB_FUNC) where {iip}
    tspan, p = _promote_trace_args(u0, tspan, p, FO_PARAMS)
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
    tspan, p = _promote_trace_args(u0, tspan, p, FO_PARAMS)
    isinplace = !(u0 isa StaticArray)
    relativistic = _is_relativistic(prob.f.f)
    _f = if f === prob.f && (relativistic || prob.f.f in (trace!, trace))
        ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(
            _fo_equations(isinplace, relativistic)
        )
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

_canonical_equations(isinplace, relativistic) = relativistic ?
    (isinplace ? trace_canonical_relativistic! : trace_canonical_relativistic) :
    (isinplace ? trace_canonical! : trace_canonical)

function TraceCanonicalProblem(
        u0, tspan, p;
        relativistic::Bool = false,
        prob_func = DEFAULT_PROB_FUNC
    )
    tspan, p = _promote_trace_args(u0, tspan, p, CANONICAL_PARAMS)
    isinplace = !(u0 isa StaticArray)
    _func = _canonical_equations(isinplace, relativistic)
    _f = ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(_func)
    return TraceCanonicalProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function TraceCanonicalProblem(f::Function, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    tspan, p = _promote_trace_args(u0, tspan, p, CANONICAL_PARAMS)
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
    tspan, p = _promote_trace_args(u0, tspan, p, CANONICAL_PARAMS)
    return TraceCanonicalProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function TraceCanonicalProblem{iip}(; f, u0, tspan, p, prob_func = DEFAULT_PROB_FUNC) where {iip}
    tspan, p = _promote_trace_args(u0, tspan, p, CANONICAL_PARAMS)
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
    tspan, p = _promote_trace_args(u0, tspan, p, CANONICAL_PARAMS)
    isinplace = !(u0 isa StaticArray)
    relativistic = prob.f.f in (trace_canonical_relativistic!, trace_canonical_relativistic)
    _f = if f === prob.f && (relativistic || prob.f.f in (trace_canonical!, trace_canonical))
        ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(
            _canonical_equations(isinplace, relativistic)
        )
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

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

# Only the particle constants follow the state's float type; other entries of `p`, such as
# user data, keep theirs.
const PARTICLE_CONSTANTS = (:q, :q2m, :m, :c, :μ)

@inline _promote_p(::Type{T}, p::NamedTuple{K}) where {T <: AbstractFloat, K} = NamedTuple{K}(
    map((k, x) -> k in PARTICLE_CONSTANTS && x isa Real ? T(x) : x, K, values(p))
)
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

_param_names(::Type{<:Union{FullOrbit, RelativisticOrbit}}) = FO_PARAMS
_param_names(::Type{<:Union{Canonical, RelativisticCanonical}}) = CANONICAL_PARAMS
_param_names(::Type{GuidingCenter}) = GC_PARAMS

# Default equations of each model. `_equations_model` is the reverse lookup, `nothing` for any
# other function.
_equations_model(f) = nothing
_equations_model(f::AbstractODEFunction) = _equations_model(SciMLBase.unwrapped_f(f.f))
for (M, iip, oop) in (
        (FullOrbit, trace!, trace),
        (RelativisticOrbit, trace_relativistic!, trace_relativistic),
        (Canonical, trace_canonical!, trace_canonical),
        (RelativisticCanonical, trace_canonical_relativistic!, trace_canonical_relativistic),
        (GuidingCenter, trace_gc!, trace_gc),
    )
    @eval begin
        _equations(::Type{$M}) = ($iip, $oop)
        _equations_model(::Union{$(typeof(iip)), $(typeof(oop))}) = $M
    end
end

_default_f(::Type{M}, isinplace) where {M} = ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(
    _equations(M)[isinplace ? 1 : 2]
)

"""
    _param_model(p)

Model that the parameters `p` are written for: the default where neither a problem nor the
equations name the model.
"""
@inline _param_model(p) = hasproperty(p, :μ) ? GuidingCenter() : FullOrbit()

"""
    _model(prob)

Model of `prob`. A plain `ODEProblem` has none of its own, so it follows from the equations
and the parameters.
"""
_model(::TraceProblem{M}) where {M} = M()
function _model(prob)
    M = _equations_model(prob.f)
    return isnothing(M) ? _param_model(prob.p) : M()
end

"""
    TraceProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    TraceProblem(f, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    TraceProblem{M}(...)

Problem specification for particle tracing in the model `M`: `FullOrbit` (the default),
`RelativisticOrbit`, `Canonical`, `RelativisticCanonical` or `GuidingCenter`, also written
`TraceCanonicalProblem` and `TraceGCProblem`. Supports the native solvers of each model and
standard ODE solvers (e.g. Tsit5, Vern7). `RelativisticOrbit` states hold `γv` in place of
`v`. Boris and hybrid solvers accept only the default `FullOrbit` equations without `F`.

When `f` is not provided, the model's equations are used: in-place for mutable arrays,
out-of-place for `StaticArray`s. The default equations of a model imply it, and `TraceProblem{M}`
throws on those of another. Custom `f` keeps `M`.

Normalized units need no separate equations: prepare `p` with `q = 1, m = 1`, and
`c = 1` for the relativistic case.
"""
TraceProblem(u0, tspan, p; kw...) = TraceProblem{FullOrbit}(u0, tspan, p; kw...)
TraceProblem(f, u0, tspan, p; kw...) =
    TraceProblem{something(_equations_model(f), FullOrbit)}(f, u0, tspan, p; kw...)

function TraceProblem{M}(
        u0, tspan, p; relativistic::Bool = false, prob_func = DEFAULT_PROB_FUNC
    ) where {M}
    N = relativistic ? _relativistic(M) : M
    return TraceProblem{N}(_default_f(N, !(u0 isa StaticArray)), u0, tspan, p; prob_func)
end

function TraceProblem{M}(f, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC) where {M}
    _f = ODEFunction{!(u0 isa StaticArray), DEFAULT_SPECIALIZATION}(f)
    return TraceProblem{M}(_f, u0, tspan, p; prob_func)
end

function TraceProblem{M}(
        f::AbstractODEFunction{iip}, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC
    ) where {M, iip}
    _equations_model(f) in (nothing, M) ||
        throw(ArgumentError("The equations of $(_equations_model(f)) do not fit the model $M."))
    tspan, p = _promote_trace_args(u0, tspan, p, _param_names(M))
    return TraceProblem{
        M, typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function SciMLBase.remake(
        prob::TraceProblem{M}; f = prob.f, u0 = prob.u0, tspan = prob.tspan, p = prob.p,
        prob_func = prob.prob_func, kwargs...
    ) where {M}
    # Default equations follow `u0` between mutable and static arrays.
    if f === prob.f && _equations_model(f) === M
        f = _default_f(M, !(u0 isa StaticArray))
    end
    return TraceProblem{something(_equations_model(f), M)}(f, u0, tspan, p; prob_func)
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

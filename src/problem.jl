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

# Default equations of each model: `(in-place, out-of-place)` pairs, the non-relativistic
# first and the relativistic, if any, second.
_equations(::Type{FullOrbit}) = ((trace!, trace), (trace_relativistic!, trace_relativistic))
_equations(::Type{Canonical}) = (
    (trace_canonical!, trace_canonical),
    (trace_canonical_relativistic!, trace_canonical_relativistic),
)
_equations(::Type{GuidingCenter}) = ((trace_gc!, trace_gc),)

_param_names(::Type{FullOrbit}) = FO_PARAMS
_param_names(::Type{Canonical}) = CANONICAL_PARAMS
_param_names(::Type{GuidingCenter}) = GC_PARAMS

function _default_f(::Type{M}, isinplace, relativistic) where {M}
    variants = _equations(M)
    relativistic && length(variants) < 2 &&
        throw(ArgumentError("$M has no relativistic equations."))
    f = variants[relativistic + 1][isinplace ? 1 : 2]
    return ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(f)
end

"""
    TraceProblem(u0, tspan, p; relativistic = false, prob_func = DEFAULT_PROB_FUNC)
    TraceProblem(f, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    TraceProblem{M}(...)

Problem specification for particle tracing in the model `M`: `FullOrbit` (the default),
`Canonical` or `GuidingCenter`, also written `TraceCanonicalProblem` and
`TraceGCProblem`. Supports the native solvers of each model and standard ODE solvers
(e.g. Tsit5, Vern7).

When `f` is not provided, the model's equations are used: in-place for mutable arrays,
out-of-place for `StaticArray`s, and relativistic with `relativistic = true`, in which
case `u0` holds `γv` in place of `v`. Boris solvers are non-relativistic.

Normalized units need no separate equations: prepare `p` with `q = 1, m = 1`, and
`c = 1` for the relativistic case.
"""
TraceProblem(args...; kw...) = TraceProblem{FullOrbit}(args...; kw...)

function TraceProblem{M}(
        u0, tspan, p; relativistic::Bool = false, prob_func = DEFAULT_PROB_FUNC
    ) where {M}
    f = _default_f(M, !(u0 isa StaticArray), relativistic)
    return TraceProblem{M}(f, u0, tspan, p; prob_func)
end

function TraceProblem{M}(f, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC) where {M}
    _f = ODEFunction{!(u0 isa StaticArray), DEFAULT_SPECIALIZATION}(f)
    return TraceProblem{M}(_f, u0, tspan, p; prob_func)
end

function TraceProblem{M}(
        f::AbstractODEFunction{iip}, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC
    ) where {M, iip}
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
    if f === prob.f
        i = findfirst(fs -> prob.f.f in fs, _equations(M))
        isnothing(i) || (f = _default_f(M, !(u0 isa StaticArray), i == 2))
    end
    return TraceProblem{M}(f, u0, tspan, p; prob_func)
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

const _GC_ODE_FUNC = ODEFunction{false, FullSpecialize}(trace_gc)

@inline function _gc_ode_problem(
        u0::SVector{4, T}, tspan::Tuple{T, T}, p,
        f::AbstractODEFunction = _GC_ODE_FUNC,
        prob_func = DEFAULT_PROB_FUNC
    ) where {T}
    return TraceProblem{
        GuidingCenter, SVector{4, T}, Tuple{T, T}, false,
        typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function _gc_initial_dt(prob::TraceGCProblem, ::Type{T}) where {T}
    q2m = prob.p.q2m
    Bmag = norm(prob.p.B(get_x(prob.u0), prob.tspan[1]))
    if Bmag == 0 || q2m == 0
        return T(1.0e-6)
    end
    tdir = sign(prob.tspan[2] - prob.tspan[1])
    tdir == 0 && (tdir = one(tdir))
    return T(tdir * 0.5 * 2π / abs(q2m * Bmag))
end

function _isoutside_callback(isoutside::F) where {F}
    isoutside === ODE_DEFAULT_ISOUTOFDOMAIN && return nothing

    condition = (u, t, integrator) -> isoutside(u, integrator.p, t)
    affect! = function (integrator)
        integrator.u = integrator.uprev
        return terminate!(integrator)
    end

    return DiscreteCallback(condition, affect!)
end

function solve(
        prob::TraceGCProblem, alg::SciMLBase.AbstractODEAlgorithm;
        dt::Union{Nothing, Real} = nothing,
        adaptive::Bool = isadaptive(alg),
        dense::Bool = false,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        callback = nothing,
        kwargs...
    ) where {F}
    T = eltype(prob.u0)
    if dt === nothing
        adaptive || throw(
            ArgumentError("$(nameof(typeof(alg))) needs a time step `dt`.")
        )
        step = _gc_initial_dt(prob, T)
    else
        ttotal = prob.tspan[2] - prob.tspan[1]
        if !adaptive && sign(dt) != sign(ttotal) && ttotal != 0
            throw(
                ArgumentError(
                    "dt must have the same sign as tspan[2] - tspan[1] for " *
                        "$(nameof(typeof(alg)))."
                )
            )
        end
        step = T(dt)
    end

    tspan = (T(prob.tspan[1]), T(prob.tspan[2]))
    p = _promote_p(T, prob.p)
    f = prob.f.f in (trace_gc!, trace_gc) ? _GC_ODE_FUNC :
        ODEFunction{false, FullSpecialize}(prob.f.f)
    ode_prob = _gc_ode_problem(SVector{4, T}(prob.u0), tspan, p, f, prob.prob_func)

    cb_out = _isoutside_callback(isoutside)
    cb = cb_out === nothing ? callback :
        callback === nothing ? cb_out : SciMLBase.CallbackSet(cb_out, callback)

    return SciMLBase.__solve(
        ode_prob, alg; dt = step, adaptive, dense, callback = cb, kwargs...
    )
end

# TraceProblem constructors and SciML interface integration

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false

# Only the particle constants follow the state's float type; other entries of `p`, such as
# user data, keep theirs.
const PARTICLE_CONSTANTS = (:q, :q2m, :m, :c, :μ)

@inline _promote_p(::Type{T}, p::NamedTuple{K}) where {T <: AbstractFloat, K} = NamedTuple{K}(
    map((k, x) -> k in PARTICLE_CONSTANTS && x isa Real ? T(x) : x, K, values(p))
)
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
        (RelativisticGuidingCenter, trace_gc_relativistic!, trace_gc_relativistic),
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
`RelativisticOrbit`, `Canonical`, `RelativisticCanonical`, `GuidingCenter` or
`RelativisticGuidingCenter`, also written `TraceCanonicalProblem` and `TraceGCProblem`.
Supports the native solvers of each model and standard ODE solvers (e.g. Tsit5, Vern7).
`RelativisticOrbit` states hold `γv` in place of `v`, and `RelativisticGuidingCenter`
states hold `γv∥` in place of `v∥`. Boris solvers accept the default `FullOrbit` and
`RelativisticOrbit` equations without `F`; `AdaptiveHybrid` accepts `FullOrbit` without `F`.

When `f` is not provided, the model's equations are used: in-place for mutable arrays,
out-of-place for `StaticArray`s. The default equations of a model imply it, and `TraceProblem{M}`
throws on those of another. Custom `f` keeps `M`.

Normalized units need no separate equations: prepare `p` with `q = 1, m = 1`, and
`c = 1` for the relativistic case.

# Examples
```jldoctest
julia> B(x) = SA[0.0, 0.0, 1.0e-8];

julia> p = prepare(; B, species = Electron);

julia> prob = TraceProblem(SA[0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], (0.0, 1.0e-3), p);

julia> solve(prob, Vern7()).retcode == TestParticle.ReturnCode.Success
true

julia> v0 = SA[0.6 * TestParticle.c, 0.0, 0.0];

julia> γ = inv(sqrt(1 - sum(abs2, v0) / TestParticle.c^2));

julia> u0_rel = SA[0.0, 0.0, 0.0, (γ * v0)...];

julia> prob_rel = TraceProblem{RelativisticOrbit}(u0_rel, (0.0, 1.0e-3), p);

julia> sol_rel = solve(prob_rel, Vern7());

julia> isapprox(sum(abs2, sol_rel.u[end][4:6]), sum(abs2, γ * v0); rtol = 1.0e-6)
true
```
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
    tspan, p = _promote_trace_args(u0, tspan, p)
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

const _GCModel = Union{GuidingCenter, RelativisticGuidingCenter}
const _GC_ODE_FUNC = ODEFunction{false, FullSpecialize}(trace_gc)
const _GC_REL_ODE_FUNC = ODEFunction{false, FullSpecialize}(trace_gc_relativistic)

_gc_default_func(::Type{GuidingCenter}) = _GC_ODE_FUNC
_gc_default_func(::Type{RelativisticGuidingCenter}) = _GC_REL_ODE_FUNC

@inline function _gc_ode_problem(
        ::Type{M}, u0::SVector{4, T}, tspan::Tuple{T, T}, p,
        f::AbstractODEFunction = _gc_default_func(M),
        prob_func = DEFAULT_PROB_FUNC
    ) where {M <: _GCModel, T}
    return TraceProblem{
        M, SVector{4, T}, Tuple{T, T}, false,
        typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

@inline _gc_ode_problem(
    u0::SVector{4, T}, tspan::Tuple{T, T}, p,
    f::AbstractODEFunction = _GC_ODE_FUNC,
    prob_func = DEFAULT_PROB_FUNC
) where {T} = _gc_ode_problem(GuidingCenter, u0, tspan, p, f, prob_func)

function _gc_initial_dt(prob::TraceProblem{M}, ::Type{T}) where {M <: _GCModel, T}
    q2m = prob.p.q2m
    Bmag = norm(prob.p.B(get_x(prob.u0), prob.tspan[1]))
    if Bmag == 0 || q2m == 0
        return T(1.0e-6)
    end
    tdir = sign(prob.tspan[2] - prob.tspan[1])
    tdir == 0 && (tdir = one(tdir))
    γ = _gc_gamma(M(), prob.u0[4], Bmag, prob.p)
    return T(tdir * 0.5 * 2π * γ / abs(q2m * Bmag))
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
        prob::TraceProblem{M}, alg::SciMLBase.AbstractODEAlgorithm;
        dt::Union{Nothing, Real} = nothing,
        adaptive::Bool = isadaptive(alg),
        dense::Bool = false,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        callback = nothing,
        kwargs...
    ) where {M <: _GCModel, F}
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

    default_f = _gc_default_func(M)
    default_eqs = _equations(M)
    ode_prob = if prob.u0 isa SVector{4, T} && prob.tspan isa Tuple{T, T} &&
            prob.f === default_f
        prob
    else
        tspan = prob.tspan isa Tuple{T, T} ? prob.tspan :
            (T(prob.tspan[1]), T(prob.tspan[2]))
        p = prob.tspan isa Tuple{T, T} ? prob.p : _promote_p(T, prob.p)
        f = prob.f.f in default_eqs ? default_f :
            ODEFunction{false, FullSpecialize}(prob.f.f)
        _gc_ode_problem(M, SVector{4, T}(prob.u0), tspan, p, f, prob.prob_func)
    end

    cb_out = _isoutside_callback(isoutside)
    cb = cb_out === nothing ? callback :
        callback === nothing ? cb_out : SciMLBase.CallbackSet(cb_out, callback)

    return SciMLBase.__solve(
        ode_prob, alg; dt = step, adaptive, dense, callback = cb, kwargs...
    )
end

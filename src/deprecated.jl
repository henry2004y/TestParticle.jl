# Kept for compatibility until the next breaking release, then removed.

export get_work, TraceHybridProblem, GCRK4, GCRK45

struct GCRK4 <: SciMLBase.AbstractODEAlgorithm end

struct GCRK45{T} <: SciMLBase.AbstractODEAlgorithm
    abstol::T
    reltol::T
end

GCRK45(; abstol = 1.0e-6, reltol = 1.0e-6) = GCRK45(promote(abstol, reltol)...)

@noinline function _loworderrk_missing()
    throw(
        ArgumentError(
            "GCRK4 is deprecated and forwards to `RK4()` from OrdinaryDiffEqLowOrderRK. " *
                "Run `using OrdinaryDiffEqLowOrderRK` and pass `RK4()` directly."
        )
    )
end

function _rk4_alg(::Any)
    try
        return default_gc_alg()
    catch
        _loworderrk_missing()
    end
end
_rk4_alg() = _rk4_alg(Val(:default))

function solve(prob::TraceGCProblem, alg::GCRK4; kwargs...)
    Base.depwarn(
        "`GCRK4()` is deprecated, use `RK4()` from `OrdinaryDiffEqLowOrderRK` " *
            "with `adaptive = false` instead.",
        :GCRK4
    )
    return solve(prob, _rk4_alg(); adaptive = false, kwargs...)
end

function solve(prob::TraceGCProblem, alg::GCRK45; kwargs...)
    Base.depwarn(
        "`GCRK45()` is deprecated, use `Vern7()` from `OrdinaryDiffEqVerner` " *
            "(or another SciML ODE solver) instead.",
        :GCRK45
    )
    return solve(
        prob, default_gc_alg();
        abstol = alg.abstol, reltol = alg.reltol, kwargs...
    )
end

TraceHybridProblem(u0, tspan, p; kw...) = TraceProblem(u0, tspan, p; kw...)

PotentialField(phi, A; kw...) = PotentialField(; phi, A, kw...)
PotentialField(A; phi = nothing, kw...) = PotentialField(; phi, A, kw...)

prepare_canonical(phi, A; kw...) = prepare_canonical(; phi, A, kw...)
prepare_canonical(A; phi = nothing, kw...) = prepare_canonical(; phi, A, kw...)

get_work(sol::AbstractODESolution) = get_work_rates(sol)

get_work_rates(xv, p, t, magnetic_properties, E = nothing) =
    get_work_rates(xv, p, t; magnetic_properties, E)

function _relativistic(::Type{M}) where {M}
    N = M === FullOrbit ? RelativisticOrbit : M === Canonical ? RelativisticCanonical :
        throw(ArgumentError("$M has no relativistic equations."))
    Base.depwarn("`relativistic = true` is deprecated, use `TraceProblem{$(nameof(N))}`.", :TraceProblem)
    return N
end

function get_B_parameters(x, t, Bfunc)
    (; B, Bmag, b̂, ∇B, JB) = get_magnetic_properties(x, t, Bfunc)
    return B, Bmag, b̂, ∇B, JB
end

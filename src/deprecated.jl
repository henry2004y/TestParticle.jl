# Kept for compatibility until the next breaking release, then removed.

export GCRK4, GCRK45

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

function _relativistic(::Type{M}) where {M}
    N = M === FullOrbit ? RelativisticOrbit : M === Canonical ? RelativisticCanonical :
        throw(ArgumentError("$M has no relativistic equations."))
    Base.depwarn("`relativistic = true` is deprecated, use `TraceProblem{$(nameof(N))}`.", :TraceProblem)
    return N
end

Base.@deprecate_binding qᵢ Proton.q false ", use `Proton.q` instead."
Base.@deprecate_binding mᵢ Proton.m false ", use `Proton.m` instead."
Base.@deprecate_binding qₑ Electron.q false ", use `Electron.q` instead."
Base.@deprecate_binding mₑ Electron.m false ", use `Electron.m` instead."
Base.@deprecate_binding μ₀ MU0_SI false ", use `PlasmaPopulations.MU0_SI` instead."
Base.@deprecate_binding ϵ₀ EPS0_SI false ", use `PlasmaPopulations.EPS0_SI` instead."
Base.@deprecate_binding kB KB_SI false ", use `PlasmaPopulations.KB_SI` instead."

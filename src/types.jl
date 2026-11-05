abstract type AbstractField{itd} <: Function end

const DEFAULT_PROB_FUNC(prob, ctx) = prob

"""
    FullOrbit, RelativisticOrbit, Canonical, RelativisticCanonical,
    GuidingCenter, RelativisticGuidingCenter

Model of a [`TraceProblem`](@ref): the equations that its state `u0` and parameters `p`
follow, and thereby the solvers and diagnostics that apply to it. `RelativisticOrbit`
states hold `γv` in place of `v`, and `RelativisticGuidingCenter` states hold `γv∥` in
place of `v∥`.
"""
struct FullOrbit end
struct RelativisticOrbit end
struct Canonical end
struct RelativisticCanonical end
struct GuidingCenter end
struct RelativisticGuidingCenter end

struct TraceProblem{M, uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <:
    AbstractODEProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    "function for setting initial conditions"
    prob_func::PF
end

"""
    TraceCanonicalProblem(u0, tspan, p; prob_func)

[`TraceProblem`](@ref) in canonical coordinates: `u0 = [x..., p...]` and `p = (; q, m, c, pf)`
from [`prepare_canonical`](@ref). The relativistic variant is
`TraceProblem{RelativisticCanonical}`.
"""
const TraceCanonicalProblem = TraceProblem{Canonical}

"""
    TraceGCProblem(u0, tspan, p; prob_func)

[`TraceProblem`](@ref) of a guiding center with `p` from [`prepare`](@ref) without `F`.
`u0` is the guiding-center state `[X..., vpar, μ]` or a full-orbit state `[x..., v...]`,
converted by [`full_to_gc`](@ref) at `tspan[1]`; the same holds for `u0` set by `prob_func`.
The relativistic variant is `TraceProblem{RelativisticGuidingCenter}`, which takes
`[X..., γvpar, μ]` or `[x..., (γv)...]`.
"""
const TraceGCProblem = TraceProblem{GuidingCenter}

# SciML solvers read the problem-level solver keywords from `kwargs`; there are none.
Base.propertynames(::TraceProblem) = (:f, :u0, :tspan, :p, :prob_func, :kwargs)

function Base.getproperty(prob::TraceProblem, sym::Symbol)
    sym === :kwargs && return (;)
    return getfield(prob, sym)
end

"""
    PotentialField{TP, TA, TGP, TGA}

A representation of scalar potential `phi` and vector potential `A`,
along with optional spatial gradient functions `grad_phi` and `grad_A`.
"""
struct PotentialField{TP, TA, TGP, TGA}
    phi::TP
    A::TA
    grad_phi::TGP
    grad_A::TGA
end

# Meshes.jl grid types dummy stubs used at API boundary
abstract type CartesianGrid end
abstract type RectilinearGrid end
abstract type StructuredGrid end

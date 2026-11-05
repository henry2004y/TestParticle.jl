abstract type AbstractField{itd} <: Function end

"""
Type for the particles: `Proton`, `Electron`.
"""
struct Species{M, Q}
    m::M
    q::Q
end

Species{T}(m, q) where {T} = Species(T(m), T(q))
Species{T}(s::Species) where {T} = Species(T(s.m), T(s.q))
Base.convert(::Type{Species{T}}, s::Species) where {T} = Species{T}(s)


const DEFAULT_PROB_FUNC(prob, ctx) = prob

"""
    FullOrbit, RelativisticOrbit, Canonical, RelativisticCanonical, GuidingCenter

Model of a [`TraceProblem`](@ref): the equations that its state `u0` and parameters `p`
follow, and thereby the solvers and diagnostics that apply to it. `RelativisticOrbit`
states hold `γv` in place of `v`.
"""
struct FullOrbit end
struct RelativisticOrbit end
struct Canonical end
struct RelativisticCanonical end
struct GuidingCenter end

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

[`TraceProblem`](@ref) of a guiding center: `u0 = [X..., vpar]` and `p = (; q, q2m, μ, E, B)`
from [`prepare_gc`](@ref). Native solvers: [`GCRK4`](@ref), [`GCRK45`](@ref).
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

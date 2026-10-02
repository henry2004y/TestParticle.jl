abstract type AbstractField{itd} <: Function end

"""
Abstract type for tracing solutions.
"""
abstract type AbstractTraceSolution{T, N, S} <: AbstractODESolution{T, N, S} end

"""
Abstract type for velocity distribution functions.
"""
abstract type VDF end

"""
Type for the particles: `Proton`, `Electron`.
"""
struct Species{M, Q}
    m::M
    q::Q
end

const DEFAULT_PROB_FUNC(prob, ctx) = prob

struct TraceProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <:
    AbstractODEProblem{uType, tType, isinplace}
    f::F
    "initial condition"
    u0::uType
    "time span"
    tspan::tType
    "(q2m, m, E, B, F)"
    p::P
    "function for setting initial conditions"
    prob_func::PF
end

Base.propertynames(::TraceProblem) = (:f, :u0, :tspan, :p, :prob_func, :kwargs)

function Base.getproperty(prob::TraceProblem, sym::Symbol)
    sym === :kwargs && return (;)
    return getfield(prob, sym)
end

# Meshes.jl grid types dummy stubs used at API boundary
abstract type CartesianGrid end
abstract type RectilinearGrid end
abstract type StructuredGrid end

abstract type AbstractField{itd} <: Function end

"""
Abstract type for tracing solutions.
"""
abstract type AbstractTraceSolution{T, N, S} <: AbstractODESolution{T, N, S} end

abstract type AbstractBoris end

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

const DEFAULT_PROB_FUNC(prob, i, repeat) = prob

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

function TraceProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    _f = ODEFunction{true, DEFAULT_SPECIALIZATION}(x -> nothing) # dummy func
    return TraceProblem{
        typeof(u0), typeof(tspan), true, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end
# For remake
function TraceProblem{iip}(; f, u0, tspan, p, prob_func) where {iip}
    return TraceProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

struct Boris{Adaptive} <: AbstractBoris
    safety::Float64
end

Boris() = Boris{false}(0.0)

"""
    MultistepBoris{N}(; n=1)

The Multistep/Hyper Boris method of order `N`.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order:
2 (standard), 4, or 6 (Hyper-Boris).
"""
struct MultistepBoris{N} <: AbstractBoris
    n::Int
end

@inline function MultistepBoris{N}(; n::Int = 1) where {N}
    if N != 2 && N != 4 && N != 6
        throw(ArgumentError("Multistep Boris order N must be 2, 4, or 6."))
    end
    return MultistepBoris{N}(n)
end

const MultistepBoris2 = MultistepBoris{2}
const MultistepBoris4 = MultistepBoris{4}
const MultistepBoris6 = MultistepBoris{6}

const AdaptiveBoris = Boris{true}

"""
    AdaptiveBoris(; safety=0.1)

Adaptive Boris method with adaptive time stepping based on
local gyroperiod.
The time step is determined by
`dt = safety * T_gyro = safety * 2π / |qB/m|`.
"""
AdaptiveBoris(; safety = 0.1) = Boris{true}(Float64(safety))

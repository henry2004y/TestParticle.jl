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
struct MultistepBoris{N, Adaptive} <: AbstractBoris
    n::Int
    safety::Float64
end

@inline function MultistepBoris{N, Adaptive}(; n::Int = 1, safety = 0.1) where {N, Adaptive}
    if N != 2 && N != 4 && N != 6
        throw(ArgumentError("Multistep Boris order N must be 2, 4, or 6."))
    end
    return MultistepBoris{N, Adaptive}(n, Float64(safety))
end

@inline function MultistepBoris{N}(; n::Int = 1) where {N}
    return MultistepBoris{N, false}(; n = n, safety = 0.0)
end

const MultistepBoris2 = MultistepBoris{2, false}
const MultistepBoris4 = MultistepBoris{4, false}
const MultistepBoris6 = MultistepBoris{6, false}

const AdaptiveBoris = Boris{true}

"""
    AdaptiveBoris(; safety=0.1)

Adaptive Boris method with adaptive time stepping based on
local gyroperiod.
The time step is determined by
`dt = safety * T_gyro = safety * 2π / |qB/m|`.
"""
AdaptiveBoris(; safety = 0.1) = Boris{true}(Float64(safety))

const AdaptiveMultistepBoris{N} = MultistepBoris{N, true}

"""
    AdaptiveMultistepBoris{N}(; n=1, safety=0.1)

Adaptive Multistep/Hyper Boris method of order `N`.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order: 2, 4, or 6.
"""
AdaptiveMultistepBoris{N}(; n::Int = 1, safety = 0.1) where {N} =
    MultistepBoris{N, true}(n, Float64(safety))

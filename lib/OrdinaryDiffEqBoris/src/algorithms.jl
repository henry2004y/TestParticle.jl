"""
    Boris()

The standard Boris method for particle pushing in electric and magnetic fields.

This solver expects a problem where `p` is structured as `(q2m, m, E, B, ...)`, which matches the signature used by TestParticle.jl's `TraceProblem`.
"""
struct Boris <: OrdinaryDiffEqAlgorithm end

"""
    MultistepBoris{N}(; n=1)

The Multistep/Hyper Boris method of order `N`.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order. `N=2` corresponds to the Multicycle solver, while `N=4` or `N=6` are the Hyper Boris solvers.
"""
struct MultistepBoris{N} <: OrdinaryDiffEqAlgorithm
    n::Int
end
MultistepBoris{N}(; n::Int = 1) where {N} = MultistepBoris{N}(n)

"""
    MultistepBoris2(; n=1)

The Multicycle Boris method (MultistepBoris with N=2).
"""
const MultistepBoris2 = MultistepBoris{2}

"""
    MultistepBoris4(; n=1)

The 4th order Hyper Boris method (MultistepBoris with N=4).
"""
const MultistepBoris4 = MultistepBoris{4}

"""
    MultistepBoris6(; n=1)

The 6th order Hyper Boris method (MultistepBoris with N=6).
"""
const MultistepBoris6 = MultistepBoris{6}

"""
    AdaptiveBoris(; safety=0.1)

Adaptive Boris method with adaptive time stepping based on the local gyroperiod.
The time step is evaluated as `dt = safety * 2π / |qB/m|`.
"""
struct AdaptiveBoris{T} <: OrdinaryDiffEqAlgorithm
    safety::T
end
AdaptiveBoris(; safety = 0.1) = AdaptiveBoris(safety)

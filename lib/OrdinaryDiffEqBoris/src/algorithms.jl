"""
    Boris()

The standard Boris method for particle pushing in electric and magnetic fields.

This solver expects a problem where `p` is structured as `(q2m, m, E, B, ...)`, which matches the signature used by TestParticle.jl's `TraceProblem`.
"""
struct Boris <: OrdinaryDiffEqAlgorithm end

"""
    MultistepBoris(; n=1, N=2)

The Multistep/Hyper Boris method.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order. `N=2` corresponds to the Multicycle solver, while `N=4` or `N=6` is the Hyper Boris solver.
"""
struct MultistepBoris <: OrdinaryDiffEqAlgorithm
    n::Int
    N::Int
end
MultistepBoris(; n::Int = 1, N::Int = 2) = MultistepBoris(n, N)

"""
    AdaptiveBoris(; safety=0.1)

Adaptive Boris method with adaptive time stepping based on the local gyroperiod.
The time step is evaluated as `dt = safety * 2π / |qB/m|`.
"""
struct AdaptiveBoris{T} <: OrdinaryDiffEqAlgorithm
    safety::T
end
AdaptiveBoris(; safety = 0.1) = AdaptiveBoris(safety)

"""
    Boris()

The standard Boris method for particle pushing in electric and magnetic fields.

This solver expects a problem where `p` is structured as `(q2m, m, E, B, ...)`, which matches the signature used by TestParticle.jl's `TraceProblem`.
"""
struct Boris{T} <: OrdinaryDiffEqAlgorithm
    safety::T
end
Boris(; safety = 0.0) = Boris(safety)

"""
    MultistepBoris{N}(; n=1)

The Multistep/Hyper Boris method of order `N`.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order. `N=2` corresponds to the Multicycle solver, while `N=4` or `N=6` are the Hyper Boris solvers.
"""
struct MultistepBoris{N, T} <: OrdinaryDiffEqAlgorithm
    n::Int
    safety::T
end
MultistepBoris{N}(; n::Int = 1, safety = 0.0) where {N} = MultistepBoris{N, typeof(safety)}(n, safety)

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

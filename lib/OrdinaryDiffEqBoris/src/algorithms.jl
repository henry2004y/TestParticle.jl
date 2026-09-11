"""
    Boris()

The standard Boris method for particle pushing in electric and magnetic fields.

The charge-to-mass ratio and the field functions are read from the problem
parameter `p` through `get_q2m`, `get_EField` and `get_BField`.

The time step is fixed and must be supplied with the `dt` keyword.
"""
struct Boris <: OrdinaryDiffEqAlgorithm end

"""
    MultistepBoris{N}(; n=1)

The Multistep/Hyper Boris method of order `N`.
`n` specifies the number of subcycles.
`N` specifies the gyrophase correction order. `N=2` corresponds to the
Multicycle solver, while `N=4` or `N=6` are the Hyper Boris solvers.

The time step is fixed and must be supplied with the `dt` keyword.
"""
struct MultistepBoris{N} <: OrdinaryDiffEqAlgorithm
    n::Int
    function MultistepBoris{N}(n::Int) where {N}
        N in (2, 4, 6) || throw(ArgumentError("Multistep Boris order N must be 2, 4, or 6."))
        return new{N}(n)
    end
end
MultistepBoris{N}(; n::Int = 1) where {N} = MultistepBoris{N}(n)

"""
    AdaptiveBoris(; safety = 0.1)

The Boris method with a time step that tracks the local gyroperiod,

```math
\\Delta t = \\text{safety} \\, \\frac{2\\pi}{|q B / m|}.
```

No error estimate is formed, so every step is accepted and `safety` is a plain
fraction of the local gyroperiod rather than a tolerance. Adaptivity is
controlled by the `adaptive` solve keyword, which defaults to `true` here and
can be turned off to recover a fixed time step.
"""
struct AdaptiveBoris{T} <: OrdinaryDiffEqAdaptiveAlgorithm
    safety::T
end
AdaptiveBoris(; safety = 0.1) = AdaptiveBoris(safety)

"""
    AdaptiveMultistepBoris{N}(; n=1, safety=0.1)

The Multistep/Hyper Boris method of order `N` with a time step that tracks the
local gyroperiod, as in `AdaptiveBoris`.
"""
struct AdaptiveMultistepBoris{N, T} <: OrdinaryDiffEqAdaptiveAlgorithm
    n::Int
    safety::T
end
function AdaptiveMultistepBoris{N}(; n::Int = 1, safety = 0.1) where {N}
    N in (2, 4, 6) || throw(ArgumentError("Multistep Boris order N must be 2, 4, or 6."))
    return AdaptiveMultistepBoris{N, typeof(safety)}(n, safety)
end

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

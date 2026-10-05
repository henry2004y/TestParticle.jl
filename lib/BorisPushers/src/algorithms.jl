"""
    AbstractBoris

Abstract supertype for all Boris particle-pusher algorithms in this package:

- `Boris()`: standard Boris method, second order, fixed step.
- `MultistepBoris{N}(; n=1)`: `n` sub-cycles per step with gyrophase correction
  of order `N` (`N = 2, 4, 6`), fixed step.
- `AdaptiveBoris(; safety=0.1)`: Boris method with gyroperiod-tracking time step.
- `AdaptiveMultistepBoris{N}(; n=1, safety=0.1)`: `MultistepBoris{N}` with
  gyroperiod-tracking time step.
"""
abstract type AbstractBoris <: OrdinaryDiffEqAlgorithm end

"""
    Boris()

The standard Boris method for particle pushing in electric and magnetic fields.

The charge-to-mass ratio and the field functions are read from the problem
parameter `p` through `get_q2m`, `get_EField` and `get_BField`.

The time step is fixed and must be supplied with the `dt` keyword.
"""
struct Boris <: AbstractBoris end

"""
    MultistepBoris{N}(; n=1)

The Boris method with `n` subcycles per step, optionally corrected to order `N`
([Zenitani & Kato, 2025](https://arxiv.org/abs/2505.02270)).

`n` splits the rotation of one step into `n` smaller rotations while the fields
are still evaluated once per step. The gyration is then resolved more finely at
almost no cost, but the order of the method is unchanged.

`N` is the order to which the gyrophase is tracked:

  - `N = 2`: no correction. Second order multicycle method.
  - `N = 4`, `N = 6`: Hyper Boris methods with gyrophase error shrinking
    as `Δt⁴` and `Δt⁶`.

The time step is fixed and must be supplied with the `dt` keyword.
"""
struct MultistepBoris{N} <: AbstractBoris
    n::Int
    function MultistepBoris{N}(n::Int) where {N}
        N in (2, 4, 6) ||
            throw(ArgumentError("Multistep Boris order N must be 2, 4, or 6; got $N."))
        n >= 1 || throw(ArgumentError("Subcycle count n must be >= 1; got $n."))
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
struct AdaptiveBoris{T} <: AbstractBoris
    safety::T
    function AdaptiveBoris(safety::T) where {T}
        safety > 0 || throw(ArgumentError("safety must be positive; got $safety."))
        return new{T}(safety)
    end
end
AdaptiveBoris(; safety = 0.1) = AdaptiveBoris(safety)

"""
    AdaptiveMultistepBoris{N}(; n=1, safety=0.1)

`MultistepBoris{N}` with a time step that tracks the local gyroperiod, as in
`AdaptiveBoris`. See `MultistepBoris{N}` for `n` and `N`.
"""
struct AdaptiveMultistepBoris{N, T} <: AbstractBoris
    n::Int
    safety::T
    function AdaptiveMultistepBoris{N, T}(n::Int, safety::T) where {N, T}
        N in (2, 4, 6) ||
            throw(ArgumentError("Multistep Boris order N must be 2, 4, or 6; got $N."))
        n >= 1 || throw(ArgumentError("Subcycle count n must be >= 1; got $n."))
        safety > 0 || throw(ArgumentError("safety must be positive; got $safety."))
        return new{N, T}(n, safety)
    end
end
function AdaptiveMultistepBoris{N}(n::Int, safety::T) where {N, T}
    return AdaptiveMultistepBoris{N, T}(n, safety)
end
function AdaptiveMultistepBoris{N}(; n::Int = 1, safety = 0.1) where {N}
    return AdaptiveMultistepBoris{N, typeof(safety)}(n, safety)
end

"""
    MultistepBoris2(; n=1)

The multicycle Boris method, `MultistepBoris` with `N = 2`. Second order.
"""
const MultistepBoris2 = MultistepBoris{2}

"""
    MultistepBoris4(; n=1)

The fourth order Hyper Boris method, `MultistepBoris` with `N = 4`.
"""
const MultistepBoris4 = MultistepBoris{4}

"""
    MultistepBoris6(; n=1)

The sixth order Hyper Boris method, `MultistepBoris` with `N = 6`.
"""
const MultistepBoris6 = MultistepBoris{6}

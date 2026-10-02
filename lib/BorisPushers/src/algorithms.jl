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

The Boris method with `n` subcycles per step, optionally corrected to order `N`
([Zenitani & Kato, 2025](https://arxiv.org/abs/2505.02270)).

`n` splits the rotation of one step into `n` smaller rotations while the fields
are still evaluated once per step. The gyration is then resolved more finely at
almost no cost, but the order of the method is unchanged.

`N` is the order to which the gyrophase is tracked:

  - `N = 2`: no correction. This is the *multicycle* method, so called because
    the rotation cycle is repeated `n` times within one step. Second order.
  - `N = 4`, `N = 6`: the *Hyper Boris* methods, which correct the electric and
    magnetic terms with the factors `f_N` and `c_N` so that the gyrophase error
    shrinks as `Δt⁴` and `Δt⁶` instead of `Δt²`.

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

`MultistepBoris{N}` with a time step that tracks the local gyroperiod, as in
`AdaptiveBoris`. See `MultistepBoris{N}` for `n` and `N`.
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

"""
    AbstractBoris

Union of the Boris particle pushers defined in this package:

- `Boris()`: the standard Boris method, second order, fixed step.
- `AdaptiveBoris(; safety)`: the same method with a step that follows the local
  gyroperiod.
- `MultistepBoris{N}(; n)`: `n` sub-cycles per step with a gyrophase correction
  of order `N` (`N = 2` is the multicycle method, `N = 4` and `N = 6` are the
  Hyper Boris methods), fixed step.
- `AdaptiveMultistepBoris{N}(; n, safety)`: the same with a gyroperiod-following
  step.

They are ordinary SciML algorithms, so they are used like any other one,
`solve(prob, Boris(); dt)`, and are accepted inside a SciML `EnsembleProblem`.
Instead of indexing into the parameter container `p`, the methods ask it for the
charge-to-mass ratio and the two field functions through `get_q2m(p)`,
`get_EField(p)` and `get_BField(p)`, so any container that answers those three
can be traced.

The algorithms are immutable and hold nothing but numbers, so one can be passed
into a GPU kernel. Keeping it that way is a constraint: a field function or any
other non-isbits value in an algorithm would stop it from crossing to the device.
"""
const AbstractBoris =
    Union{Boris, AdaptiveBoris, MultistepBoris, AdaptiveMultistepBoris}

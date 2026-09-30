abstract type AbstractField{itd} <: Function end

"""
Abstract type for tracing solutions.
"""
abstract type AbstractTraceSolution{T, N, S} <: AbstractODESolution{T, N, S} end

"""
    AbstractBoris

Union of the Boris particle pushers, which are SciML algorithms defined in the
`OrdinaryDiffEqBoris` package and re-exported here:

- `Boris()`: the standard Boris method, second order, fixed step.
- `AdaptiveBoris(; safety)`: the same method with a step that follows the local
  gyroperiod.
- `MultistepBoris{N}(; n)`: `n` sub-cycles per step with a gyrophase correction
  of order `N` (`N = 2` is the multicycle method, `N = 4` and `N = 6` are the
  Hyper Boris methods), fixed step.
- `AdaptiveMultistepBoris{N}(; n, safety)`: the same with a gyroperiod-following
  step.

They are used like any other SciML algorithm, `solve(prob, Boris(); dt)`, and
are also accepted inside a SciML `EnsembleProblem`. Instead of indexing into the
parameter container `p`, the methods ask it for the charge-to-mass ratio and the
two field functions through `get_q2m(p)`, `get_EField(p)` and `get_BField(p)`,
so any container that answers those three can be traced.
"""
const AbstractBoris =
    Union{Boris, AdaptiveBoris, MultistepBoris, AdaptiveMultistepBoris}

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

# Meshes.jl grid types dummy stubs used at API boundary
abstract type CartesianGrid end
abstract type RectilinearGrid end
abstract type StructuredGrid end

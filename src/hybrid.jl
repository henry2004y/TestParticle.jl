# Adaptive hybrid GC-FO solver

const MIN_DT = 1.0e-14

"""
    AdaptiveHybrid{T}

Hybrid algorithm that switches between adaptive Boris (FO) and RK45 (GC).

Mode switching uses hysteresis with two adiabaticity thresholds to avoid
rapid GC ↔ FO chattering near a single bound:
- `threshold_gc_to_fo` (α): switch GC → FO when `ε ≥ α`.
- `threshold_fo_to_gc` (β): switch FO → GC when `ε < β`.
With `α ≥ β` the band `β ≤ ε < α` forms a buffer where the current mode is
retained, so a particle hovering near the boundary does not oscillate.

The adiabaticity parameter `ε` compared against the thresholds is selected by
the `adiabaticity` field, which may be:
- `:curvature` (default) → `ε_curv = ρ_L / R_c` (curvature drift; reproduces
  the legacy `get_adiabaticity` behavior exactly),
- `:gradB` → `ε_gradB = ρ_L / L_B`, with `L_B = |B| / |∇B|`,
- `:both` → OR of the two criteria: the solver switches to the full orbit
  whenever *either* `ε_curv ≥ α` *or* `ε_gradB ≥ α` (equivalently
  `max(ε_curv, ε_gradB) ≥ α`),
- `:jacobian` → `ε_jac = ρ_L · ‖JB‖_F / |B|`, CHIMP's original all-in-one
  criterion (Bates & Fuentes, 2020), where `JB` is the Jacobian of `B` and
  `‖JB‖_F` its Frobenius norm. It captures curvature, grad-B, torsion and
  shear at once, so it switches to full orbit in non-adiabatic regions that
  `:curvature`/`:gradB` miss (e.g. a field whose direction rotates in space
  while keeping constant magnitude).

Each check point records only the adiabaticity value used by the selected
`adiabaticity` mode (a single scalar per check) together with the active mode
in `sol.stats.adiabaticity`. This diagnostic is only saved when
`save_adiabaticity = true` (the default); set it to `false` to skip the
buffers and reduce memory/allocation overhead.
"""
struct AdaptiveHybrid{T}
    threshold_gc_to_fo::T
    threshold_fo_to_gc::T
    dtmin::T
    dtmax::T
    safety_fo::T
    abstol::T
    reltol::T
    maxiters::Int
    check_interval::Int
    adiabaticity::Symbol
    save_adiabaticity::Bool
end

function AdaptiveHybrid(;
        threshold = 0.1,
        threshold_gc_to_fo = threshold,
        threshold_fo_to_gc = threshold,
        dtmax, dtmin = 1.0e-2 * dtmax, safety_fo = 0.1,
        abstol = 1.0e-6, reltol = 1.0e-6, maxiters = 10000, check_interval = 10,
        adiabaticity = :curvature, save_adiabaticity = true
    )
    dtmax > 0 || throw(ArgumentError("dtmax must be positive."))
    dtmin > 0 || throw(ArgumentError("dtmin must be positive."))
    dtmin <= dtmax || throw(ArgumentError("dtmin must be <= dtmax."))
    safety_fo > 0 || throw(ArgumentError("safety_fo must be positive."))
    maxiters > 0 || throw(ArgumentError("maxiters must be positive."))
    threshold_fo_to_gc >= 0 ||
        throw(ArgumentError("threshold_fo_to_gc must be non-negative."))
    abstol > 0 || throw(ArgumentError("abstol must be positive."))
    reltol > 0 || throw(ArgumentError("reltol must be positive."))
    check_interval > 0 || throw(ArgumentError("check_interval must be positive."))
    threshold_gc_to_fo >= threshold_fo_to_gc ||
        throw(
        ArgumentError(
            "threshold_gc_to_fo must be >= threshold_fo_to_gc for hysteresis."
        )
    )
    adiabaticity in (:curvature, :gradB, :both, :jacobian) ||
        throw(
        ArgumentError(
            "adiabaticity must be one of :curvature, :gradB, :both, :jacobian."
        )
    )
    T = promote_type(
        typeof(threshold_gc_to_fo), typeof(threshold_fo_to_gc),
        typeof(dtmin), typeof(dtmax), typeof(safety_fo), typeof(abstol),
        typeof(reltol)
    )
    return AdaptiveHybrid{T}(
        T(threshold_gc_to_fo), T(threshold_fo_to_gc), T(dtmin), T(dtmax), T(safety_fo),
        T(abstol), T(reltol), maxiters, check_interval, adiabaticity, save_adiabaticity
    )
end

"""
    TraceHybridProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF}

Problem type for hybrid GC-FO tracing. Initial condition `u0` should be the full orbit state (6-element).
"""
struct TraceHybridProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <:
    AbstractODEProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    prob_func::PF
end

function TraceHybridProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    length(u0) == 6 ||
        throw(ArgumentError("u0 must have 6 elements (full orbit state)."))
    tspan[2] > tspan[1] ||
        throw(ArgumentError("tspan[2] must be greater than tspan[1]."))
    tspan, p = _promote_trace_args(u0, tspan, p, FO_PARAMS)
    _f = ODEFunction{true, DEFAULT_SPECIALIZATION}(x -> nothing) # dummy func
    return TraceHybridProblem{
        typeof(u0), typeof(tspan), true, typeof(p), typeof(_f),
        typeof(prob_func),
    }(
        _f,
        u0,
        tspan,
        p,
        prob_func
    )
end

# For remake
function TraceHybridProblem{iip}(; f, u0, tspan, p, prob_func) where {iip}
    p = _named_params(p, FO_PARAMS)
    return TraceHybridProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f),
        typeof(prob_func),
    }(
        f,
        u0,
        tspan,
        p,
        prob_func
    )
end

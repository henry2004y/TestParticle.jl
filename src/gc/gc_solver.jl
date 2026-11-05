# Native GC solver

"""
    GCRK4()

Fixed-step RK4 native solver of [`TraceGCProblem`](@ref); needs `dt`.
"""
struct GCRK4 <: SciMLBase.AbstractODEAlgorithm end

"""
    GCRK45(; abstol = 1e-6, reltol = 1e-6)

Adaptive RK45 native solver of [`TraceGCProblem`](@ref); `dt` is the optional initial step.
"""
struct GCRK45{T} <: SciMLBase.AbstractODEAlgorithm
    abstol::T
    reltol::T
end

GCRK45(; abstol = 1.0e-6, reltol = 1.0e-6) = GCRK45(promote(abstol, reltol)...)

const GCAlgorithm = Union{GCRK4, GCRK45}

"""
    solve(prob::TraceGCProblem, alg::Union{GCRK4, GCRK45};
        dt = nothing, saveat = (), isoutside = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start = true, save_end = true, save_everystep = true, maxiters = 10000)

Trace a guiding center with a native solver and return its `ODESolution`.

`saveat` names the times at which to report the state, as a collection or as an
interval; the state is interpolated linearly between the two steps that bracket
each time, leaving the trajectory unchanged.
"""
function solve(
        prob::TraceGCProblem, alg::GCAlgorithm;
        saveat = (),
        dt::Union{Real, Nothing} = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        maxiters = 10000
    ) where {F}
    plan = SavingPlan(saveat, prob.tspan, _span_direction(prob.tspan), _gc_time_type(prob, dt))

    return _trace_gc_native(
        alg, prob, dt, plan, isoutside, save_start, save_end, save_everystep, maxiters
    )
end

"""
    _gc_time_type(prob, dt) -> Type

The type of the time values a solve will produce.
"""
function _gc_time_type(prob::TraceGCProblem, dt)
    dt_guess = isnothing(dt) ? one(eltype(prob.u0)) : dt
    return typeof(prob.tspan[1] + dt_guess)
end

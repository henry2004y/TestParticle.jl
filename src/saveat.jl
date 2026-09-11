# `saveat` for the native solvers.
#
# `saveat` replaces the deprecated `savestepinterval`. Where that keyword picked
# every k-th step and so tied the output to the step grid, `saveat` names the
# times to report and the state is interpolated linearly between the two steps
# that bracket each of them. The integration is untouched, so the trajectory is
# the same as a run without `saveat`, exactly as for the Boris solvers.

"""
    SavingPlan(saveat, savestepinterval, tspan, dir, ::Type{T})

How a native solve chooses the times at which it reports the state: every
`interval`-th step when `interval` is positive, or an explicit list of `times`
when it is zero.

`saveat` and `savestepinterval` describe the output in incompatible ways, so
passing both is an error. `savestepinterval` is deprecated and warns once per
session.
"""
struct SavingPlan{T}
    interval::Int
    times::Vector{T}
    dir::Int
end

use_saveat(plan::SavingPlan) = plan.interval == 0

const _DEPRECATED_SAVESTEPINTERVAL =
    "The savestepinterval keyword is deprecated; use saveat instead."

function SavingPlan(saveat, savestepinterval, tspan, dir, ::Type{T}) where {T}
    if !isempty(saveat)
        isnothing(savestepinterval) || throw(
            ArgumentError(
                "saveat and savestepinterval select the output times in different " *
                    "ways and cannot be combined"
            )
        )
        return SavingPlan{T}(0, _saveat_times(saveat, tspan, dir, T), dir)
    end

    interval = 1
    if !isnothing(savestepinterval)
        # `maxlog` keeps the notice to one per session; the repository's own
        # tests exercise the keyword often enough that repeating it would bury
        # everything else.
        @warn _DEPRECATED_SAVESTEPINTERVAL maxlog = 1
        interval = savestepinterval
    end

    return SavingPlan{T}(interval, T[], dir)
end

"""
    _span_direction(tspan) -> Int

`+1` for a forward time span and `-1` for a backward one.
"""
_span_direction(tspan) = tspan[2] >= tspan[1] ? 1 : -1

"""
    _saveat_times(saveat, tspan, dir, ::Type{T}) -> Vector{T}

The interior times at which to report the state, ordered along the direction of
integration `dir`. Times outside `tspan` are dropped, and the two ends are left
out because `save_start` and `save_end` already report them.
"""
function _saveat_times(saveat, tspan, dir, ::Type{T}) where {T}
    t0, t1 = T(tspan[1]), T(tspan[2])

    raw = if saveat isa Number
        step = abs(T(saveat))
        step == 0 && throw(ArgumentError("saveat must not be zero"))
        n = floor(Int, abs(t1 - t0) / step)
        T[t0 + dir * k * step for k in 1:n]
    else
        T[T(t) for t in saveat]
    end

    times = T[t for t in raw if dir * (t - t0) > 0 && dir * (t - t1) < 0]
    return sort!(times; rev = dir < 0)
end

"""
    _saveat_reached(target, t, dir) -> Bool

Whether a solve that has advanced to `t` has reached `target`.
"""
@inline _saveat_reached(target, t, dir) = dir * (target - t) <= 0

"""
    _saveat_interpolate(t_prev, xv_prev, t_next, xv_next, t) -> xv

Linear interpolation of the state between two consecutive steps, which is the
dense output these solvers admit: it reproduces the step values exactly at the
step times and is as accurate as the method in between.
"""
@inline function _saveat_interpolate(t_prev, xv_prev, t_next, xv_next, t)
    θ = (t - t_prev) / (t_next - t_prev)
    return xv_prev + θ * (xv_next - xv_prev)
end

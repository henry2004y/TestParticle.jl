# `saveat` for the native solvers.
#
# A native solve reports every step by default. `saveat` names the times to
# report instead, and the state is interpolated linearly between the two steps
# that bracket each of them. The integration is untouched, so the trajectory is
# the same as a run without `saveat`, which is also how the SciML-backed Boris
# solvers behave.

"""
    SavingPlan(saveat, tspan, dir, ::Type{T})

The times at which a native solve reports the state, ordered along the direction
of integration `dir`. An empty plan means `saveat` was not used and every step is
reported.
"""
struct SavingPlan{T}
    times::Vector{T}
    dir::Int
end

use_saveat(plan::SavingPlan) = !isempty(plan.times)

function SavingPlan(saveat, tspan, dir, ::Type{T}) where {T}
    return SavingPlan{T}(_saveat_times(saveat, tspan, dir, T), dir)
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

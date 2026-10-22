struct SavingPlan{T}
    times::Vector{T}
    dir::Int
    has_saveat::Bool
end

use_saveat(plan::SavingPlan) = plan.has_saveat

function SavingPlan(saveat, tspan, dir, ::Type{T}) where {T}
    has_saveat = !isempty(saveat)
    return SavingPlan{T}(_saveat_times(saveat, tspan, dir, T), dir, has_saveat)
end

function SavingPlan(saveat, tspan, ::Type{T}) where {T}
    return SavingPlan(saveat, tspan, _span_direction(tspan), T)
end

_span_direction(tspan) = tspan[2] >= tspan[1] ? 1 : -1

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

@inline saveat_reached(target, t, dir) = dir * (target - t) <= 0

@inline function saveat_interpolate(t_prev, xv_prev, t_next, xv_next, t)
    θ = (t - t_prev) / (t_next - t_prev)
    return xv_prev + θ * (xv_next - xv_prev)
end

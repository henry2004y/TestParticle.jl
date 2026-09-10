mutable struct BorisConstantCache{vType, tType} <: OrdinaryDiffEqConstantCache
    v_half::vType
    dt_prev::tType
end

mutable struct BorisCache{uType, rateType, vType, tType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
    v_half::vType
    dt_prev::tType
end

mutable struct MultistepBorisConstantCache{vType, tType} <: OrdinaryDiffEqConstantCache
    v_half::vType
    dt_prev::tType
end

mutable struct MultistepBorisCache{uType, rateType, vType, tType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
    v_half::vType
    dt_prev::tType
end

@inline _empty_half_velocity(u) = zero(SVector(u[1], u[2], u[3]))

function alg_cache(
        alg::Union{Boris, AdaptiveBoris}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return BorisConstantCache(_empty_half_velocity(u), zero(dt))
end

function alg_cache(
        alg::Union{Boris, AdaptiveBoris}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return BorisCache(
        u, uprev, similar(u), similar(rate_prototype),
        _empty_half_velocity(u), zero(dt)
    )
end

function alg_cache(
        alg::Union{MultistepBoris{N}, AdaptiveMultistepBoris{N}}, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{false}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, N}
    return MultistepBorisConstantCache(_empty_half_velocity(u), zero(dt))
end

function alg_cache(
        alg::Union{MultistepBoris{N}, AdaptiveMultistepBoris{N}}, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{true}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, N}
    return MultistepBorisCache(
        u, uprev, similar(u), similar(rate_prototype),
        _empty_half_velocity(u), zero(dt)
    )
end

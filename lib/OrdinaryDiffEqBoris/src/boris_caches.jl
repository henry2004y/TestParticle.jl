struct BorisConstantCache <: OrdinaryDiffEqConstantCache end

struct BorisCache{uType, rateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
end

function OrdinaryDiffEqCore.alg_cache(
        alg::Boris, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return BorisConstantCache()
end

function OrdinaryDiffEqCore.alg_cache(
        alg::Boris, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return BorisCache(u, uprev, similar(u), similar(rate_prototype))
end

struct MultistepBorisConstantCache <: OrdinaryDiffEqConstantCache end

struct MultistepBorisCache{uType, rateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
end

function OrdinaryDiffEqCore.alg_cache(
        alg::MultistepBoris, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MultistepBorisConstantCache()
end

function OrdinaryDiffEqCore.alg_cache(
        alg::MultistepBoris, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MultistepBorisCache(u, uprev, similar(u), similar(rate_prototype))
end

struct AdaptiveBorisConstantCache <: OrdinaryDiffEqConstantCache end

struct AdaptiveBorisCache{uType, rateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
end

function OrdinaryDiffEqCore.alg_cache(
        alg::AdaptiveBoris, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return AdaptiveBorisConstantCache()
end

function OrdinaryDiffEqCore.alg_cache(
        alg::AdaptiveBoris, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return AdaptiveBorisCache(u, uprev, similar(u), similar(rate_prototype))
end

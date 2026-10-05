mutable struct BorisConstantCache{vType, tType, fType} <: OrdinaryDiffEqConstantCache
    v_half::vType
    dt_prev::tType
    fields::fType
end

mutable struct BorisCache{uType, rateType, vType, tType, fType} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
    v_half::vType
    dt_prev::tType
    fields::fType
end

@inline _position(u) = SVector(u[1], u[2], u[3])
@inline _empty_half_velocity(u) = zero(_position(u))

# A Boris step forms no derivative stages, so the mutable caches have no first
# and last stage to hand to the integrator, unlike the caches of a Runge-Kutta
# method.
get_fsalfirstlast(::BorisCache, u) = (nothing, nothing)

function alg_cache(
        alg::AbstractBoris, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt,
        reltol, p, calck, ::Val{false}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return BorisConstantCache(
        _empty_half_velocity(u), zero(dt), _node_fields(p, _position(u), t)
    )
end

function alg_cache(
        alg::AbstractBoris, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt,
        reltol, p, calck, ::Val{true}, args...; kwargs...
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return BorisCache(
        u, uprev, similar(u), similar(rate_prototype),
        _empty_half_velocity(u), zero(dt), _node_fields(p, _position(u), t)
    )
end

"""
    BorisController(safety)

Step size controller that tracks the local gyroperiod rather than a local error
estimate, proposing `dt = safety * 2π / |q B / m|`. Since no error estimate is
formed, steps are never rejected and `safety` is a fraction of the gyroperiod
rather than a tolerance.
"""
struct BorisController{QT} <: AbstractController
    safety::QT
end

mutable struct BorisControllerCache{QT, UT} <: AbstractControllerCache
    controller::BorisController{QT}
    q::QT
    atmp::UT
end

# The controller entry points were reworked between the OrdinaryDiffEqCore v3
# "legacy" controllers and the v7 controllers, which changed the signatures of
# `default_controller` and `setup_controller_cache`. The methods below are
# written to be valid under both layouts.
function default_controller(alg::Union{AdaptiveBoris, AdaptiveMultistepBoris}, cache, qoldinit, beta1, beta2)
    return BorisController(alg.safety)
end

function default_controller(::Type{QT}, alg::Union{AdaptiveBoris, AdaptiveMultistepBoris}) where {QT}
    return BorisController(QT(alg.safety))
end

function setup_controller_cache(alg, atmp, controller::BorisController{QT}, args...) where {QT}
    return BorisControllerCache{QT, typeof(atmp)}(controller, one(QT), atmp)
end

"""
    gyroperiod_dt(integrator, safety)

Return `safety * 2π / |q B / m|` evaluated at the position and time the last
step landed on, or `nothing` when the magnetic field vanishes there.
"""
function gyroperiod_dt(integrator, safety)
    q2m = get_q2m(integrator.p)
    Bfunc = get_BField(integrator.p)
    u = integrator.u
    Bmag = norm(Bfunc(SVector(u[1], u[2], u[3]), integrator.t + integrator.dt))
    iszero(Bmag) && return nothing

    return integrator.tdir * (2π * safety) / (abs(q2m) * Bmag)
end

function stepsize_controller!(integrator, cache::BorisControllerCache, alg)
    dt_proposed = gyroperiod_dt(integrator, cache.controller.safety)
    cache.q = dt_proposed === nothing ? one(cache.q) : integrator.dt / dt_proposed
    return cache.q
end

accept_step_controller(integrator, ::BorisControllerCache) = true

accept_step_controller(integrator, ::BorisControllerCache, alg) = true

function step_accept_controller!(integrator, cache::BorisControllerCache, alg, q)
    return integrator.dt / cache.q
end

function step_reject_controller!(integrator, ::BorisControllerCache, alg)
    return nothing
end

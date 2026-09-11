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

mutable struct BorisControllerCache{QT, EEstT, UT} <: AbstractControllerCache
    controller::BorisController{QT}
    q::QT
    EEst::EEstT
    atmp::UT
end

function default_controller(
        ::Type{QT}, alg::Union{AdaptiveBoris, AdaptiveMultistepBoris}
    ) where {QT}
    return BorisController(QT(alg.safety))
end

function setup_controller_cache(
        alg, atmp, controller::BorisController{QT}, ::Type{EEstT}, args...
    ) where {QT, EEstT}
    return BorisControllerCache{QT, EEstT, typeof(atmp)}(
        controller, one(QT), oneunit(EEstT), atmp
    )
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

# No error estimate is formed, so a step is never rejected.
accept_step_controller(integrator, ::BorisControllerCache, alg) = true

function step_accept_controller!(integrator, cache::BorisControllerCache, alg, q)
    return integrator.dt / cache.q
end

function step_reject_controller!(integrator, ::BorisControllerCache, alg)
    return nothing
end

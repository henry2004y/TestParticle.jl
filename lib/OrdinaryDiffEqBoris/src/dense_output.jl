# Dense output.
#
# A Boris step forms no derivative stages. The position advances linearly with
# the half-step velocity, `r_{n+1} = r_n + dt * v_{n+1/2}`, and the node
# velocity is reconstructed from that half-step velocity, so there are no
# endpoint derivatives to feed a Hermite interpolant. The dense output implied
# by the scheme is therefore linear between the endpoints of a step, which
# reproduces the stored position exactly and the velocity to the accuracy the
# scheme itself provides.
#
# OrdinaryDiffEqCore already falls back to linear interpolation whenever the
# stage array is empty, so the only thing needed here is to keep it empty. Left
# to the generic `_ode_addsteps!` it would instead try to fill `k[1]` and `k[2]`
# by evaluating `f`, which these methods do not define, and the solve would fail
# as soon as `saveat` asked for a value inside a step.

const BorisAlgorithm = Union{Boris, AdaptiveBoris, MultistepBoris, AdaptiveMultistepBoris}

const BorisCacheTypes = Union{
    BorisConstantCache, BorisCache,
    MultistepBorisConstantCache, MultistepBorisCache,
}

"""
    default_linear_interpolation(alg::BorisAlgorithm, prob)

Boris methods carry no derivative stages, so only linear interpolation is
available. Reporting that here keeps `dense` off, instead of storing a stage
history that is never filled.
"""
default_linear_interpolation(::BorisAlgorithm, prob) = true

@inline function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::BorisCacheTypes,
        always_calc_begin = false, allow_calc_end = true, force_calc_end = false
    )
    return nothing
end

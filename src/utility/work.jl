"""
    get_work_rates(xv, p, t; magnetic_properties = nothing, E = nothing) -> SVector{4}

Calculate the leading-order adiabatic work rates [W] `(P_par, P_fermi, P_grad, P_betatron)`
acting on a full-orbit particle or a guiding center:

| Term | Formula | Physics |
|:---|:---|:---|
| `P_par` | `q v∥ (E⋅b̂)` | direct acceleration by `E∥` |
| `P_fermi` | `(m v∥²/B) (b̂×κ)⋅E` | work by `E` on the curvature drift |
| `P_grad` | `(μ/B) (b̂×∇B)⋅E` | convective part of the betatron channel, `μ v_E⋅∇B` |
| `P_betatron` | `μ ∂B/∂t` | local part of the betatron channel |

Here `v_E = (E×b̂)/B` is the ExB drift, so `P_grad` is the perpendicular heating a
particle picks up by being convected across `∇B` at fixed `μ`, while `P_betatron` is
the heating from the local change of `|B|`. Together they are `μ (∂_t + v_E⋅∇) B`.
See [`get_work_channels`](@ref) for the grouping into physical acceleration channels.

The state is taken as a guiding center when `p` carries a magnetic moment `μ`. Saved
guiding-center states can carry diagnostics appended by `save_fields`/`save_work`, so
their length is not a reliable discriminator.

`magnetic_properties` and `E` are optional precomputed values of
[`get_magnetic_properties`](@ref) and of the electric field at `(get_x(xv), t)`, used
to avoid recomputing the field Jacobian when the caller already has them.

For a full-orbit state these are gyro-averaged quantities: they do not reproduce the
instantaneous `q v⋅E`, which oscillates at the gyrofrequency. They are therefore best
compared against the energy change accumulated over many gyroperiods.

For the guiding-center equations solved by [`trace_gc!`](@ref), the four terms sum to
`dK/dt` with `K = ½ m v∥² + μB` up to `O(ε)` corrections. The sum is *exact* when
`b̂⋅(∇×b̂) = 0`, as in axisymmetric fields without a `B_φ` component. In this case
`B*_∥ = |B|` and the residual twist term vanishes. Planar field lines or translation
invariance alone do not imply this condition. When it holds, any remaining mismatch
is numerical (e.g. time-integration or quadrature error).

All four terms are zero in a magnetic null.

Reference: Northrop, T. G. (1963), *The Adiabatic Motion of Charged Particles*.
"""
@inline function get_work_rates(xv, p, t; magnetic_properties = nothing, E = nothing)
    r = get_x(xv)
    (; ∇B, κ, b̂, Bmag) = @something magnetic_properties get_magnetic_properties(r, t, p.B)
    iszero(Bmag) && return zero(SVector{4, eltype(xv)})

    E = @something E p.E(r, t)
    q, m, v_par, μ = _work_rate_particle(xv, p, b̂, Bmag)
    dBdt = derivative_t(p.B, r, t) ⋅ b̂

    P_par = q * v_par * (E ⋅ b̂)
    P_fermi = (m * v_par^2 / Bmag) * ((b̂ × κ) ⋅ E)
    P_grad = (μ / Bmag) * ((b̂ × ∇B) ⋅ E)
    P_betatron = μ * dBdt

    return SVector{4}(P_par, P_fermi, P_grad, P_betatron)
end

# Charge, mass, parallel velocity and magnetic moment of the traced particle.
@inline function _work_rate_particle(xv, p, b̂, Bmag)
    if hasproperty(p, :μ)
        (; q, q2m, μ) = p
        return q, q / q2m, xv[4], μ
    end
    (; q2m, m) = p
    v = get_v(xv)
    v_perp = _get_v_perp(v, b̂)
    return q2m * m, m, v ⋅ b̂, m * (v_perp ⋅ v_perp) / (2 * Bmag)
end

"""
    get_work_channels(xv, p, t; kwargs...) -> SVector{3}

Group the work rates into the three physical acceleration channels
`(P_parallel, P_fermi, P_betatron)` [W].

The mapping between the four leading-order terms of [`get_work_rates`](@ref) and the
channels is:

- `P_parallel = P_par`: direct acceleration by `E∥`;
- `P_fermi = P_fermi`: work by `E` on the curvature drift;
- `P_betatron = P_betatron + P_grad`: perpendicular heating under conservation of `μ`.
  Grad-B drift work is *not* a separate mechanism; since
  `(μ/B)(b̂×∇B)⋅E = μ v_E⋅∇B`, it is the convective (spatial) part of the betatron
  channel, while `μ ∂B/∂t` is its local (temporal) part.

The state and the keywords are as in `get_work_rates`.
"""
function get_work_channels(xv, p, t; kwargs...)
    P_par, P_fermi, P_grad, P_betatron = get_work_rates(xv, p, t; kwargs...)

    return SVector{3}(P_par, P_fermi, P_grad + P_betatron)
end

"""
    get_work(sol::AbstractODESolution)

Return the work done by the electric field from the solution `sol`.
"""
function get_work(sol::AbstractODESolution)
    p = sol.prob.p
    return map((u, t) -> get_work_rates(u, p, t), sol.u, sol.t)
end

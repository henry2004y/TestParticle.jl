"""
    _work_rate_terms(q, m, v_par, μ, E, ∇B, κ, b̂, Bmag, dBdt) -> SVector{4}

Leading-order guiding-center work rates [W] shared by the full-orbit and the
guiding-center flavors. Returns `(P_par, P_fermi, P_grad, P_betatron)`:

| Term | Formula | Physics |
|:---|:---|:---|
| `P_par` | `q v∥ (E⋅b̂)` | direct acceleration by `E∥` |
| `P_fermi` | `(m v∥²/B) (b̂×κ)⋅E` | work by `E` on the curvature drift |
| `P_grad` | `(μ/B) (b̂×∇B)⋅E` | convective part of the betatron channel, `μ v_E⋅∇B` |
| `P_betatron` | `μ ∂B/∂t` | local part of the betatron channel |

Here `v_E = (E×b̂)/B` is the ExB drift, so `P_grad` is the perpendicular heating a
particle picks up by being convected across `∇B` at fixed `μ`, while `P_betatron` is
the heating from the local change of `|B|`. Together they are `μ (∂_t + v_E⋅∇) B`.

Reference: Northrop, T. G. (1963), *The Adiabatic Motion of Charged Particles*.
"""
@inline function _work_rate_terms(q, m, v_par::Number, μ, E, ∇B, κ, b̂, Bmag, dBdt)
    # 1. Parallel work: q v∥ (E ⋅ b̂)
    P_par = q * v_par * (E ⋅ b̂)
    # 2. Fermi work: m v∥² / B (b̂ × κ) ⋅ E
    P_fermi = (m * v_par^2 / Bmag) * ((b̂ × κ) ⋅ E)
    # 3. Gradient drift work: μ / B (b̂ × ∇B) ⋅ E = μ v_E ⋅ ∇B
    P_grad = (μ / Bmag) * ((b̂ × ∇B) ⋅ E)
    # 4. Betatron work: μ ∂B/∂t
    P_betatron = μ * dBdt

    return SVector{4}(P_par, P_fermi, P_grad, P_betatron)
end

"""
    get_work_rates(xv, p, t[, magnetic_properties, E_field]) -> SVector{4}

Calculate the work rates [W] acting on a full-orbit particle or a guiding center.
Returns the four leading-order adiabatic terms `(P_par, P_fermi, P_grad, P_betatron)`
documented in [`_work_rate_terms`](@ref); see [`get_work_channels`](@ref) for the
grouping into physical acceleration channels.

The state is taken as a guiding center when `p` carries a magnetic moment `μ`. Saved
guiding-center states can carry diagnostics appended by `save_fields`/`save_work`, so
their length is not a reliable discriminator.

`magnetic_properties` and `E_field` are optional precomputed values of
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
"""
@inline function get_work_rates(xv, p, t, magnetic_properties = nothing, E_field = nothing)
    return hasproperty(p, :μ) ?
        _get_work_rates_gc(xv, p, t, magnetic_properties, E_field) :
        _get_work_rates_fo(xv, p, t, magnetic_properties, E_field)
end

@inline function _get_work_rates_fo(xu, p, t, magnetic_properties = nothing, E_field = nothing)
    (; q2m, m) = p
    Efunc, Bfunc = p.E, p.B
    r = get_x(xu)
    q = q2m * m

    E = E_field === nothing ? Efunc(r, t) : E_field

    if magnetic_properties === nothing
        _, ∇B, κ, b̂, Bmag = get_magnetic_properties(r, t, Bfunc)
    else
        _, ∇B, κ, b̂, Bmag = magnetic_properties
    end

    if Bmag == 0
        return SVector{4, eltype(xu)}(0, 0, 0, 0)
    end

    v = get_v(xu)
    v_par = v ⋅ b̂
    v_perp = _get_v_perp(v, b̂)
    μ = m * (v_perp ⋅ v_perp) / (2 * Bmag)
    dBdt = derivative_t(Bfunc, r, t) ⋅ b̂

    return _work_rate_terms(q, m, v_par, μ, E, ∇B, κ, b̂, Bmag, dBdt)
end

@inline function _get_work_rates_gc(xv, p, t, magnetic_properties = nothing, E_field = nothing)
    (; q, q2m, μ) = p
    Efunc, Bfunc = p.E, p.B
    r = get_x(xv)
    v_par = xv[4]

    E = E_field === nothing ? Efunc(r, t) : E_field

    if magnetic_properties === nothing
        _, ∇B, κ, b̂, Bmag = get_magnetic_properties(r, t, Bfunc)
    else
        _, ∇B, κ, b̂, Bmag = magnetic_properties
    end

    if Bmag == 0
        return SVector{4, eltype(xv)}(0, 0, 0, 0)
    end

    m = q / q2m
    dBdt = derivative_t(Bfunc, r, t) ⋅ b̂

    return _work_rate_terms(q, m, v_par, μ, E, ∇B, κ, b̂, Bmag, dBdt)
end

"""
    get_work_channels(xv, p, t) -> SVector{3}

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

Both a full-orbit and a guiding-center state are accepted, as in `get_work_rates`.
"""
function get_work_channels(xv, p, t)
    P_par, P_fermi, P_grad, P_betatron = get_work_rates(xv, p, t)

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

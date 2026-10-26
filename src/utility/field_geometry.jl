@inline function _get_curvature(B, Bmag, b̂, JB)
    # Grad-B from Jacobian
    ∇B = JB' * b̂
    # Curvature vector κ = (b̂ ⋅ ∇) b̂
    return (JB * b̂ + b̂ * (-∇B ⋅ b̂)) / Bmag, ∇B
end

@inline function _get_B_jacobian(x, t, Bfunc)
    result = DiffResults.JacobianResult(x)
    result = ForwardDiff.jacobian!(result, r -> Bfunc(r, t), x)
    JB = DiffResults.jacobian(result)
    B = DiffResults.value(result)
    return B, JB
end

@inline function _get_B_jacobian(x, t, ::ZeroField)
    return zero(x), zeros(eltype(x), 3, 3)
end

"""
    get_magnetic_properties(x, t, Bfunc)

Calculate magnetic field properties at position `x` and time `t`.
Returns the NamedTuple `(; B, ∇B, κ, b̂, Bmag)`:
- `B`: Magnetic field vector
- `∇B`: Gradient of magnetic field magnitude
- `κ`: Curvature vector
- `b̂`: Unit magnetic field vector
- `Bmag`: Magnitude of B
"""
@inline function get_magnetic_properties(x, t, Bfunc)
    B, JB = _get_B_jacobian(x, t, Bfunc)

    Bmag = norm(B)
    if Bmag == 0
        vzero = zero(x)
        return (; B = vzero, ∇B = vzero, κ = vzero, b̂ = vzero, Bmag)
    end
    b̂ = B / Bmag

    κ, ∇B = _get_curvature(B, Bmag, b̂, JB)

    return (; B, ∇B, κ, b̂, Bmag)
end

"""
    get_curvature_radius(x, t, Bfunc)

Calculate the radius of curvature of the magnetic field at position `x` and time `t`.
Returns `Inf` if the field is zero or the field lines are straight.
"""
@inline function get_curvature_radius(x, t, Bfunc)
    (; κ) = get_magnetic_properties(x, t, Bfunc)

    k_mag = norm(κ)
    return iszero(k_mag) ? Inf : inv(k_mag)
end

"""
    get_adiabaticity(r, Bfunc, q, m, μ, t=0.0)
    get_adiabaticity(r, Bfunc, μ, t=0.0; species = Proton)

Calculate the adiabaticity parameter `ϵ = ρ / Rc` at position `r` and time `t`.
`ρ` is the gyroradius and `Rc` is the radius of curvature of the magnetic field.
"""
@inline function get_adiabaticity(r, Bfunc, q, m, μ, t = 0.0)
    B, JB = _get_B_jacobian(r, t, Bfunc)
    Bmag = norm(B)
    iszero(Bmag) && return Inf

    ρ = sqrt(2 * μ * m / Bmag) / abs(q) # Gyroradius

    κ, _ = _get_curvature(B, Bmag, B / Bmag, JB)

    k_mag = norm(κ)
    invRc = iszero(k_mag) ? zero(k_mag) : k_mag

    return ρ * invRc
end

function get_adiabaticity(r, Bfunc, μ, t::Number = 0.0; species = Proton)
    return get_adiabaticity(r, Bfunc, species.q, species.m, μ, t)
end

function get_adiabaticity(xv, p::NamedTuple, t)
    if hasproperty(p, :μ)
        (; q, q2m, μ) = p
        return get_adiabaticity(xv[SA[1, 2, 3]], p.B, q, q / q2m, μ, t)
    else
        (; q2m, m) = p
        Bfunc = p.B
        q = q2m * m

        x = xv[SA[1, 2, 3]]
        v = xv[SA[4, 5, 6]]

        B = Bfunc(x, t)
        Bmag = norm(B)
        b̂ = B / Bmag
        vpar_val = v ⋅ b̂
        vperp_vec = v - vpar_val * b̂
        μ = m * (vperp_vec ⋅ vperp_vec) / (2 * Bmag)

        return get_adiabaticity(x, Bfunc, q, m, μ, t)
    end
end

"""
    get_adiabaticity(sol::AbstractODESolution)

Calculate the adiabaticity parameter `ϵ` from a solution `sol` at all time steps.
The parameters `q`, `m`, `Bfunc` are extracted from `sol.prob.p`.
"""
function get_adiabaticity(sol::AbstractODESolution)
    n = length(sol.t)
    ε = Vector{Float64}(undef, n)
    p = sol.prob.p

    @inbounds for i in 1:n
        ε[i] = get_adiabaticity(sol.u[i], p, sol.t[i])
    end

    return ε
end

"""
    get_adiabaticity(sol::AbstractODESolution, t)

Calculate the adiabaticity parameter `ϵ` from a solution `sol` at time `t`.
`t` must be a scalar.
"""
function get_adiabaticity(sol::AbstractODESolution, t)
    return get_adiabaticity(sol(t), sol.prob.p, t)
end

"""
    adiabaticity_components(x, Bfunc, q, m, μ, t = 0.0)
    adiabaticity_components(x, Bfunc, μ, t = 0.0; species = Proton)

Return the tuple `(ε_curv, ε_gradB, ε_jac)` of first-order adiabaticity
parameters at position `x` and time `t`:
- `ε_curv  = ρ / R_c` (curvature drift; the definition used by
  [`get_adiabaticity`](@ref)),
- `ε_gradB = ρ / L_B` (grad-B drift; CHIMP's `eGC`), with
  `L_B = |B| / |∇B|` the magnetic-field length scale,
- `ε_jac   = ρ · ‖JB‖_F / |B|` (CHIMP's all-in-one criterion; `JB` is the
  Jacobian of `B` and `‖JB‖_F` its Frobenius norm). Because
  `‖JB‖_F² = |∇B|² + |B|²|κ|² + |B|²τ² + shear`, this folds curvature,
  grad-B, torsion and shear into one number, so it switches to full orbit in
  non-adiabatic regions that `ε_curv`/`ε_gradB` miss (e.g. a field whose
  direction rotates in space with constant magnitude). It always satisfies
  `ε_jac ≥ max(ε_curv, ε_gradB)`.

`ρ = √(2 μ m / |B|) / |q|` is the gyroradius. When `|B| = 0` all three
components are `Inf`.
"""
@inline function adiabaticity_components(x, Bfunc, q, m, μ, t = 0.0)
    B, JB = _get_B_jacobian(x, t, Bfunc)
    Bmag = norm(B)
    iszero(Bmag) && return (Inf, Inf, Inf)

    ρ = sqrt(2 * μ * m / Bmag) / abs(q)
    b̂ = B / Bmag

    # grad-B length scale L_B = |B| / |∇B|, with ∇B = ∇|B|
    ∇B = JB' * b̂
    gradB_mag = norm(∇B)
    L_B = iszero(gradB_mag) ? Inf : Bmag / gradB_mag
    ε_gradB = ρ / L_B

    # curvature drift: ε_curv = ρ · |κ|, κ = (b̂ ⋅ ∇) b̂
    κ, _ = _get_curvature(B, Bmag, b̂, JB)
    ε_curv = ρ * norm(κ)

    # CHIMP's full criterion: ε = ρ · ‖JB‖_F / |B|, where ‖JB‖_F is the
    # Frobenius norm of the B-Jacobian. This dominates both ε_curv and ε_gradB.
    ε_jac = ρ * sqrt(sum(abs2, JB)) / Bmag

    return (ε_curv, ε_gradB, ε_jac)
end

function adiabaticity_components(x, Bfunc, μ, t::Number = 0.0; species = Proton)
    return adiabaticity_components(x, Bfunc, species.q, species.m, μ, t)
end

"""
    get_B_parameters(x, t, Bfunc)

Evaluate the magnetic field `Bfunc` at `(x, t)` and return the field vector `B`,
its magnitude `Bmag`, the unit vector `b̂ = B / Bmag`, the gradient `∇B`, and the
Jacobian `JB` of the field.
"""
@inline function get_B_parameters(x, t, Bfunc)
    B, JB = _get_B_jacobian(x, t, Bfunc)
    Bmag = norm(B)
    b̂ = B / Bmag
    # Grad-B from Jacobian
    ∇B = JB' * b̂

    return B, Bmag, b̂, ∇B, JB
end

@inline function get_E_parameters(x, t, Efunc)
    E = Efunc(x, t)
    JE = jacobian(Efunc, x, t)

    return E, JE
end

# Guiding center.

"""
    _get_gc_parameters(xv, E, B, q, m, t) -> (X, vpar, q, q / m, μ, phase)

Extract the guiding-center position `X`, parallel velocity `vpar`, magnetic
moment `μ`, and gyrophase `phase` from a full-particle state `xv` at time `t`.
`E` and `B` must be callable as `E(x, t)` / `B(x, t)`; static fields simply
ignore `t`.

#TODO : support external forces
"""
function _get_gc_parameters(xv, E, B, q, m, t)
    x, v = xv[SA[1:3...]], xv[SA[4:6...]]

    bparticle = B(x, t)
    Bmag_particle = norm(bparticle)
    b̂particle = bparticle / Bmag_particle
    # vector of Larmor radius
    ρ = (b̂particle × v) / (q / m * Bmag_particle)
    # Get the guiding center location
    X = x - ρ
    # Get EM field at guiding center
    b = B(X, t)
    Bmag = norm(b)
    b̂ = b / Bmag
    vpar = b̂ ⋅ v

    vperp = v - vpar * b̂
    e = E(X, t)
    vE = (e × b̂) / Bmag
    w = vperp - vE
    μ = m * (w ⋅ w) / (2 * Bmag)

    e1, e2 = get_perp_vector(b̂)
    phase = atan(w ⋅ e2, w ⋅ e1)

    return X, vpar, q, q / m, μ, phase
end

"""
    prepare_gc(xv, xrange, yrange, zrange, E, B;
        species = Proton, q = nothing, m = nothing, order::Int = 1, bc = FillExtrap(NaN))
    prepare_gc(xv, E, B; species = Proton, q = nothing, m = nothing)

Prepare the guiding center parameters for a particle.
# TODO: support external forces
"""
function prepare_gc(
        xv, xrange::T, yrange::T, zrange::T, E::TE, B::TB;
        species = Proton, q = nothing, m = nothing, t = 0,
        order::Int = 1, bc = FillExtrap(NaN)
    ) where {T <: AbstractRange, TE, TB}
    q = @something q species.q
    m = @something m species.m

    E = TE <: AbstractArray ?
        build_interpolator(CartesianGrid, E, xrange, yrange, zrange, order, bc) :
        E
    B = TB <: AbstractArray ?
        build_interpolator(CartesianGrid, B, xrange, yrange, zrange, order, bc) :
        B
    fE = Field(E)
    fB = Field(B)

    X, vpar, q, q2m, μ, _ = _get_gc_parameters(xv, fE, fB, q, m, t)

    stateinit_gc = [X..., vpar]

    return stateinit_gc, (q, q2m, μ, fE, fB)
end

function prepare_gc(xv, E, B; species = Proton, q = nothing, m = nothing, t = 0)
    q = @something q species.q
    m = @something m species.m
    fE = Field(E)
    fB = Field(B)

    X, vpar, q, q2m, μ, _ = _get_gc_parameters(xv, fE, fB, q, m, t)

    stateinit_gc = [X..., vpar]

    return stateinit_gc, (q, q2m, μ, fE, fB)
end

"""
    full_to_gc(xu, param)

Convert full particle state `xu` to guiding center state `state_gc` and magnetic moment `μ`.
Returns `(state_gc, μ)`, where `state_gc = [R..., vpar]`.
"""
function full_to_gc(xu, param, t = 0)
    q2m = get_q2m(param)
    m = param[2]
    q = q2m * m
    E = get_EField(param)
    B = get_BField(param)

    X, vpar, _, _, μ, _ = _get_gc_parameters(xu, E, B, q, m, t)
    state_gc = SVector{4}(X[1], X[2], X[3], vpar)

    return state_gc, μ
end

"""
    gc_to_full(state_gc, param, μ, phase=0)

Convert guiding center state `state_gc` to full particle state `xu`.
Returns `xu = [x, y, z, vx, vy, vz]`.
"""
function gc_to_full(state_gc, param, μ, phase = 0; t = 0)
    q2m = get_q2m(param)
    m = param[2]
    q = q2m * m
    E_field = get_EField(param)
    B_field = get_BField(param)
    return Vector(_gc_to_full(state_gc, E_field, B_field, q, m, μ, t, phase))
end

"""
    _gc_to_full(state_gc, E_field, B_field, q, m, μ, t, phase = 0)

Internal GC → full-orbit conversion taking the field functions and charges
directly, so it can be shared by the GC and hybrid solvers. `E_field` and
`B_field` must be callable as `E_field(x, t)` / `B_field(x, t)`; static fields
simply ignore `t`.
"""
function _gc_to_full(state_gc, E_field, B_field, q, m, μ, t, phase = 0)
    R = get_x(state_gc)
    vpar = state_gc[4]

    E = E_field(R, t)
    B = B_field(R, t)
    Bmag = norm(B)
    b̂ = B / Bmag

    v_E = (E × b̂) / Bmag

    # perp speed, μ = m * w^2 / (2B)
    w = sqrt(2 * μ * Bmag / m)

    # perp vector
    e1, e2 = get_perp_vector(b̂)
    v_gyr = w * (cos(phase) * e1 + sin(phase) * e2)
    v_perp = v_gyr + v_E

    v = vpar * b̂ + v_perp

    # gyroradius
    Ω = q * Bmag / m
    ρ_vec = (b̂ × v_perp) / Ω

    x = R + ρ_vec

    return vcat(x, v)
end

"""
    get_gc(xu, param; use_vE = false)
    get_gc(x, y, z, vx, vy, vz, bx, by, bz, q2m)

Calculate the coordinates of the guiding center according to the phase space coordinates of a particle.
Reference: [wiki](https://en.wikipedia.org/wiki/Guiding_center)

Nonrelativistic definition:

```math
\\mathbf{X}=\\mathbf{x}-m\\frac{\\mathbf{b}\\times\\mathbf{v}}{qB}
```

If `use_vE=true`, it calculates the guiding center relative to the E × B drift frame:

```math
\\mathbf{X}=\\mathbf{x}-m\\frac{\\mathbf{b}\\times(\\mathbf{v}-\\mathbf{v}_E)}{qB}
```

which is essential for correctly identifying the polarization drift in time-varying fields.
"""
function get_gc(xu, param; use_vE::Bool = false)
    # TODO: support external forces
    q2m = get_q2m(param)
    B_field = get_BField(param)
    t = length(xu) == 7 ? xu[end] : zero(eltype(xu))
    v = xu[SA[4:6...]]
    B = B_field(xu, t)
    B² = B[1]^2 + B[2]^2 + B[3]^2

    if use_vE
        E_field = get_EField(param)
        E = E_field(xu, t)
        v = v - (E × B) / B²
    end

    # vector of Larmor radius
    ρ = B × v ./ (q2m * B²)

    return X = @views xu[1:3] - ρ
end

function get_gc(x, y, z, vx, vy, vz, bx, by, bz, q2m)
    l = SVector{3}(x, y, z)
    B = SVector{3}(bx, by, bz)
    v = SVector{3}(vx, vy, vz)

    B² = bx^2 + by^2 + bz^2
    # vector of Larmor radius
    ρ = B × v ./ (q2m * B²)

    return X = l - ρ
end

function get_gc(
        x::T,
        y::T,
        z::T,
        vx::T,
        vy::T,
        vz::T,
        bx::U,
        by::U,
        bz::U,
        q2m
    ) where {T <: AbstractVector, U <: AbstractVector}
    X = Vector{SVector{3, eltype(x)}}(undef, length(x))
    @inbounds for i in eachindex(x, y, z, vx, vy, vz, bx, by, bz)
        X[i] = get_gc(x[i], y[i], z[i], vx[i], vy[i], vz[i], bx[i], by[i], bz[i], q2m)
    end

    return X
end

function get_gc(x, y, z, vx, vy, vz, bx, by, bz, q, m)
    return get_gc(x, y, z, vx, vy, vz, bx, by, bz, q / m)
end

"""
    get_gc_func(param)

Return the function for plotting the orbit of guiding center.

# Example

```julia
param = prepare(E, B; species = Proton)
# The definitions of stateinit, tspan, E and B are ignored.
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern7(); dt = 2e-11)

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1], aspect = :data)
gc = param |> get_gc_func
gc_plot(x, y, z, vx, vy, vz) = (gc(SA[x, y, z, vx, vy, vz])...,)
lines!(ax, sol, idxs = (gc_plot, 1, 2, 3, 4, 5, 6))
```
"""
get_gc_func(param; use_vE::Bool = false) = gc(xu) = get_gc(xu, param; use_vE)

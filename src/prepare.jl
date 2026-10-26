# Construction of tracing parameters.

"""
    is_time_dependent(f::Function)

Judge whether the field function is time dependent.
"""
is_time_dependent(f::Function) = applicable(f, zeros(3), 0.0) || applicable(f, zeros(6), 0.0)

is_time_dependent(::AbstractField{itd}) where {itd} = itd
is_time_dependent(::AbstractFieldInterpolator) = false

"""
    Field{itd, F} <: AbstractField{itd}

A representation of a field function `f`, defined by:

time-independent field

```math
\\mathbf{F} = F(\\mathbf{x}),
```

time-dependent field

```math
\\mathbf{F} = F(\\mathbf{x}, t).
```

# Arguments

  - `field_function::Function`: the function of field.
  - `itd::Bool`: whether the field function is time dependent.
  - `F`: the type of `field_function`.
"""
struct Field{itd, F} <: AbstractField{itd}
    field_function::F
    function Field{itd, F}(field_function::F) where {itd, F}
        return isa(itd, Bool) ? new(field_function) : throw(ArgumentError("itd must be a boolean."))
    end
end

Field(f::Function) = Field{is_time_dependent(f), typeof(f)}(f)

@inline (f::Field{true})(xu, t) = f.field_function(xu, t)
@inline function (f::Field{true})(xu)
    throw(ArgumentError("Time-dependent field function must have a time argument."))
end
@inline (f::Field{false})(xu, t) = f.field_function(xu)
@inline (f::Field{false})(xu) = f.field_function(xu)

function jacobian(f::Field{true}, xu, t)
    return jacobian(f.field_function, xu, t)
end

function jacobian(f::Field{false}, xu, t)
    return jacobian(f.field_function, xu)
end

function jacobian(f::Function, xu, t)
    return ForwardDiff.jacobian(r -> f(r, t), xu)
end

function derivative_t(f::Field{true}, xu, t)
    return derivative_t(f.field_function, xu, t)
end

function derivative_t(f::Field{false}, xu, t)
    return zero(xu)
end

function derivative_t(f::Function, xu, t)
    return ForwardDiff.derivative(τ -> f(xu, τ), t)
end

function Base.show(io::IO, f::Field)
    println(io, "Field with interpolation support")
    return println(io, "Time-dependent: ", is_time_dependent(f.field_function))
end


prepare_field(f, args...; kwargs...) = Field(f)
prepare_field(f::ZeroField, args...; kwargs...) = f
prepare_field(::Nothing, args...; kwargs...) = ZeroField()


function prepare_field(f::AbstractArray, x...; gridtype, order, bc, kw...)
    return Field(build_interpolator(gridtype, f, x..., order, bc; kw...))
end

function _prepare(
        E::TE, B::TB, F::TF, args...; species = Proton, q = nothing,
        m = nothing, c = nothing, gridtype::Type{G} = CartesianGrid, type = nothing, kw...
    ) where {TE, TB, TF, G}
    q = @something q species.q
    m = @something m species.m
    if !isnothing(type)
        q, m = type(q), type(m)
    end
    E, B, F = map(f -> prepare_field(f, args...; gridtype, kw...), (E, B, F))
    p = (; q2m = q / m, m, E, B, F)
    isnothing(c) && return p
    return (; p..., c = isnothing(type) ? c : type(c))
end

"""
    prepare([grid...;] E = nothing, B = nothing, F = nothing, kwargs...) -> (; q2m, m, E, B, F[, c])
    prepare([grid...,] E, B, F = nothing; kwargs...)
    prepare(B; E = nothing, F = nothing, kwargs...)

Tracing parameters for a particle in the electric field `E`, the magnetic field `B` and
the external force `F`. Omitted fields are zero. Prefer the keyword form, which names each
field; the positional form takes **`E` before `B`**, and a single positional field is `B`.

```julia
prepare(; B)              # magnetic field only
prepare(; E, B)           # analytic fields
prepare(x, y, z; E, B)    # arrays interpolated on a 3D grid
prepare(grid; B, F)       # Meshes.jl grid, with an external force
```

Each field is either a function `f(x)` or `f(x, t)` returning a 3-vector, or an array of
size `(3, size(grid)...)`, interpolated on the grid. The grid is given as
  - `x` (1D, along the axis `dir = 1, 2, 3`), `x, y` (2D) or `x, y, z` (3D): sorted vectors,
    or `r, θ, ϕ` for a spherical grid (`gridtype = StructuredGrid`, implied by a `LogRange`
    `r`), with field components `(Br, Bθ, Bϕ)`;
  - a `CartesianGrid` or `RectilinearGrid` from Meshes.jl.

# Keywords

  - `species = Proton`: `Proton`, `Electron` or `Ion(m̄, q̄)` with mass and charge numbers.
  - `q`, `m`: particle charge and mass, overriding `species`.
  - `c`: speed of light for the relativistic equations, e.g. `c = 1` in normalized units;
    `c` in SI when omitted.
  - `type`: number type to convert `q2m`, `m` (and `c`) to, e.g. `Float32`.
  - `order = 1`: interpolation order, 0, 1 or 3.
  - `bc`: boundary condition from `FastInterpolations.jl`; `FillExtrap(NaN)` by default,
    `ClampExtrap()` in 1D and `WrapExtrap()` on spherical grids.
  - `store = StorePolicy()`: interpolation storage policy; `StorePolicy(copy = false)` skips
    copying field and grid arrays that stay unchanged for the interpolator's lifetime.
  - `gridtype`: `CartesianGrid`, `RectilinearGrid` or `StructuredGrid`.
"""
function prepare(
        x::AbstractVector, y::AbstractVector, E, B, F = nothing; order = 1,
        bc = FillExtrap(NaN), kw...
    )
    @assert issorted(x) "Grid vector `x` must be sorted."
    @assert issorted(y) "Grid vector `y` must be sorted."
    return _prepare(E, B, F, x, y; gridtype = CartesianGrid, order, bc, kw...)
end

function prepare(
        x::AbstractVector, y::AbstractVector, z::AbstractVector, E, B, F = nothing; order = 1,
        bc = FillExtrap(NaN), gridtype = CartesianGrid, kw...
    )
    @assert issorted(x) "Grid vector `x` must be sorted."
    @assert issorted(y) "Grid vector `y` must be sorted."
    @assert issorted(z) "Grid vector `z` must be sorted."
    return _prepare(E, B, F, x, y, z; gridtype, order, bc, kw...)
end

function prepare(
        x::Base.LogRange, y::AbstractVector, z::AbstractVector, E, B, F = nothing;
        order = 1, bc = WrapExtrap(), kw...
    )
    @assert issorted(x) "Grid vector `x` must be sorted."
    @assert issorted(y) "Grid vector `y` must be sorted."
    @assert issorted(z) "Grid vector `z` must be sorted."
    return _prepare(E, B, F, x, y, z; gridtype = StructuredGrid, order, bc, kw...)
end

function prepare(x::AbstractVector, E, B, F = nothing; order = 1, bc = ClampExtrap(), dir = 1, kw...)
    @assert issorted(x) "Grid vector `x` must be sorted."
    return _prepare(E, B, F, x; gridtype = CartesianGrid, order, bc, dir, kw...)
end
prepare(E, B, F = nothing; kw...) = _prepare(E, B, F; kw...)
prepare(B; E = nothing, F = nothing, kw...) = _prepare(E, B, F; kw...)

# Keyword forms: the field roles are named, so they cannot be swapped.
prepare(; E = nothing, B = nothing, F = nothing, kw...) = _prepare(E, B, F; kw...)
prepare(x::AbstractVector; E = nothing, B = nothing, F = nothing, kw...) =
    prepare(x, E, B, F; kw...)
prepare(
    x::AbstractVector, y::AbstractVector;
    E = nothing, B = nothing, F = nothing, kw...
) = prepare(x, y, E, B, F; kw...)
prepare(
    x::AbstractVector, y::AbstractVector, z::AbstractVector;
    E = nothing, B = nothing, F = nothing, kw...
) = prepare(x, y, z, E, B, F; kw...)

"""
    PotentialField(; phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing)

Scalar potential `phi(x[, t])` and vector potential `A(x[, t])` for canonical tracing;
`nothing` is an absent potential. Gradients that are not given are computed with
`ForwardDiff`. `grad_A(x[, t])[i, j]` is ``∂A_j/∂x_i``, so that `grad_A * v` is
``∇(\\mathbf{A} ⋅ \\mathbf{v})`` at fixed ``\\mathbf{v}``.
"""
function PotentialField(; phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing)
    f_phi, f_A = prepare_field(phi), prepare_field(A)
    return PotentialField(f_phi, f_A, _grad_phi(grad_phi, f_phi), _grad_A(grad_A, f_A))
end

# Positional forms, kept for compatibility until the next breaking release.
PotentialField(phi, A; kw...) = PotentialField(; phi, A, kw...)
PotentialField(A; phi = nothing, kw...) = PotentialField(; phi, A, kw...)

_grad_phi(g, phi) = prepare_field(g)
_grad_phi(::Nothing, ::ZeroField) = ZeroField()
_grad_phi(::Nothing, phi) = (x, t) -> SVector{3}(ForwardDiff.gradient(r -> phi(r, t), x))

_grad_A(g, A) = prepare_field(g)
_grad_A(::Nothing, ::ZeroField) = ZeroField()
_grad_A(::Nothing, A) = (x, t) -> SMatrix{3, 3}(ForwardDiff.jacobian(r -> A(r, t), x)')

"""
    prepare(pf::PotentialField; species = Proton, q = nothing, m = nothing, c = c, type = nothing)

Return parameters `(; q, m, c, pf)` for `TraceCanonicalProblem`.
"""
function prepare(
        pf::PotentialField;
        species = Proton, q = nothing, m = nothing, c = c, type = nothing, kw...
    )
    q = @something q species.q
    m = @something m species.m
    p = (; q, m, c, pf)
    return isnothing(type) ? p : _promote_p(type, p)
end

"""
    prepare_canonical(; phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing, kwargs...)

[`PotentialField`](@ref) and [`prepare`](@ref) in one call, for `TraceCanonicalProblem`.
"""
prepare_canonical(; phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing, kw...) =
    prepare(PotentialField(; phi, A, grad_phi, grad_A); kw...)

# Positional forms, kept for compatibility until the next breaking release.
prepare_canonical(phi, A; kw...) = prepare_canonical(; phi, A, kw...)
prepare_canonical(A; phi = nothing, kw...) = prepare_canonical(; phi, A, kw...)

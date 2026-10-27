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

@inline (f::Field{true})(xu, t) = getfield(f, :field_function)(xu, t)
@inline function (f::Field{true})(xu)
    throw(ArgumentError("Time-dependent field function must have a time argument."))
end
@inline (f::Field{false})(xu, t) = getfield(f, :field_function)(xu)
@inline (f::Field{false})(xu) = getfield(f, :field_function)(xu)

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

function Base.getproperty(f::Field, sym::Symbol)
    if sym === :field_function
        fn = getfield(f, :field_function)
        return fn isa ErasedFunction ? fn.raw : fn
    end
    return getfield(f, sym)
end

include("erased.jl")

prepare_field(f, args...; erased = false, kwargs...) = erased ? ErasedField(f) : Field(f)
prepare_field(
    f::Field{itd, ErasedFunction}, args...; erased = false, kwargs...
) where {itd} = f
prepare_field(f::Field, args...; erased = false, kwargs...) = erased ? ErasedField(f) : f
prepare_field(f::ZeroField, args...; kwargs...) = f
prepare_field(::Nothing, args...; kwargs...) = ZeroField()

function prepare_field(f::AbstractArray, x...; gridtype, order, bc, erased = false, kw...)
    itp = build_interpolator(gridtype, f, x..., order, bc; kw...)
    return erased ? ErasedField(itp) : Field(itp)
end

function _prepare(
        E::TE, B::TB, F::TF, args...; species = Proton, q = nothing,
        m = nothing, c = nothing, gridtype::Type{G} = CartesianGrid, type = nothing,
        erased::Bool = false, kw...
    ) where {TE, TB, TF, G}
    q = @something q species.q
    m = @something m species.m
    if !isnothing(type)
        q, m = type(q), type(m)
    end
    E, B, F = map(f -> prepare_field(f, args...; gridtype, erased, kw...), (E, B, F))
    p = (; q2m = q / m, m, E, B, F)
    isnothing(c) && return p
    return (; p..., c = isnothing(type) ? c : type(c))
end

"""
    prepare([grid...;] E = nothing, B = nothing, F = nothing, kwargs...) -> (; q2m, m, E, B, F[, c])

Tracing parameters for a particle in the electric field `E`, the magnetic field `B` and
the external force `F`. Omitted fields are zero.

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
    `ClampExtrap()` in 1D. Ignored on spherical grids, which fill `r` and `θ` with `NaN`
    and wrap `ϕ`.
  - `store = StorePolicy()`: interpolation storage policy; `StorePolicy(copy = false)` skips
    copying field and grid arrays that stay unchanged for the interpolator's lifetime.
  - `gridtype`: `CartesianGrid`, `RectilinearGrid` or `StructuredGrid`.
"""
function prepare(
        grid::AbstractVector...; E = nothing, B = nothing, F = nothing,
        gridtype = _default_gridtype(grid), order = 1, bc = _default_bc(grid),
        erased::Bool = false, kw...
    )
    all(issorted, grid) || throw(ArgumentError("Grid vectors must be sorted."))
    return _prepare(E, B, F, grid...; gridtype, order, bc, erased, kw...)
end

"""
    prepare_erased([grid...;] E = nothing, B = nothing, F = nothing, kwargs...)

Tracing parameters with type-erased fields for uniform dispatch and static compilation.
Equivalent to `prepare(args...; erased = true, kwargs...)`.
"""
prepare_erased(args...; kw...) = prepare(args...; erased = true, kw...)

_default_gridtype(::Tuple{Base.LogRange, Any, Any}) = StructuredGrid
_default_gridtype(_) = CartesianGrid
_default_bc(::Tuple{Any}) = ClampExtrap()
_default_bc(_) = FillExtrap(NaN)

"""
    PotentialField(;
        phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing, erased = false
    )

Scalar potential `phi(x[, t])` and vector potential `A(x[, t])` for canonical tracing;
`nothing` is an absent potential. Gradients that are not given are computed with
`ForwardDiff`. `grad_A(x[, t])[i, j]` is ``∂A_j/∂x_i``, so that `grad_A * v` is
``∇(\\mathbf{A} ⋅ \\mathbf{v})`` at fixed ``\\mathbf{v}``.
"""
_prepare_potential_scalar(f, erased) = erased ? ErasedField(f) : Field(f)
_prepare_potential_scalar(f::ZeroField, erased) = f
_prepare_potential_scalar(::Nothing, erased) = ZeroField()

function PotentialField(;
        phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing,
        erased::Bool = false, kw...
    )
    f_phi, f_A = _prepare_potential_scalar(phi, erased), prepare_field(A; erased, kw...)
    return PotentialField(
        f_phi, f_A, _grad_phi(grad_phi, f_phi; erased, kw...),
        _grad_A(grad_A, f_A; erased, kw...)
    )
end

_grad_phi(g, phi; erased = false, kw...) = prepare_field(g; erased, kw...)
_grad_phi(::Nothing, ::ZeroField; kw...) = ZeroField()
_grad_phi(::Nothing, phi; kw...) =
    (x, t) -> SVector{3}(ForwardDiff.gradient(r -> phi(r, t), x))

_grad_A(g::ZeroField, A; kw...) = g
_grad_A(g, A; erased = false, kw...) = erased ? ErasedField(g) : Field(g)
_grad_A(::Nothing, ::ZeroField; kw...) = ZeroField()
_grad_A(::Nothing, A; kw...) =
    (x, t) -> SMatrix{3, 3}(ForwardDiff.jacobian(r -> A(r, t), x)')

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
prepare_canonical(;
    phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing,
    erased::Bool = false, kw...
) =
    prepare(PotentialField(; phi, A, grad_phi, grad_A, erased); kw...)

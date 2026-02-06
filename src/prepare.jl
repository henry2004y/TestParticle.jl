# Construction of tracing parameters.

"""
Judge whether the field function is time dependent.
"""
is_time_dependent(f::Function) = applicable(f, zeros(6), 0.0) ? true : false

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

is_time_dependent(::AbstractField{itd}) where {itd} = itd
# Note: Without it, AbstractFieldInterpolators are treated as time-dependent (due to their (x,t) method).
# This causes TestParticle to wrap them in Field{true}, which forbids calling f(x) (without time).
# Several components/tests rely on f(x) for static fields, leading to ArgumentError.
#TODO: Have a proper treatment for the time dependency with LazyTimeInterpolator.
is_time_dependent(::AbstractFieldInterpolator) = false # Always treat as static by default

@inline (f::Field{true})(xu, t) = f.field_function(xu, t)
@inline function (f::Field{true})(xu)
    throw(ArgumentError("Time-dependent field function must have a time argument."))
end
@inline (f::Field{false})(xu, t) = SVector{3}(f.field_function(xu))
@inline (f::Field{false})(xu) = SVector{3}(f.field_function(xu))

function Base.show(io::IO, f::Field)
    println(io, "Field with interpolation support")
    return println(io, "Time-dependent: ", is_time_dependent(f.field_function))
end


get_q2m(param) = param[1]
get_BField(param) = param[4]
get_EField(param) = param[3]

prepare_field(f, args...; kwargs...) = Field(f)
prepare_field(f::ZeroField, args...; kwargs...) = f

function prepare_field(f::AbstractArray, x...; gridtype, order, bc, kw...)
    return Field(getinterp(gridtype, f, x..., order, bc; kw...))
end

function _prepare(
        E, B, F, args...; species = Proton, q = nothing,
        m = nothing, gridtype = CartesianGrid, kw...
    )
    q = @something q species.q
    m = @something m species.m
    q2m = q / m
    fE = prepare_field(E, args...; gridtype, kw...)
    fB = prepare_field(B, args...; gridtype, kw...)
    fF = prepare_field(F, args...; gridtype, kw...)
    return q2m, m, fE, fB, fF
end

"""
    prepare(args...; kwargs...) -> (q2m, m, E, B, F)
    prepare(E, B, F = ZeroField(); kwargs...)
    prepare(grid::CartesianGrid, E, B, F = ZeroField(); kwargs...)
    prepare(x, E, B, F = ZeroField(); dir = 1, kwargs...)
    prepare(x, y, E, B, F = ZeroField(); kwargs...)
    prepare(x, y, z, E, B, F = ZeroField(); kwargs...)
    prepare(B; E = ZeroField(), F = ZeroField(), kwargs...)

Return a tuple consists of particle charge-mass ratio for a prescribed `species` of charge `q` and mass `m`,
mass `m` for a prescribed `species`, analytic/interpolated EM field functions, and external force `F`.

Prescribed `species` are `Electron` and `Proton`;
other species can be manually specified with `m` and `q` keywords or `species = Ion(m̄, q̄)`,
where `m̄` and `q̄` are the mass and charge numbers respectively.

Direct range input for uniform grid in 1/2/3D is supported.
For 1D grid, an additional keyword `dir` is used for specifying the spatial direction, 1 -> x, 2 -> y, 3 -> z.
For 3D grid, the default grid type is `CartesianGrid`. To use `StructuredGrid` (spherical) grid, an additional keyword `gridtype` is needed.
For `StructuredGrid` (spherical) grid, dimensions of field arrays should be `(Br, Bθ, Bϕ)`.

# Keywords

  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic.
  - `species=Proton`: particle species.
  - `q=nothing`: particle charge.
  - `m=nothing`: particle mass.
  - `gridtype`: `CartesianGrid`, `RectilinearGrid`, `StructuredGrid`.
"""
function prepare(grid::CartesianGrid, E, B, F = ZeroField(); order = 1, bc = 1, kw...)
    return _prepare(E, B, F, makegrid(grid)...; gridtype = CartesianGrid, order, bc, kw...)
end

function prepare(grid::RectilinearGrid, E, B, F = ZeroField(); order = 1, bc = 1, kw...)
    return _prepare(E, B, F, makegrid(grid)...; gridtype = RectilinearGrid, order, bc, kw...)
end

function prepare(
        x::T, y::T, E, B, F = ZeroField(); order = 1,
        bc = 1, kw...
    ) where {T <: AbstractRange}
    return _prepare(E, B, F, x, y; gridtype = CartesianGrid, order, bc, kw...)
end

function prepare(
        x::T, y::T, z::T, E, B, F = ZeroField(); order = 1,
        bc = 1, gridtype = CartesianGrid, kw...
    ) where {T <: AbstractRange}
    return _prepare(E, B, F, x, y, z; gridtype, order, bc, kw...)
end

function prepare(
        x::Base.LogRange, y::T, z::T, E, B, F = ZeroField();
        order = 1, bc = 2, kw...
    ) where {T <: AbstractRange}
    return _prepare(E, B, F, x, y, z; gridtype = StructuredGrid, order, bc, kw...)
end

function prepare(x::AbstractRange, E, B, F = ZeroField(); order = 1, bc = 3, dir = 1, kw...)
    return _prepare(E, B, F, x; gridtype = CartesianGrid, order, bc, dir, kw...)
end
prepare(E, B, F = ZeroField(); kw...) = _prepare(E, B, F; kw...)
prepare(B; E = ZeroField(), F = ZeroField(), kw...) = _prepare(E, B, F; kw...)

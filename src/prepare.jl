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

adapt_field_to_gpu(field::Field, ::CPU) = field
adapt_field_to_gpu(field::ZeroField, ::CPU) = field
adapt_field_to_gpu(field::ZeroField, ::Backend) = field

function adapt_field_to_gpu(field::Field, backend::Backend)
    backend isa CPU && return field

    adapted_func = adapt_field_to_gpu(field.field_function, backend)
    return Field{is_time_dependent(field), typeof(adapted_func)}(adapted_func)
end

adapt_field_to_gpu(f::Function, backend::Backend) = Adapt.adapt(backend, f)

function adapt_field_to_gpu(fi::FieldInterpolator, backend::Backend)
    backend isa CPU && return fi
    return _to_gpu_grid(fi.itp, backend)
end

function adapt_field_to_gpu(fi::FieldInterpolator2D, backend::Backend)
    backend isa CPU && return fi
    return _to_gpu_grid(fi.itp, backend)
end

function adapt_field_to_gpu(fi::FieldInterpolator1D, backend::Backend)
    backend isa CPU && return fi
    return _to_gpu_grid(fi.itp, backend; dir = fi.dir)
end

function adapt_field_to_gpu(fi::SphericalFieldInterpolator, backend::Backend)
    backend isa CPU && return fi
    return _to_gpu_spherical_grid(fi.itp, backend)
end

adapt_field_to_gpu(g::GPUGrid3D, backend::Backend) =
    backend isa CPU ? g : Adapt.adapt(backend, g)
adapt_field_to_gpu(g::GPUGrid2D, backend::Backend) =
    backend isa CPU ? g : Adapt.adapt(backend, g)
adapt_field_to_gpu(g::GPUGrid1D, backend::Backend) =
    backend isa CPU ? g : Adapt.adapt(backend, g)
adapt_field_to_gpu(g::GPUSphericalGrid, backend::Backend) =
    backend isa CPU ? g : Adapt.adapt(backend, g)

function prepare_field(f::AbstractArray, x...; gridtype, order, bc, kw...)
    return Field(build_interpolator(gridtype, f, x..., order, bc; kw...))
end

function _prepare(
        E, B, F, args...; species = Proton, q = nothing,
        m = nothing, gridtype = CartesianGrid, type = nothing, kw...
    )
    if type !== nothing
        T = type
        sp = species isa Species ? Species{T}(species) : species
        q_val = isnothing(q) ? sp.q : T(q)
        m_val = isnothing(m) ? sp.m : T(m)
        q2m = q_val / m_val
        fE = prepare_field(E, args...; gridtype, kw...)
        fB = prepare_field(B, args...; gridtype, kw...)
        fF = prepare_field(F, args...; gridtype, kw...)
        return q2m, m_val, fE, fB, fF
    else
        q = @something q species.q
        m = @something m species.m
        q2m = q / m
        fE = prepare_field(E, args...; gridtype, kw...)
        fB = prepare_field(B, args...; gridtype, kw...)
        fF = prepare_field(F, args...; gridtype, kw...)
        return q2m, m, fE, fB, fF
    end
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

Direct range input for uniform grid in 1/2/3D is supported. The grid vectors must be sorted.
For 1D grid, an additional keyword `dir` is used for specifying the spatial direction, 1 -> x, 2 -> y, 3 -> z.
For 3D grid, the default grid type is `CartesianGrid`. To use `StructuredGrid` (spherical) grid, an additional keyword `gridtype` is needed.
For `StructuredGrid` (spherical) grid, dimensions of field arrays should be `(Br, Bθ, Bϕ)`.

# Keywords

  - `order::Int=1`: order of interpolation in [0,1,3].
  - `bc=FillExtrap(NaN)`: boundary condition type from `FastInterpolations.jl`.
  - `store=StorePolicy()`: interpolation storage policy. Use
    `StorePolicy(copy=false)` for zero-copy construction when field and grid inputs
    will remain unchanged for the interpolator's lifetime.
  - `species=Proton`: particle species.
  - `q=nothing`: particle charge.
  - `m=nothing`: particle mass.
  - `gridtype`: `CartesianGrid`, `RectilinearGrid`, `StructuredGrid`.
"""
function prepare(
        x::AbstractVector, y::AbstractVector, E, B, F = ZeroField(); order = 1,
        bc = FillExtrap(NaN), kw...
    )
    @assert issorted(x) "Grid vector `x` must be sorted."
    @assert issorted(y) "Grid vector `y` must be sorted."
    return _prepare(E, B, F, x, y; gridtype = CartesianGrid, order, bc, kw...)
end

function prepare(
        x::AbstractVector, y::AbstractVector, z::AbstractVector, E, B, F = ZeroField(); order = 1,
        bc = FillExtrap(NaN), gridtype = CartesianGrid, kw...
    )
    @assert issorted(x) "Grid vector `x` must be sorted."
    @assert issorted(y) "Grid vector `y` must be sorted."
    @assert issorted(z) "Grid vector `z` must be sorted."
    return _prepare(E, B, F, x, y, z; gridtype, order, bc, kw...)
end

function prepare(
        x::Base.LogRange, y::AbstractVector, z::AbstractVector, E, B, F = ZeroField();
        order = 1, bc = WrapExtrap(), kw...
    )
    @assert issorted(x) "Grid vector `x` must be sorted."
    @assert issorted(y) "Grid vector `y` must be sorted."
    @assert issorted(z) "Grid vector `z` must be sorted."
    return _prepare(E, B, F, x, y, z; gridtype = StructuredGrid, order, bc, kw...)
end

function prepare(x::AbstractVector, E, B, F = ZeroField(); order = 1, bc = ClampExtrap(), dir = 1, kw...)
    @assert issorted(x) "Grid vector `x` must be sorted."
    return _prepare(E, B, F, x; gridtype = CartesianGrid, order, bc, dir, kw...)
end
prepare(E, B, F = ZeroField(); kw...) = _prepare(E, B, F; kw...)
prepare(B; E = ZeroField(), F = ZeroField(), kw...) = _prepare(E, B, F; kw...)

"""
    PotentialField(phi, A; grad_phi = nothing, grad_A = nothing)
    PotentialField(A; phi = ZeroField(), grad_phi = nothing, grad_A = nothing)

Construct a `PotentialField` for canonical tracing. When `grad_phi` or `grad_A` are omitted,
they are computed automatically via `ForwardDiff`.
"""
function PotentialField(
        phi, A;
        grad_phi = nothing, grad_A = nothing
    )
    f_phi = phi isa ZeroField ? phi : Field(phi)
    f_A = A isa ZeroField ? A : Field(A)

    g_phi = if grad_phi !== nothing
        applicable(grad_phi, zeros(3), 0.0) ? grad_phi : ((x, t) -> grad_phi(x))
    elseif f_phi isa ZeroField
        ZeroField()
    else
        (x, t) -> SVector{3}(ForwardDiff.gradient(r -> f_phi(r, t), x))
    end

    g_A = if grad_A !== nothing
        applicable(grad_A, zeros(3), 0.0) ? grad_A : ((x, t) -> grad_A(x))
    elseif f_A isa ZeroField
        (x, t) -> zero(SMatrix{3, 3, eltype(x), 9})
    else
        (x, t) -> SMatrix{3, 3}(ForwardDiff.jacobian(r -> f_A(r, t), x)')
    end

    return PotentialField(f_phi, f_A, g_phi, g_A)
end

PotentialField(A; phi = ZeroField(), grad_phi = nothing, grad_A = nothing) =
    PotentialField(phi, A; grad_phi, grad_A)

"""
    prepare(pf::PotentialField; species = Proton, q = nothing, m = nothing, c = c, type = nothing)

Return parameters `(q, m, c, pf)` for `TraceCanonicalProblem`.
"""
function prepare(
        pf::PotentialField;
        species = Proton, q = nothing, m = nothing, c = c, type = nothing, kw...
    )
    if type !== nothing
        T = type
        sp = species isa Species ? Species{T}(species) : species
        q_val = isnothing(q) ? sp.q : T(q)
        m_val = isnothing(m) ? sp.m : T(m)
        c_val = T(c)
        return q_val, m_val, c_val, pf
    else
        q_val = @something q species.q
        m_val = @something m species.m
        c_val = c
        return q_val, m_val, c_val, pf
    end
end

"""
    prepare_canonical(phi, A; kwargs...)
    prepare_canonical(A; phi = ZeroField(), kwargs...)

Convenience wrapper to construct a `PotentialField` and prepare parameters for `TraceCanonicalProblem`.
"""
prepare_canonical(phi, A; grad_phi = nothing, grad_A = nothing, kw...) =
    prepare(PotentialField(phi, A; grad_phi, grad_A); kw...)

prepare_canonical(A; phi = ZeroField(), grad_phi = nothing, grad_A = nothing, kw...) =
    prepare(PotentialField(phi, A; grad_phi, grad_A); kw...)

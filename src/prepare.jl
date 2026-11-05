# Construction of tracing parameters.

"""
    is_time_dependent(f::Function)

Judge whether the field function is time dependent.
"""
is_time_dependent(f) = applicable(f, zeros(3), 0.0) || applicable(f, zeros(6), 0.0)

is_time_dependent(::AbstractField{itd}) where {itd} = itd
is_time_dependent(::AbstractFieldInterpolator) = false

"""
    Field(f; time_dependent=is_time_dependent(f))

Wrap `f(x)` or `f(x, t)` as a field. Automatic detection prefers the two-argument
signature; set `time_dependent` explicitly for ambiguous or input-type-specific callables.
"""
struct Field{itd, F} <: AbstractField{itd}
    field_function::F
    function Field{itd, F}(field_function::F) where {itd, F}
        return isa(itd, Bool) ? new(field_function) : throw(ArgumentError("itd must be a boolean."))
    end
end

Field(f; time_dependent::Bool = is_time_dependent(f)) =
    Field{time_dependent, typeof(f)}(f)

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

function jacobian(f, xu, t)
    return ForwardDiff.jacobian(r -> f(r, t), xu)
end

function derivative_t(f::Field{true}, xu, t)
    return derivative_t(f.field_function, xu, t)
end

function derivative_t(f::Field{false}, xu, t)
    return zero(xu)
end

function derivative_t(f, xu, t)
    return ForwardDiff.derivative(τ -> f(xu, τ), t)
end

function Base.show(io::IO, f::Field)
    println(io, "Field with interpolation support")
    return println(io, "Time-dependent: ", is_time_dependent(f))
end

function Base.getproperty(f::Field, sym::Symbol)
    if sym === :field_function
        fn = getfield(f, :field_function)
        return fn isa ErasedFunction ? fn.raw : fn
    end
    return getfield(f, sym)
end

include("erased.jl")

prepare_field(f; erased = false, kwargs...) = erased ? ErasedField(f) : Field(f)
prepare_field(
    f::Field{itd, ErasedFunction}; erased = false, kwargs...
) where {itd} = f
prepare_field(f::Field; erased = false, kwargs...) = erased ? ErasedField(f) : f
prepare_field(f::ZeroField; kwargs...) = f
prepare_field(::Nothing; kwargs...) = ZeroField()


"""
    prepare(; E=nothing, B=nothing, F=nothing, species=Proton, q, m, c, type, erased=false)

Build particle parameters from callable electric, magnetic, and external force fields.
Omitted fields are zero. Each callable accepts `x` or `(x, t)` and returns three components.
Construct numerical fields with [`build_interpolator`](@ref) before calling `prepare`.
See [`Field`](@ref) to select the calling convention. `q` and `m` override `species`;
`c` sets the relativistic speed of light. `type` converts particle constants.
"""
function prepare(;
        E = nothing, B = nothing, F = nothing, species = Proton,
        q = nothing, m = nothing, c = nothing, type = nothing, erased::Bool = false
    )
    any(f -> f isa AbstractArray, (E, B, F)) && throw(
        ArgumentError(
            "Array fields require build_interpolator(A, grids; bc)."
        )
    )
    q = @something q species.q
    m = @something m species.m
    if !isnothing(type)
        q, m = type(q), type(m)
    end
    E, B, F = map(f -> prepare_field(f; erased), (E, B, F))
    p = (; q2m = q / m, m, E, B, F)
    isnothing(c) && return p
    return (; p..., c = isnothing(type) ? c : type(c))
end

"""
    prepare_erased(; E = nothing, B = nothing, F = nothing, kwargs...)

Tracing parameters with type-erased fields for uniform dispatch and static compilation.
Equivalent to `prepare(; erased = true, kwargs...)`.
"""
prepare_erased(; kw...) = prepare(; erased = true, kw...)

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

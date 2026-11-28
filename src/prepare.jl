# Construction of tracing parameters.

"""
    is_time_dependent(f::Function)

Judge whether the field function is time dependent.
"""
is_time_dependent(f) = applicable(f, zeros(3), 0.0) || applicable(f, zeros(6), 0.0)

is_time_dependent(::AbstractField{itd}) where {itd} = itd
is_time_dependent(::AbstractFieldInterpolator) = false
is_time_dependent(::Union{LazyTimeInterpolator, _LazyFieldComponent}) = true

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

prepare_field(f; erased = false) = erased ? ErasedField(f) : Field(f)
prepare_field(f::Field{itd, ErasedFunction}; erased = false) where {itd} = f
prepare_field(f::Field; erased = false) = erased ? ErasedField(f) : f
prepare_field(f::ZeroField; erased = false) = f
prepare_field(::Nothing; erased = false) = ZeroField()


"""
    prepare(; E=nothing, B=nothing, F=nothing, species=Proton, q, m, c, erased=false)

Build particle parameters `(; q, q2m, m, E, B, F)` from callable electric, magnetic, and
external force fields, for every model. Omitted fields are zero. Each callable accepts `x` or
`(x, t)` and returns three components. Construct numerical fields with
[`build_interpolator`](@ref) before calling `prepare`. See [`Field`](@ref) to select the
calling convention. `q` and `m` override `species`; `c` sets the relativistic speed of light.
Particle constants take the float type of the state when traced. `erased = true` wraps each
field in an [`ErasedField`](@ref), whose type does not depend on the field function, for
uniform dispatch and static compilation.

# Examples
```jldoctest
julia> B(x) = SA[0.0, 0.0, 1.0e-8];

julia> p = prepare(; B, species = Proton);

julia> P = TestParticle.get_gyroperiod(1.0e-8);

julia> u0 = SA[1.0, 0.0, 0.0, 0.0, 1.0, 0.0];

julia> prob = TraceProblem(u0, (0.0, P), p);

julia> sol1 = solve(prob, Vern7(); reltol = 1.0e-9, abstol = 1.0e-9);

julia> sol2 = solve(prob, Boris(); dt = P / 1000);

julia> isapprox(sol1.u[end], u0; atol = 1.0e-6)
true

julia> isapprox(sol2.u[end], u0; atol = 1.0e-4)
true
```
"""
function prepare(;
        E = nothing, B = nothing, F = nothing, species = Proton,
        q = nothing, m = nothing, c = nothing, erased::Bool = false
    )
    any(f -> f isa AbstractArray, (E, B, F)) && throw(
        ArgumentError(
            "Array fields require build_interpolator(A, grids; bc)."
        )
    )
    q = @something q species.q
    m = @something m species.m
    E, B, F = map(f -> prepare_field(f; erased), (E, B, F))
    p = (; q, q2m = q / m, m, E, B, F)
    isnothing(c) && return p
    return (; p..., c)
end

"""
    PotentialField(;
        phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing, erased = false
    )

Scalar potential `phi(x[, t])` and vector potential `A(x[, t])` for canonical tracing;
`nothing` is an absent potential. Gradients that are not given are computed with
`ForwardDiff`. `grad_A(x[, t])[i, j]` is ``∂A_j/∂x_i``, so that `grad_A * v` is
``∇(\\mathbf{A} ⋅ \\mathbf{v})`` at fixed ``\\mathbf{v}``. `erased` applies to `A` and
a given `grad_phi` only.
"""
function PotentialField(;
        phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing,
        erased::Bool = false
    )
    # `ErasedFunction` holds 3-vector fields; scalar `phi` and matrix `grad_A` stay plain.
    f_phi, f_A = prepare_field(phi), prepare_field(A; erased)
    return PotentialField(
        f_phi, f_A, _grad_phi(grad_phi, f_phi; erased), _grad_A(grad_A, f_A)
    )
end

_grad_phi(g, phi; erased) = prepare_field(g; erased)
_grad_phi(::Nothing, ::ZeroField; erased) = ZeroField()
_grad_phi(::Nothing, phi; erased) =
    (x, t) -> SVector{3}(ForwardDiff.gradient(r -> phi(r, t), x))

_grad_A(g, A) = prepare_field(g)
_grad_A(::Nothing, ::ZeroField) = ZeroField()
_grad_A(::Nothing, A) =
    (x, t) -> SMatrix{3, 3}(ForwardDiff.jacobian(r -> A(r, t), x)')

"""
    prepare(pf::PotentialField; species = Proton, q = nothing, m = nothing, c = c)

Return parameters `(; q, m, c, pf)` for `TraceCanonicalProblem`.
"""
function prepare(
        pf::PotentialField; species = Proton, q = nothing, m = nothing, c = c
    )
    q = @something q species.q
    m = @something m species.m
    return (; q, m, c, pf)
end

"""
    prepare_canonical(; phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing, kwargs...)

[`PotentialField`](@ref) and [`prepare`](@ref) in one call, for `TraceCanonicalProblem`.

# Examples
```jldoctest
julia> A(x) = SA[0.0, 1.0e-8 * x[1], 0.0];

julia> p = prepare_canonical(; A, species = Proton);

julia> x0 = SA[0.0, 0.0, 0.0]; v0 = SA[0.0, 1.0, 0.0];

julia> u0 = TestParticle.velocity_to_canonical(x0, v0, p);

julia> prob = TraceCanonicalProblem(u0, (0.0, 1.0), p);

julia> sol = solve(prob, Vern7(); reltol = 1.0e-9, abstol = 1.0e-9);

julia> sol.retcode == TestParticle.ReturnCode.Success
true

julia> isapprox(sol.u[end][5], u0[5]; rtol = 1.0e-9)
true
```
"""
prepare_canonical(;
    phi = nothing, A = nothing, grad_phi = nothing, grad_A = nothing,
    erased::Bool = false, kw...
) =
    prepare(PotentialField(; phi, A, grad_phi, grad_A, erased); kw...)

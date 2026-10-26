using FunctionWrappers: FunctionWrapper

"""
    ErasedFunction

Type-erased callable for 3D vector fields wrapping arbitrary user functions or
interpolators. Holds precompiled `FunctionWrapper` handles for both `Float64` and
`Float32` inputs, along with the raw callable for autodiff (`ForwardDiff.Dual`)
and custom types.
"""
struct ErasedFunction
    f64::FunctionWrapper{SVector{3, Float64}, Tuple{SVector{3, Float64}, Float64}}
    f32::FunctionWrapper{SVector{3, Float32}, Tuple{SVector{3, Float32}, Float32}}
    raw::Any
    itd::Bool
end

@inline (ef::ErasedFunction)(xu::SVector{3, Float64}, t::Float64) = ef.f64(xu, t)
@inline (ef::ErasedFunction)(xu::SVector{3, Float32}, t::Float32) = ef.f32(xu, t)
@inline (ef::ErasedFunction)(xu::SVector{3, Float64}) = ef.f64(xu, 0.0)
@inline (ef::ErasedFunction)(xu::SVector{3, Float32}) = ef.f32(xu, 0.0f0)

@inline (ef::ErasedFunction)(xu::AbstractVector{Float64}, t::Float64) =
    length(xu) == 3 ? ef.f64(SVector{3, Float64}(xu), t) :
    (ef.itd ? ef.raw(xu, t) : ef.raw(xu))
@inline (ef::ErasedFunction)(xu::AbstractVector{Float32}, t::Float32) =
    length(xu) == 3 ? ef.f32(SVector{3, Float32}(xu), t) :
    (ef.itd ? ef.raw(xu, t) : ef.raw(xu))
@inline (ef::ErasedFunction)(xu::AbstractVector{Float64}) =
    length(xu) == 3 ? ef.f64(SVector{3, Float64}(xu), 0.0) :
    (ef.itd ? ef.raw(xu, 0.0) : ef.raw(xu))
@inline (ef::ErasedFunction)(xu::AbstractVector{Float32}) =
    length(xu) == 3 ? ef.f32(SVector{3, Float32}(xu), 0.0f0) :
    (ef.itd ? ef.raw(xu, 0.0f0) : ef.raw(xu))

@inline (ef::ErasedFunction)(xu, t) = ef.itd ? ef.raw(xu, t) : ef.raw(xu)
@inline (ef::ErasedFunction)(xu) = ef.itd ? ef.raw(xu, 0.0) : ef.raw(xu)

function ErasedFunction(f, itd::Bool = is_time_dependent(f))
    f_2arg = itd ? f : (r, t) -> f(r)

    f64 = FunctionWrapper{SVector{3, Float64}, Tuple{SVector{3, Float64}, Float64}}(
        (r, t) -> SVector{3, Float64}(f_2arg(r, t))
    )
    f32 = FunctionWrapper{SVector{3, Float32}, Tuple{SVector{3, Float32}, Float32}}(
        (r, t) -> SVector{3, Float32}(f_2arg(r, t))
    )
    return ErasedFunction(f64, f32, f, itd)
end

ErasedFunction(ef::ErasedFunction, itd::Bool = ef.itd) = ef

"""
    ErasedField{itd}

Type alias for a `Field` wrapping an `ErasedFunction`.
"""
const ErasedField{itd} = Field{itd, ErasedFunction}

function ErasedField(f::Function)
    itd = is_time_dependent(f)
    return Field{itd, ErasedFunction}(ErasedFunction(f, itd))
end

function ErasedField(f::Field{itd}) where {itd}
    fn = getfield(f, :field_function)
    fn isa ErasedFunction && return f
    return Field{itd, ErasedFunction}(ErasedFunction(fn, itd))
end

ErasedField(f::Field{itd, ErasedFunction}) where {itd} = f
ErasedField(f::AbstractFieldInterpolator) =
    Field{false, ErasedFunction}(ErasedFunction(f, false))

@inline (f::Field{false, ErasedFunction})(xu::SVector{3, Float64}, t::Float64) =
    getfield(f, :field_function).f64(xu, t)
@inline (f::Field{true, ErasedFunction})(xu::SVector{3, Float64}, t::Float64) =
    getfield(f, :field_function).f64(xu, t)
@inline (f::Field{false, ErasedFunction})(xu::SVector{3, Float32}, t::Float32) =
    getfield(f, :field_function).f32(xu, t)
@inline (f::Field{true, ErasedFunction})(xu::SVector{3, Float32}, t::Float32) =
    getfield(f, :field_function).f32(xu, t)

@inline (f::Field{false, ErasedFunction})(xu::SVector{3, Float64}) =
    getfield(f, :field_function).f64(xu, 0.0)
@inline (f::Field{true, ErasedFunction})(xu::SVector{3, Float64}) =
    throw(ArgumentError("Time-dependent field function must have a time argument."))
@inline (f::Field{false, ErasedFunction})(xu::SVector{3, Float32}) =
    getfield(f, :field_function).f32(xu, 0.0f0)
@inline (f::Field{true, ErasedFunction})(xu::SVector{3, Float32}) =
    throw(ArgumentError("Time-dependent field function must have a time argument."))

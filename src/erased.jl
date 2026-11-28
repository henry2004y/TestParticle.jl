using FunctionWrappers: FunctionWrapper

"""
    ErasedFunction

Type-erased callable for 3D vector fields wrapping arbitrary user functions or
interpolators. Holds precompiled `FunctionWrapper` handles for both `Float64` and
`Float32` inputs (supporting both 3-component and full 6-component state inputs),
along with the raw callable for autodiff (`ForwardDiff.Dual`) and custom types.
"""
struct ErasedFunction
    f64::FunctionWrapper{SVector{3, Float64}, Tuple{SVector{3, Float64}, Float64}}
    f32::FunctionWrapper{SVector{3, Float32}, Tuple{SVector{3, Float32}, Float32}}
    f64_full::FunctionWrapper{SVector{3, Float64}, Tuple{SVector{6, Float64}, Float64}}
    f32_full::FunctionWrapper{SVector{3, Float32}, Tuple{SVector{6, Float32}, Float32}}
    raw::Any
    itd::Bool
end

_wrappers(ef::ErasedFunction, ::Type{Float64}) = (ef.f64, ef.f64_full)
_wrappers(ef::ErasedFunction, ::Type{Float32}) = (ef.f32, ef.f32_full)

@inline function (ef::ErasedFunction)(
        xu::AbstractVector{T}, t::T = zero(T)
    )::SVector{3, T} where {T <: Union{Float32, Float64}}
    f3, f6 = _wrappers(ef, T)
    n = length(xu)
    n == 3 && return f3(convert(SVector{3, T}, xu), t)
    n == 6 && return f6(convert(SVector{6, T}, xu), t)
    return SVector{3, T}(ef.itd ? ef.raw(xu, t) : ef.raw(xu))
end

# Autodiff (`ForwardDiff.Dual`) and other element types bypass the wrappers.
@inline (ef::ErasedFunction)(xu, t) = ef.itd ? ef.raw(xu, t) : ef.raw(xu)
@inline (ef::ErasedFunction)(xu) = ef.itd ? ef.raw(xu, 0.0) : ef.raw(xu)

# `f` sees the full state `[x; v]` when it accepts six components, else only `x`.
function _wrap(f, itd, ::Type{T}) where {T}
    f2 = itd ? f : (r, t) -> f(r)
    f3 = FunctionWrapper{SVector{3, T}, Tuple{SVector{3, T}, T}}(
        (r, t) -> SVector{3, T}(f2(r, t))
    )
    takes_full = itd ? applicable(f, zero(SVector{6, T}), zero(T)) :
        applicable(f, zero(SVector{6, T}))
    full = takes_full ? ((y, t) -> SVector{3, T}(f2(y, t))) :
        ((y, t) -> SVector{3, T}(f2(SVector{3, T}(y[1], y[2], y[3]), t)))
    return f3, FunctionWrapper{SVector{3, T}, Tuple{SVector{6, T}, T}}(full)
end

function ErasedFunction(f, itd::Bool = is_time_dependent(f))
    f64, f64_full = _wrap(f, itd, Float64)
    f32, f32_full = _wrap(f, itd, Float32)
    return ErasedFunction(f64, f32, f64_full, f32_full, f, itd)
end

ErasedFunction(ef::ErasedFunction, itd::Bool = ef.itd) = ef

"""
    ErasedField{itd}

Type alias for a `Field` wrapping an `ErasedFunction`.
"""
const ErasedField{itd} = Field{itd, ErasedFunction}

function ErasedField(f)
    itd = is_time_dependent(f)
    return Field{itd, ErasedFunction}(ErasedFunction(f, itd))
end

ErasedField(f::Field{itd}) where {itd} =
    Field{itd, ErasedFunction}(ErasedFunction(getfield(f, :field_function), itd))

function Adapt.adapt_structure(to, f::Field{itd, ErasedFunction}) where {itd}
    raw = getfield(f, :field_function).raw
    adapted_raw = Adapt.adapt(to, raw)
    return Field{itd, typeof(adapted_raw)}(adapted_raw)
end

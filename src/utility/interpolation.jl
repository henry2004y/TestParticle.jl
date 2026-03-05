# Field interpolations using the Interpolations.jl backend.

"""
    TupleCallAdaptor{T}

Wraps an Interpolations.jl interpolant so it can be called with a single
tuple argument, matching the call convention used by `FieldInterpolator`.
"""
struct TupleCallAdaptor{T}
    itp::T
end

Adapt.adapt_structure(to, a::TupleCallAdaptor) = TupleCallAdaptor(Adapt.adapt(to, a.itp))

@inline (a::TupleCallAdaptor{T})(t::NTuple{1}) where {T} = a.itp(t[1])
@inline (a::TupleCallAdaptor{T})(t::NTuple{2}) where {T} = a.itp(t[1], t[2])
@inline (a::TupleCallAdaptor{T})(t::NTuple{3}) where {T} = a.itp(t[1], t[2], t[3])
@inline (a::TupleCallAdaptor{T})(t::AbstractVector) where {T} = a.itp(t[1], t[2], t[3])

function _get_bspline(order::Int, periodic::Bool)
    gt = Interpolations.OnCell()
    interp_type = if order == 1
        Interpolations.Linear
    elseif order == 2
        Interpolations.Quadratic
    elseif order == 3
        Interpolations.Cubic
    else
        throw(ArgumentError("Unsupported interpolation order!"))
    end

    if periodic
        return Interpolations.BSpline(interp_type(Interpolations.Periodic(gt)))
    else
        if interp_type == Interpolations.Linear
            return Interpolations.BSpline(Interpolations.Linear())
        else
            return Interpolations.BSpline(interp_type(Interpolations.Flat(gt)))
        end
    end
end

function _get_interp_object(A, order::Int, bc::Int)
    bspline = _get_bspline(order, bc == 2)

    bctype = if bc == 1
        if eltype(A) <: SVector
            SVector{3, eltype(eltype(A))}(NaN, NaN, NaN)
        else
            eltype(A)(NaN)
        end
    elseif bc == 2
        Interpolations.Periodic()
    else
        Interpolations.Flat()
    end

    return Interpolations.extrapolate(Interpolations.interpolate(A, bspline), bctype)
end

function build_interpolator(
        b::InterpolationsBackend, ::Type{<:CartesianGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 1
    ) where {T}
    @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(b, CartesianGrid, As, gridx, gridy, gridz, order, bc)
end

function build_interpolator(
        ::InterpolationsBackend, ::Type{<:CartesianGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 1
    ) where {T}
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    end
    itp = _get_interp_object(A, order, bc)
    interp = Interpolations.scale(itp, gridx, gridy, gridz)
    return FieldInterpolator(TupleCallAdaptor(interp))
end

function build_interpolator(
        b::InterpolationsBackend, ::Type{<:RectilinearGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 1
    ) where {T}
    @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(b, RectilinearGrid, As, gridx, gridy, gridz, order, bc)
end

function build_interpolator(
        ::InterpolationsBackend, ::Type{<:RectilinearGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 1
    ) where {T}
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    end
    if order != 1
        throw(ArgumentError("RectilinearGrid (CartesianNonUniform) only supports order=1 (Linear) interpolation."))
    end

    bctype = if bc == 1
        if T <: SVector
            fill(eltype(T)(NaN), 3)
        else
            T(NaN)
        end
    elseif bc == 2
        Interpolations.Periodic()
    else
        Interpolations.Flat()
    end

    itp = Interpolations.extrapolate(
        Interpolations.interpolate!((gridx, gridy, gridz), A, Interpolations.Gridded(Interpolations.Linear())),
        bctype
    )
    return FieldInterpolator(TupleCallAdaptor(itp))
end

function build_interpolator(
        b::InterpolationsBackend, ::Type{<:StructuredGrid}, A::AbstractArray{T, 4},
        gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1
    ) where {T}
    @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(b, StructuredGrid, As, gridr, gridθ, gridϕ, order, bc)
end

function build_interpolator(
        ::InterpolationsBackend, ::Type{<:StructuredGrid}, A::AbstractArray{T, 3},
        gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1
    ) where {T}
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    end
    r_min, r_max = extrema(gridr)
    θ_min, θ_max = extrema(gridθ)
    ϕ_min, ϕ_max = extrema(gridϕ)

    @assert r_min >= 0 "r must be non-negative."
    @assert θ_min >= 0 && θ_max <= π "θ must be within [0, π]."
    @assert ϕ_min >= 0 && ϕ_max <= 2π "ϕ must be within [0, 2π]."

    has_0 = isapprox(ϕ_min, 0, atol = 1.0e-5)
    has_2pi = isapprox(ϕ_max, 2π, atol = 1.0e-5)

    ϕ_bc = if has_0 && has_2pi
        Interpolations.Periodic(Interpolations.OnGrid())
    else
        Interpolations.Periodic(Interpolations.OnCell())
    end

    bctype = (Interpolations.Flat(), Interpolations.Flat(), ϕ_bc)

    if order == 1
        itp = Interpolations.extrapolate(
            Interpolations.interpolate!(
                (gridr, gridθ, gridϕ), A,
                Interpolations.Gridded(Interpolations.Linear())
            ),
            bctype
        )
    else
        interp_type = if order == 2
            Interpolations.Quadratic
        elseif order == 3
            Interpolations.Cubic
        else
            throw(ArgumentError("Unsupported interpolation order!"))
        end
        itp_type = (
            Interpolations.BSpline(interp_type(Interpolations.Flat(Interpolations.OnCell()))),
            Interpolations.BSpline(interp_type(Interpolations.Flat(Interpolations.OnCell()))),
            Interpolations.BSpline(interp_type(ϕ_bc)),
        )
        itp_obj = Interpolations.extrapolate(Interpolations.interpolate(A, itp_type), bctype)
        itp = Interpolations.scale(itp_obj, gridr, gridθ, gridϕ)
    end

    return SphericalFieldInterpolator(TupleCallAdaptor(itp))
end

function build_interpolator(
        ::InterpolationsBackend, ::Type{<:CartesianGrid}, A,
        gridx::AbstractVector, gridy::AbstractVector, order::Int = 1, bc::Int = 2
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 2 "Inconsistent 2D force field and grid! Expected 2D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 3 "Inconsistent 2D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    itp = _get_interp_object(As, order, bc)
    interp = Interpolations.scale(itp, gridx, gridy)
    return FieldInterpolator2D(TupleCallAdaptor(interp))
end

function build_interpolator(
        ::InterpolationsBackend, ::Type{<:CartesianGrid}, A, gridx::AbstractVector,
        order::Int = 1, bc::Int = 3; dir = 1
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 1 "Inconsistent 1D force field and grid! Expected 1D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 2 "Inconsistent 1D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    itp = _get_interp_object(As, order, bc)
    interp = Interpolations.scale(itp, gridx)
    return FieldInterpolator1D(TupleCallAdaptor(interp), dir)
end

# Field interpolations.

"""
Type for grid.
"""
abstract type Grid end
"""
Cartesian grid.
"""
struct Cartesian <: Grid end
"""
Spherical grid with uniform r, θ and ϕ.
"""
struct Spherical <: Grid end
"""
Spherical grid with non-uniform r and uniform θ, ϕ.
"""
struct SphericalNonUniformR <: Grid end

function getinterp(A, grid1, grid2, grid3, args...)
   getinterp(Cartesian(), A, grid1, grid2, grid3, args...)
end
getinterp(A, grid1, grid2, args...) = getinterp(Cartesian(), A, grid1, grid2, args...)
getinterp(A, grid1, args...) = getinterp(Cartesian(), A, grid1, args...)

function getinterp_scalar(A, grid1, grid2, grid3, args...)
   getinterp_scalar(Cartesian(), A, grid1, grid2, grid3, args...)
end

"""
     getinterp(::Grid, A, gridx, gridy, gridz, order::Int=1, bc::Int=1)

Return a function for interpolating field array `A` on the grid given by `gridx`, `gridy`,
and `gridz`.

# Arguments

  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic, 3 -> Flat.
  - `dir::Int`: 1/2/3, representing x/y/z direction.
"""
function getinterp(gridtype::Cartesian, A, gridx, gridy, gridz, order::Int = 1, bc::Int = 1)
   if eltype(A) <: SVector
      @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
   else
      @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
   end
   return get_interpolator(gridtype, A, gridx, gridy, gridz, order, bc)
end

function getinterp(
      gridtype::Union{Spherical, SphericalNonUniformR}, A, gridr, gridθ, gridϕ,
      order::Int = 1, bc::Int = 3)
   if eltype(A) <: SVector
      @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
   else
      @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
   end
   return get_interpolator(gridtype, A, gridr, gridθ, gridϕ, order, bc)
end

function getinterp(::Cartesian, A, gridx, gridy, order::Int = 1, bc::Int = 2)
   if eltype(A) <: SVector
      @assert ndims(A) == 2 "Inconsistent 2D force field and grid! Expected 2D array of SVectors."
      As = A
   else
      @assert size(A, 1) == 3 && ndims(A) == 3 "Inconsistent 2D force field and grid!"
      As = reinterpret(reshape, SVector{3, eltype(A)}, A)
   end

   itp = _get_interp_object(As, order, bc)

   interp = scale(itp, gridx, gridy)

   # Return field value at a given location.
   function get_field(xu)
      r = @view xu[1:2]

      return interp(r...)
   end

   return get_field
end

function getinterp(::Cartesian, A, gridx, order::Int = 1, bc::Int = 3; dir = 1)
   if eltype(A) <: SVector
      @assert ndims(A) == 1 "Inconsistent 1D force field and grid! Expected 1D array of SVectors."
      As = A
   else
      @assert size(A, 1) == 3 && ndims(A) == 2 "Inconsistent 1D force field and grid!"
      As = reinterpret(reshape, SVector{3, eltype(A)}, A)
   end

   itp = _get_interp_object(As, order, bc)

   interp = scale(itp, gridx)

   # Return field value at a given location.
   function get_field(xu)
      r = xu[dir]

      return interp(r)
   end

   return get_field
end

function _get_bspline(order::Int, periodic::Bool)
   gt = OnCell()

   interp_type = if order == 1
      Linear
   elseif order == 2
      Quadratic
   elseif order == 3
      Cubic
   else
      throw(ArgumentError("Unsupported interpolation order!"))
   end

   if periodic
      return BSpline(interp_type(Periodic(gt)))
   else
      # For non-periodic, Linear() is special as it doesn't take an argument.
      if interp_type == Linear
         return BSpline(Linear())
      else
         return BSpline(interp_type(Flat(gt)))
      end
   end
end

function _getinterp(Ax, Ay, Az, order::Int, bc::Int)
   itpx = _get_interp_object(Ax, order, bc)
   itpy = _get_interp_object(Ay, order, bc)
   itpz = _get_interp_object(Az, order, bc)

   itpx, itpy, itpz
end

function _getinterp(gridtype::Spherical, Ax, Ay, Az, order::Int, bc::Int)
   itpr = _get_interp_object(gridtype, Ax, order, bc)
   itpθ = _get_interp_object(gridtype, Ay, order, bc)
   itpϕ = _get_interp_object(gridtype, Az, order, bc)

   itpr, itpθ, itpϕ
end

"""
     getinterp_scalar(::Grid, A, gridx, gridy, gridz, order::Int=1, bc::Int=1)

Return a function for interpolating scalar array `A` on the grid given by `gridx`, `gridy`,
and `gridz`. Currently only 3D arrays are supported.

# Arguments

  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic, 3 -> Flat.
  - `dir::Int`: 1/2/3, representing x/y/z direction.
"""
function getinterp_scalar(
      gridtype::Cartesian, A, gridx, gridy, gridz, order::Int = 1, bc::Int = 1)
   return get_interpolator(gridtype, A, gridx, gridy, gridz, order, bc)
end

function getinterp_scalar(gridtype::Union{Spherical, SphericalNonUniformR}, A,
      gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 3)
   return get_interpolator(gridtype, A, gridr, gridθ, gridϕ, order, bc)
end

"""
     get_interpolator(A, gridx, gridy, gridz, order::Int=1, bc::Int=1)
     get_interpolator(gridtype, A, grid1, grid2, grid3, order::Int=1, bc::Int=1)

Return a function for interpolating field array `A` on the grid given by `gridx`, `gridy`,
and `gridz`.

# Arguments

  - `gridtype`: `Cartesian`, `Spherical` or `SphericalNonUniformR`.
  - `A`: field array. For vector field, the first dimension should be 3.
  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic, 3 -> Flat.
"""
function get_interpolator(::Cartesian, A::AbstractArray{T, 4},
      gridx, gridy, gridz, order::Int = 1, bc::Int = 1) where T
   As = reinterpret(reshape, SVector{3, T}, A)
   return get_interpolator(Cartesian(), As, gridx, gridy, gridz, order, bc)
end

function get_interpolator(::Cartesian, A::AbstractArray{T, 3},
      gridx, gridy, gridz, order::Int = 1, bc::Int = 1) where T
   itp = _get_interp_object(A, order, bc)

   interp = scale(itp, gridx, gridy, gridz)

   # Return field value at a given location.
   function get_field(xu)
      r = @view xu[1:3]

      return interp(r...)
   end

   return get_field
end

function get_interpolator(gridtype::Union{Spherical, SphericalNonUniformR},
      A::AbstractArray{T, 4}, gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1) where T
   As = reinterpret(reshape, SVector{3, T}, A)
   return get_interpolator(gridtype, As, gridr, gridθ, gridϕ, order, bc)
end

function _create_spherical_vector_field_interpolator(interpr, interpθ, interpϕ)
   function get_field(xu)
      r_val, θ_val, ϕ_val = cart2sph(xu)

      Br = interpr(r_val, θ_val, ϕ_val)
      Bθ = interpθ(r_val, θ_val, ϕ_val)
      Bϕ = interpϕ(r_val, θ_val, ϕ_val)

      Bvec = sph_to_cart_vector(Br, Bθ, Bϕ, θ_val, ϕ_val)

      return Bvec
   end
   return get_field
end

function _create_spherical_vector_field_interpolator(itp)
   function get_field(xu)
      r_val, θ_val, ϕ_val = cart2sph(xu)

      B_local = itp(r_val, θ_val, ϕ_val)
      Br, Bθ, Bϕ = B_local

      Bvec = sph_to_cart_vector(Br, Bθ, Bϕ, θ_val, ϕ_val)

      return Bvec
   end
   return get_field
end

function get_interpolator(gridtype::Union{Spherical, SphericalNonUniformR},
      A::AbstractArray{T, 3}, gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1) where T
   if gridtype isa Spherical
      itp_unscaled = _get_interp_object(gridtype, A, order, bc)
      itp = scale(itp_unscaled, gridr, gridθ, gridϕ)
   else # SphericalNonUniformR
      #TODO: Respect the passed boundary conditions!
      bctype = (Flat(), Flat(), Periodic())
      itp = extrapolate(interpolate((gridr, gridθ, gridϕ), A, Gridded(Linear())), bctype)
   end

   if T <: SVector
      return _create_spherical_vector_field_interpolator(itp)
   else
      return _create_spherical_scalar_field_interpolator(itp)
   end
end

function _create_spherical_scalar_field_interpolator(interp)
   function get_field(xu)
      r_val, θ_val, ϕ_val = cart2sph(xu[1], xu[2], xu[3])
      return interp(r_val, θ_val, ϕ_val)
   end
   return get_field
end

function _get_interp_object(A, order::Int, bc::Int)
   bspline = _get_bspline(order, bc == 2)

   bctype = if bc == 1
      NaN
   elseif bc == 2
      Periodic()
   else
      Flat()
   end

   itp = extrapolate(interpolate(A, bspline), bctype)
end

function _get_interp_object(::Spherical, A, order::Int, bc::Int)
   bspline_r = _get_bspline(order, false)
   bspline_θ = _get_bspline(order, false)
   bspline_ϕ = _get_bspline(order, true)

   itp_type = (bspline_r, bspline_θ, bspline_ϕ)

   bctype = if bc == 1
      (NaN, NaN, NaN)
   elseif bc == 2
      (Periodic(), Periodic(), Periodic())
   else
      (Flat(), Flat(), Periodic()) # Default to periodic in ϕ
   end

   itp = extrapolate(interpolate(A, itp_type), bctype)
end

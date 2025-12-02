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
   @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
   return get_interpolator(gridtype, A, gridx, gridy, gridz, order, bc)
end

function getinterp(
      gridtype::Union{Spherical, SphericalNonUniformR}, A, gridr, gridθ, gridϕ,
      order::Int = 1, bc::Int = 3)
   @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
   return get_interpolator(gridtype, A, gridr, gridθ, gridϕ, order, bc)
end

function getinterp(::Cartesian, A, gridx, gridy, order::Int = 1, bc::Int = 2)
   @assert size(A, 1) == 3 && ndims(A) == 3 "Inconsistent 2D force field and grid!"

   Ax = @view A[1, :, :]
   Ay = @view A[2, :, :]
   Az = @view A[3, :, :]

   itpx, itpy, itpz = _getinterp(Ax, Ay, Az, order, bc)

   interpx = scale(itpx, gridx, gridy)
   interpy = scale(itpy, gridx, gridy)
   interpz = scale(itpz, gridx, gridy)

   # Return field value at a given location.
   function get_field(xu)
      r = @view xu[1:2]

      return SA[interpx(r...), interpy(r...), interpz(r...)]
   end

   return get_field
end

function getinterp(::Cartesian, A, gridx, order::Int = 1, bc::Int = 3; dir = 1)
   @assert size(A, 1) == 3 && ndims(A) == 2 "Inconsistent 1D force field and grid!"

   Ax = @view A[1, :]
   Ay = @view A[2, :]
   Az = @view A[3, :]

   itpx, itpy, itpz = _getinterp(Ax, Ay, Az, order, bc)

   interpx = scale(itpx, gridx)
   interpy = scale(itpy, gridx)
   interpz = scale(itpz, gridx)

   # Return field value at a given location.
   function get_field(xu)
      r = xu[dir]

      return SA[interpx(r), interpy(r), interpz(r)]
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
   itpx, itpy,
   itpz = _getinterp(
      view(A,1,:,:,:), view(A,2,:,:,:), view(A,3,:,:,:), order, bc)

   interpx = scale(itpx, gridx, gridy, gridz)
   interpy = scale(itpy, gridx, gridy, gridz)
   interpz = scale(itpz, gridx, gridy, gridz)

   # Return field value at a given location.
   function get_field(xu)
      r = @view xu[1:3]

      return SA[interpx(r...), interpy(r...), interpz(r...)]
   end

   return get_field
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
   is_periodic = false

   if gridtype isa Spherical
      dϕ = gridϕ[2] - gridϕ[1]
      if isapprox(gridϕ[end] + dϕ, 2π)
         gridϕ = range(gridϕ[1], step=dϕ, length=length(gridϕ)+1)
         A = cat(A, selectdim(A, 4, 1), dims=4)
      end
      is_periodic = true

      Ar = @view A[1, :, :, :]
      Aθ = @view A[2, :, :, :]
      Aϕ = @view A[3, :, :, :]

      itpr_u, itpθ_u, itpϕ_u = _getinterp(gridtype, Ar, Aθ, Aϕ, order, bc)

      interpr = scale(itpr_u, gridr, gridθ, gridϕ)
      interpθ = scale(itpθ_u, gridr, gridθ, gridϕ)
      interpϕ = scale(itpϕ_u, gridr, gridθ, gridϕ)
   else # SphericalNonUniformR
      Ar = @view A[1, :, :, :]
      Aθ = @view A[2, :, :, :]
      Aϕ = @view A[3, :, :, :]
      if order != 1
         throw(ArgumentError("Only linear interpolation is supported for non-uniform spherical grids!"))
      end
      grid = (gridr, gridθ, gridϕ)
      bctype_r, bctype_θ, bctype_ϕ = if bc == 1
         (NaN, NaN, NaN)
      elseif bc == 2
         (Periodic(), Periodic(), Periodic())
      else
         (Flat(), Flat(), Periodic())
      end
      interpr = extrapolate(interpolate(grid, Ar, Gridded(Linear())), bctype_r)
      interpθ = extrapolate(interpolate(grid, Aθ, Gridded(Linear())), bctype_θ)
      interpϕ = extrapolate(interpolate(grid, Aϕ, Gridded(Linear())), bctype_ϕ)
   end

   return _create_spherical_vector_field_interpolator(interpr, interpθ, interpϕ, is_periodic)
end

function _create_spherical_vector_field_interpolator(interpr, interpθ, interpϕ, periodic::Bool=false)
   function get_field(xu)
      r_val, θ_val, ϕ_val = cart2sph(xu)

      if periodic
         ϕ_val = mod2pi(ϕ_val)
      end

      Br = interpr(r_val, θ_val, ϕ_val)
      Bθ = interpθ(r_val, θ_val, ϕ_val)
      Bϕ = interpϕ(r_val, θ_val, ϕ_val)

      Bvec = sph_to_cart_vector(Br, Bθ, Bϕ, θ_val, ϕ_val)

      return Bvec
   end
   return get_field
end

function get_interpolator(gridtype::Union{Spherical, SphericalNonUniformR},
      A::AbstractArray{T, 3}, gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1) where T
   is_periodic = false

   if gridtype isa Spherical
      dϕ = gridϕ[2] - gridϕ[1]
      if isapprox(gridϕ[end] + dϕ, 2π)
         gridϕ = range(gridϕ[1], step=dϕ, length=length(gridϕ)+1)
         A = cat(A, selectdim(A, 3, 1), dims=3)
      end
      is_periodic = true

      itp_unscaled = _get_interp_object(gridtype, A, order, bc)
      itp = scale(itp_unscaled, gridr, gridθ, gridϕ)
   else # SphericalNonUniformR
      #TODO: Respect the passed boundary conditions!
      bctype = (Flat(), Flat(), Periodic())
      itp = extrapolate(interpolate((gridr, gridθ, gridϕ), A, Gridded(Linear())), bctype)
   end
   return _create_spherical_scalar_field_interpolator(itp, is_periodic)
end

function _create_spherical_scalar_field_interpolator(interp, periodic::Bool=false)
   function get_field(xu)
      r_val, θ_val, ϕ_val = cart2sph(xu[1], xu[2], xu[3])
      if periodic
         ϕ_val = mod2pi(ϕ_val)
      end
      return interp(r_val, θ_val, ϕ_val)
   end
   return get_field
end

function _get_interp_object(A, order::Int, bc::Int)
   bspline = _get_bspline(order, bc == 2)
   itp = interpolate(A, bspline)
end

function _get_interp_object(::Spherical, A, order::Int, bc::Int)
   bspline_r = _get_bspline(order, false)
   bspline_θ = _get_bspline(order, false)
   bspline_ϕ = _get_bspline(order, false)

   itp_type = (bspline_r, bspline_θ, bspline_ϕ)

   itp = interpolate(A, itp_type)
end

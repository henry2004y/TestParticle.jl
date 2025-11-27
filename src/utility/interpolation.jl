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
   if bc == 1
      _getinterp_2d_impl(A, gridx, gridy, order, Val(1))
   elseif bc == 2
      _getinterp_2d_impl(A, gridx, gridy, order, Val(2))
   else
      _getinterp_2d_impl(A, gridx, gridy, order, Val(3))
   end
end

function _getinterp_2d_impl(A, gridx, gridy, order, val_bc::Val{BC}) where BC
   Ax = @view A[1, :, :]
   Ay = @view A[2, :, :]
   Az = @view A[3, :, :]

   itpx, itpy, itpz = _getinterp(Ax, Ay, Az, order, val_bc)

   interpx = scale(itpx, gridx, gridy)
   interpy = scale(itpy, gridx, gridy)
   interpz = scale(itpz, gridx, gridy)

   # Return field value at a given location.
   function get_field(xu)
      x, y = xu[1], xu[2]

      return SA[interpx(x, y), interpy(x, y), interpz(x, y)]
   end

   return get_field
end

function getinterp(::Cartesian, A, gridx, order::Int = 1, bc::Int = 3; dir = 1)
   @assert size(A, 1) == 3 && ndims(A) == 2 "Inconsistent 1D force field and grid!"
   if bc == 1
      _getinterp_1d_impl(A, gridx, order, Val(1), dir)
   elseif bc == 2
      _getinterp_1d_impl(A, gridx, order, Val(2), dir)
   else
      _getinterp_1d_impl(A, gridx, order, Val(3), dir)
   end
end

function _getinterp_1d_impl(A, gridx, order, val_bc::Val{BC}, dir) where BC
   Ax = @view A[1, :]
   Ay = @view A[2, :]
   Az = @view A[3, :]

   itpx, itpy, itpz = _getinterp(Ax, Ay, Az, order, val_bc)

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

function _getinterp(Ax, Ay, Az, order::Int, val_bc::Val{BC}) where BC
   itpx = _get_interp_object(Ax, order, val_bc)
   itpy = _get_interp_object(Ay, order, val_bc)
   itpz = _get_interp_object(Az, order, val_bc)

   itpx, itpy, itpz
end

function _getinterp(gridtype::Spherical, Ax, Ay, Az, order::Int, val_bc::Val{BC}) where BC
   itpr = _get_interp_object(gridtype, Ax, order, val_bc)
   itpθ = _get_interp_object(gridtype, Ay, order, val_bc)
   itpϕ = _get_interp_object(gridtype, Az, order, val_bc)

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
function get_interpolator(gridtype::Cartesian, A::AbstractArray{T, 4},
      gridx, gridy, gridz, order::Int = 1, bc::Int = 1) where T
   if bc == 1
      _get_interpolator_3d_vec_impl(gridtype, A, gridx, gridy, gridz, order, Val(1))
   elseif bc == 2
      _get_interpolator_3d_vec_impl(gridtype, A, gridx, gridy, gridz, order, Val(2))
   else
      _get_interpolator_3d_vec_impl(gridtype, A, gridx, gridy, gridz, order, Val(3))
   end
end

function _get_interpolator_3d_vec_impl(::Cartesian, A, gridx, gridy, gridz, order, val_bc::Val{BC}) where BC
   itpx, itpy,
   itpz = _getinterp(
      view(A,1,:,:,:), view(A,2,:,:,:), view(A,3,:,:,:), order, val_bc)

   interpx = scale(itpx, gridx, gridy, gridz)
   interpy = scale(itpy, gridx, gridy, gridz)
   interpz = scale(itpz, gridx, gridy, gridz)

   # Return field value at a given location.
   function get_field(xu)
      x, y, z = xu[1], xu[2], xu[3]

      return SA[interpx(x, y, z), interpy(x, y, z), interpz(x, y, z)]
   end

   return get_field
end

function get_interpolator(gridtype::Cartesian, A::AbstractArray{T, 3},
      gridx, gridy, gridz, order::Int = 1, bc::Int = 1) where T
   if bc == 1
      _get_interpolator_3d_scalar_impl(gridtype, A, gridx, gridy, gridz, order, Val(1))
   elseif bc == 2
      _get_interpolator_3d_scalar_impl(gridtype, A, gridx, gridy, gridz, order, Val(2))
   else
      _get_interpolator_3d_scalar_impl(gridtype, A, gridx, gridy, gridz, order, Val(3))
   end
end

function _get_interpolator_3d_scalar_impl(::Cartesian, A, gridx, gridy, gridz, order, val_bc::Val{BC}) where BC
   itp = _get_interp_object(A, order, val_bc)

   interp = scale(itp, gridx, gridy, gridz)

   # Return field value at a given location.
   function get_field(xu)
      x, y, z = xu[1], xu[2], xu[3]

      return interp(x, y, z)
   end

   return get_field
end

function get_interpolator(gridtype::Union{Spherical, SphericalNonUniformR},
      A::AbstractArray{T, 4}, gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1) where T
   if bc == 1
      _get_interpolator_spherical_vec_impl(gridtype, A, gridr, gridθ, gridϕ, order, Val(1))
   elseif bc == 2
      _get_interpolator_spherical_vec_impl(gridtype, A, gridr, gridθ, gridϕ, order, Val(2))
   else
      _get_interpolator_spherical_vec_impl(gridtype, A, gridr, gridθ, gridϕ, order, Val(3))
   end
end

function _get_interpolator_spherical_vec_impl(gridtype, A, gridr, gridθ, gridϕ, order, val_bc::Val{BC}) where BC
   Ar = @view A[1, :, :, :]
   Aθ = @view A[2, :, :, :]
   Aϕ = @view A[3, :, :, :]

   if gridtype isa Spherical
      itpr_u, itpθ_u, itpϕ_u = _getinterp(gridtype, Ar, Aθ, Aϕ, order, val_bc)

      interpr = scale(itpr_u, gridr, gridθ, gridϕ)
      interpθ = scale(itpθ_u, gridr, gridθ, gridϕ)
      interpϕ = scale(itpϕ_u, gridr, gridθ, gridϕ)
   else # SphericalNonUniformR
      if order != 1
         throw(ArgumentError("Only linear interpolation is supported for non-uniform spherical grids!"))
      end
      grid = (gridr, gridθ, gridϕ)
      bctype_r, bctype_θ, bctype_ϕ = if BC == 1
         (NaN, NaN, NaN)
      elseif BC == 2
         (Periodic(), Periodic(), Periodic())
      else
         (Flat(), Flat(), Periodic())
      end
      interpr = extrapolate(interpolate(grid, Ar, Gridded(Linear())), bctype_r)
      interpθ = extrapolate(interpolate(grid, Aθ, Gridded(Linear())), bctype_θ)
      interpϕ = extrapolate(interpolate(grid, Aϕ, Gridded(Linear())), bctype_ϕ)
   end

   return _create_spherical_vector_field_interpolator(interpr, interpθ, interpϕ)
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

function get_interpolator(gridtype::Union{Spherical, SphericalNonUniformR},
      A::AbstractArray{T, 3}, gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1) where T
   if bc == 1
      _get_interpolator_spherical_scalar_impl(gridtype, A, gridr, gridθ, gridϕ, order, Val(1))
   elseif bc == 2
      _get_interpolator_spherical_scalar_impl(gridtype, A, gridr, gridθ, gridϕ, order, Val(2))
   else
      _get_interpolator_spherical_scalar_impl(gridtype, A, gridr, gridθ, gridϕ, order, Val(3))
   end
end

function _get_interpolator_spherical_scalar_impl(gridtype, A, gridr, gridθ, gridϕ, order, val_bc::Val{BC}) where BC
   if gridtype isa Spherical
      itp_unscaled = _get_interp_object(gridtype, A, order, val_bc)
      itp = scale(itp_unscaled, gridr, gridθ, gridϕ)
   else # SphericalNonUniformR
      #TODO: Respect the passed boundary conditions!
      bctype = (Flat(), Flat(), Periodic())
      itp = extrapolate(interpolate((gridr, gridθ, gridϕ), A, Gridded(Linear())), bctype)
   end
   return _create_spherical_scalar_field_interpolator(itp)
end

function _create_spherical_scalar_field_interpolator(interp)
   function get_field(xu)
      r_val, θ_val, ϕ_val = cart2sph(xu[1], xu[2], xu[3])
      return interp(r_val, θ_val, ϕ_val)
   end
   return get_field
end

function _get_interp_object(A, order::Int, ::Val{BC}) where BC
   bspline = _get_bspline(order, BC == 2)

   bctype = if BC == 1
      NaN
   elseif BC == 2
      Periodic()
   else
      Flat()
   end

   itp = extrapolate(interpolate(A, bspline), bctype)
end

function _get_interp_object(::Spherical, A, order::Int, ::Val{BC}) where BC
   bspline_r = _get_bspline(order, false)
   bspline_θ = _get_bspline(order, false)
   bspline_ϕ = _get_bspline(order, true)

   itp_type = (bspline_r, bspline_θ, bspline_ϕ)

   bctype = if BC == 1
      (NaN, NaN, NaN)
   elseif BC == 2
      (Periodic(), Periodic(), Periodic())
   else
      (Flat(), Flat(), Periodic()) # Default to periodic in ϕ
   end

   itp = extrapolate(interpolate(A, itp_type), bctype)
end

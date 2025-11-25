# Field interpolations.

"Type for grid."
abstract type Grid end
"Cartesian grid."
struct Cartesian <: Grid end
"Spherical grid."
struct Spherical <: Grid end

"""
     getinterp(::Cartesian, A, gridx, gridy, gridz, order::Int=1, bc::Int=1)
     getinterp(::Cartesian, A, gridx, gridy, order::Int=1, bc::Int=1)
     getinterp(::Cartesian, A, gridx, order::Int=1, bc::Int=1, dir::Int=1)
     getinterp(::Spherical, A, gridr, gridθ, gridϕ, order::Int=1, bc::Int=1)

Return a function for interpolating field array `A` on the grid given by `gridx`, `gridy`,
and `gridz`.

# Arguments

  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic, 3 -> Flat.
  - `dir::Int`: 1/2/3, representing x/y/z direction.
"""
function getinterp(::Cartesian, A, gridx, gridy, gridz, order::Int = 1, bc::Int = 1)
   @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"

   Ax = @view A[1, :, :, :]
   Ay = @view A[2, :, :, :]
   Az = @view A[3, :, :, :]

   itpx, itpy, itpz = _getinterp(Ax, Ay, Az, order, bc)

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

function getinterp(::Spherical, A, gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1)
   @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"

   Ar = @view A[1, :, :, :]
   Aθ = @view A[2, :, :, :]
   Aϕ = @view A[3, :, :, :]

   itpr, itpθ, itpϕ = _getinterp(Ar, Aθ, Aϕ, order, bc)

   interpr = scale(itpr, gridr, gridθ, gridϕ)
   interpθ = scale(itpθ, gridr, gridθ, gridϕ)
   interpϕ = scale(itpϕ, gridr, gridθ, gridϕ)

   # Return field value at a given location.
   function get_field(xu)
      r, ϕ, θ = cart2sph(xu[1], xu[2], xu[3])

      Br = interpr(r, θ, ϕ)
      Bθ = interpθ(r, θ, ϕ)
      Bϕ = interpϕ(r, θ, ϕ)

      Bx, By, Bz = sph_to_cart_vector(Br, Bϕ, Bθ, ϕ, θ)

      return SA[Bx, By, Bz]
   end

   return get_field
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

function _getinterp(Ax, Ay, Az, order::Int, bc::Int)
   gt = OnCell()

   if bc == 1
      if order == 1
         itpx = extrapolate(interpolate(Ax, BSpline(Linear())), NaN)
         itpy = extrapolate(interpolate(Ay, BSpline(Linear())), NaN)
         itpz = extrapolate(interpolate(Az, BSpline(Linear())), NaN)
      elseif order == 2
         itpx = extrapolate(interpolate(Ax, BSpline(Quadratic(Flat(gt)))), NaN)
         itpy = extrapolate(interpolate(Ay, BSpline(Quadratic(Flat(gt)))), NaN)
         itpz = extrapolate(interpolate(Az, BSpline(Quadratic(Flat(gt)))), NaN)
      elseif order == 3
         itpx = extrapolate(interpolate(Ax, BSpline(Cubic(Flat(gt)))), NaN)
         itpy = extrapolate(interpolate(Ay, BSpline(Cubic(Flat(gt)))), NaN)
         itpz = extrapolate(interpolate(Az, BSpline(Cubic(Flat(gt)))), NaN)
      end
   elseif bc == 2
      bctype = Periodic()
      if order == 1
         itpx = extrapolate(interpolate(Ax, BSpline(Linear(Periodic(gt)))), bctype)
         itpy = extrapolate(interpolate(Ay, BSpline(Linear(Periodic(gt)))), bctype)
         itpz = extrapolate(interpolate(Az, BSpline(Linear(Periodic(gt)))), bctype)
      elseif order == 2
         itpx = extrapolate(interpolate(Ax, BSpline(Quadratic(Periodic(gt)))), bctype)
         itpy = extrapolate(interpolate(Ay, BSpline(Quadratic(Periodic(gt)))), bctype)
         itpz = extrapolate(interpolate(Az, BSpline(Quadratic(Periodic(gt)))), bctype)
      elseif order == 3
         itpx = extrapolate(interpolate(Ax, BSpline(Cubic(Periodic(gt)))), bctype)
         itpy = extrapolate(interpolate(Ay, BSpline(Cubic(Periodic(gt)))), bctype)
         itpz = extrapolate(interpolate(Az, BSpline(Cubic(Periodic(gt)))), bctype)
      end
   else
      bctype = Flat()
      if order == 1
         itpx = extrapolate(interpolate(Ax, BSpline(Linear())), bctype)
         itpy = extrapolate(interpolate(Ay, BSpline(Linear())), bctype)
         itpz = extrapolate(interpolate(Az, BSpline(Linear())), bctype)
      elseif order == 2
         itpx = extrapolate(interpolate(Ax, BSpline(Quadratic(Flat(gt)))), bctype)
         itpy = extrapolate(interpolate(Ay, BSpline(Quadratic(Flat(gt)))), bctype)
         itpz = extrapolate(interpolate(Az, BSpline(Quadratic(Flat(gt)))), bctype)
      elseif order == 3
         itpx = extrapolate(interpolate(Ax, BSpline(Cubic(Flat(gt)))), bctype)
         itpy = extrapolate(interpolate(Ay, BSpline(Cubic(Flat(gt)))), bctype)
         itpz = extrapolate(interpolate(Az, BSpline(Cubic(Flat(gt)))), bctype)
      end
   end

   itpx, itpy, itpz
end

"""
     getinterp_scalar(::Cartesian, A, gridx, gridy, gridz, order::Int=1, bc::Int=1)
     getinterp_scalar(::Spherical, A, gridr, gridθ, gridϕ, order::Int=1, bc::Int=1)

Return a function for interpolating scalar array `A` on the grid given by `gridx`, `gridy`,
and `gridz`. Currently only 3D arrays are supported.

# Arguments

  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic, 3 -> Flat.
  - `dir::Int`: 1/2/3, representing x/y/z direction.
"""
function getinterp_scalar(::Cartesian, A, gridx, gridy, gridz, order::Int = 1, bc::Int = 1)
   itp = _getinterp_scalar(A, order, bc)

   interp = scale(itp, gridx, gridy, gridz)

   # Return field value at a given location.
   function get_field(xu)
      r = @view xu[1:3]

      return interp(r...)
   end

   return get_field
end

function getinterp_scalar(::Spherical, A, gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1)
   itp = _getinterp_scalar(A, order, bc)

   interp = scale(itp, gridr, gridθ, gridϕ)

   # Return field value at a given location.
   function get_field(xu)
      r, ϕ, θ = cart2sph(xu[1], xu[2], xu[3])

      return interp(r, θ, ϕ)
   end

   return get_field
end

function _getinterp_scalar(A, order::Int, bc::Int)
   gt = OnCell()

   if bc == 1
      if order == 1
         itp = extrapolate(interpolate(A, BSpline(Linear())), NaN)
      elseif order == 2
         itp = extrapolate(interpolate(A, BSpline(Quadratic(Flat(gt)))), NaN)
      elseif order == 3
         itp = extrapolate(interpolate(A, BSpline(Cubic(Flat(gt)))), NaN)
      end
   elseif bc == 2
      bctype = Periodic()
      if order == 1
         itp = extrapolate(interpolate(A, BSpline(Linear(Periodic(gt)))), bctype)
      elseif order == 2
         itp = extrapolate(interpolate(A, BSpline(Quadratic(Periodic(gt)))), bctype)
      elseif order == 3
         itp = extrapolate(interpolate(A, BSpline(Cubic(Periodic(gt)))), bctype)
      end
   else
      bctype = Flat()
      if order == 1
         itp = extrapolate(interpolate(A, BSpline(Linear())), bctype)
      elseif order == 2
         itp = extrapolate(interpolate(A, BSpline(Quadratic(Flat(gt)))), bctype)
      elseif order == 3
         itp = extrapolate(interpolate(A, BSpline(Cubic(Flat(gt)))), bctype)
      end
   end

   itp
end

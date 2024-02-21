# Field interpolations.

"""
    getinterp(A, gridx, gridy, gridz, order::Int=1, bc::Int=1)
    getinterp(A, gridx, gridy, order::Int=1, bc::Int=1)

Return a function for interpolating array `A` on the grid given by `gridx`, `gridy`, and
`gridz`.

# Arguments
- `order::Int=1`: order of interpolation in [1,2,3].
- `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic.
"""
function getinterp(A, gridx, gridy, gridz, order::Int=1, bc::Int=1)
   @assert size(A,1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"

   Ax = @view A[1,:,:,:]
   Ay = @view A[2,:,:,:]
   Az = @view A[3,:,:,:]

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

function getinterp(A, gridx, gridy, order::Int=1, bc::Int=2)
   @assert size(A,1) == 3 && ndims(A) == 3 "Inconsistent 2D force field and grid!"

   Ax = @view A[1,:,:]
   Ay = @view A[2,:,:]
   Az = @view A[3,:,:]

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

function _getinterp(Ax, Ay, Az, order::Int, bc::Int)
   gtype = OnCell()

   if bc == 1
      if order == 1
         itpx = extrapolate(interpolate(Ax, BSpline(Linear())), NaN)
         itpy = extrapolate(interpolate(Ay, BSpline(Linear())), NaN)
         itpz = extrapolate(interpolate(Az, BSpline(Linear())), NaN)
      elseif order == 2
         itpx = extrapolate(interpolate(Ax, BSpline(Quadratic())), NaN)
         itpy = extrapolate(interpolate(Ay, BSpline(Quadratic())), NaN)
         itpz = extrapolate(interpolate(Az, BSpline(Quadratic())), NaN)
      elseif order == 3
         itpx = extrapolate(interpolate(Ax, BSpline(Cubic(Line(gtype)))), NaN)
         itpy = extrapolate(interpolate(Ay, BSpline(Cubic(Line(gtype)))), NaN)
         itpz = extrapolate(interpolate(Az, BSpline(Cubic(Line(gtype)))), NaN)
      end
   else
      bctype = Periodic()
      if order == 1
         itpx = extrapolate(interpolate(Ax, BSpline(Linear(bctype))), bctype)
         itpy = extrapolate(interpolate(Ay, BSpline(Linear(bctype))), bctype)
         itpz = extrapolate(interpolate(Az, BSpline(Linear(bctype))), bctype)
      elseif order == 2
         itpx = extrapolate(interpolate(Ax, BSpline(Quadratic(Periodic(gtype)))), bctype)
         itpy = extrapolate(interpolate(Ay, BSpline(Quadratic(Periodic(gtype)))), bctype)
         itpz = extrapolate(interpolate(Az, BSpline(Quadratic(Periodic(gtype)))), bctype)
      elseif order == 3
         itpx = extrapolate(interpolate(Ax, BSpline(Cubic(Periodic(gtype)))), bctype)
         itpy = extrapolate(interpolate(Ay, BSpline(Cubic(Periodic(gtype)))), bctype)
         itpz = extrapolate(interpolate(Az, BSpline(Cubic(Periodic(gtype)))), bctype)
      end
   end

   itpx, itpy, itpz
end
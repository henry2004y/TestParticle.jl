module TestParticle

using LinearAlgebra: norm,  ×
using Meshes
using Interpolations

export prepare, trace_numeric!, trace_analytic!

include("constants.jl")

"""
    prepare(grid, E, B, species="proton")

Return a tuple consists of particle charge, mass for a prescribed `species` and
interpolated EM field functions. 
"""
function prepare(grid, E, B; species="proton")

   if species == "proton"
      q = qᵢ
      m = mᵢ
   elseif species == "electron"
      q = qₑ
      m = mₑ
   end

   gridmin = coordinates(minimum(grid))
   gridmax = coordinates(maximum(grid))
   Δx = spacing(grid)
   
   gridx = range(gridmin[1], gridmax[1], step=Δx[1])
   gridy = range(gridmin[2], gridmax[2], step=Δx[2])
   gridz = range(gridmin[3], gridmax[3], step=Δx[3])
   
   if size(E,1) == 3
      Ex = @view E[1,:,:,:]
      Ey = @view E[2,:,:,:]
      Ez = @view E[3,:,:,:]
      Bx = @view B[1,:,:,:]
      By = @view B[2,:,:,:]
      Bz = @view B[3,:,:,:]
   end

   itp = interpolate(Ex, BSpline(Cubic(Line(OnGrid()))))
   interpEx = scale(itp, gridx, gridy, gridz)

   itp = interpolate(Ey, BSpline(Cubic(Line(OnGrid()))))
   interpEy = scale(itp, gridx, gridy, gridz)

   itp = interpolate(Ez, BSpline(Cubic(Line(OnGrid()))))
   interpEz = scale(itp, gridx, gridy, gridz)

   itp = interpolate(Bx, BSpline(Cubic(Line(OnGrid()))))
   interpBx = scale(itp, gridx, gridy, gridz)

   itp = interpolate(By, BSpline(Cubic(Line(OnGrid()))))
   interpBy = scale(itp, gridx, gridy, gridz)

   itp = interpolate(Bz, BSpline(Cubic(Line(OnGrid()))))
   interpBz = scale(itp, gridx, gridy, gridz)

   q, m, (interpEx, interpEy, interpEz), (interpBx, interpBy, interpBz)
end

"""
    prepare(E, B, species="proton")

Return a tuple consists of particle charge, mass for a prescribed `species` and
analytic EM field functions. 
"""
function prepare(E, B; species="proton")

   if species == "proton"
      q = qᵢ
      m = mᵢ
   elseif species == "electron"
      q = qₑ
      m = mₑ
   end

   q, m, E, B
end

# ODE equations for charged particle moving in static numerical EM field.
function trace_numeric!(dy, y, p, t)
   q, m, interpE, interpB = p
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*(getE(y, interpE) + y[4:6] × getB(y, interpB))
end

# ODE equations for charged particle moving in static analytical EM field.
function trace_analytic!(dy, y, p, t)
   q, m, E, B = p
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*(E(y) + y[4:6] × (B(y[1:3])))
end

# Return eletric field at a given location.
function getE(xu, interpE)
   x = @view xu[1:3]
   u = @view xu[4:6]

   [interpE[1](x...), interpE[2](x...), interpE[3](x...)]
end

# Return magnetic field at a given location.
function getB(xu, interpB)
   x = @view xu[1:3]
   u = @view xu[4:6]

   [interpB[1](x...), interpB[2](x...), interpB[3](x...)]
end

"Calculates the magnetic field from a dipole with magnetic moment `M` at `r`."
function dipole(rIn, M)
   x, y, z = rIn
   r = sqrt(x^2 + y^2 + z^2)
   Coef = μ₀/(4*π*r^5)

   B = [3*x^2-r^2 3*x*y     3*x*z;
        3*y*x     3*y^2-r^2 3*y*z;
        3*z*x     3*z*y     3*z^2-r^2
       ] * M * Coef
end

end
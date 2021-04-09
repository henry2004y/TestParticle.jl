module TestParticle

using LinearAlgebra: norm,  ×
using Meshes
using Interpolations

export prepare, trace_numeric!, trace_analytic!, trace_analytic_relativistic!

include("utility/utility.jl")

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

   itp = extrapolate(interpolate(Ex,
      BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpEx = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Ey,
      BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpEy = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Ez,
      BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpEz = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Bx,
      BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpBx = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(By,
      BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpBy = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Bz,
      BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpBz = scale(itp, gridx, gridy, gridz)

   q, m, (interpEx, interpEy, interpEz), (interpBx, interpBy, interpBz)
end

"""
    prepare(E, B; species="proton", q=1.0, m=1.0)

Return a tuple consists of particle charge, mass for a prescribed `species` and
analytic EM field functions. Prescribed `species` are "electron" and "proton";
other species can be manually specified with `species="other"`, `q` and `m`.
"""
function prepare(E, B; species="proton", q=1.0, m=1.0)

   if species == "proton"
      q = qᵢ
      m = mᵢ
   elseif species == "electron"
      q = qₑ
      m = mₑ
   end

   q, m, E, B
end

"ODE equations for charged particle moving in static numerical EM field."
function trace_numeric!(dy, y, p, t)
   q, m, interpE, interpB = p
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*(getE(y, interpE) + y[4:6] × getB(y, interpB))
end

"ODE equations for charged particle moving in static analytical EM field."
function trace_analytic!(dy, y, p, t)
   q, m, E, B = p
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*(E(y) + y[4:6] × (B(y[1:3])))
end

"ODE equations for relativistic charged particle moving in static analytical EM
field."
function trace_analytic_relativistic!(dy, y, p, t)
   q, m, E, B = p
   γInv = √(1.0 - (y[4]*y[4] + y[5]*y[5] + y[6]*y[6])/c^2) 
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*γInv*(E(y) + y[4:6] × (B(y[1:3])))
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

end
module TestParticle

using LinearAlgebra: norm,  ×
using Meshes
using Interpolations

export prepare, trace_numeric!, trace_analytic!, trace_analytic_relativistic!
export Proton, Electron, Ion

include("utility/utility.jl")

@enum Species Proton Electron Ion User

function getchargemass(species, q, m)

   if species == Proton
      q = qᵢ
      m = mᵢ
   elseif species == Electron
      q = qₑ
      m = mₑ 
   end
   q, m
end

function makegrid(grid)

   gridmin = coordinates(minimum(grid))
   gridmax = coordinates(maximum(grid))
   Δx = spacing(grid)
   
   gridx = range(gridmin[1], gridmax[1], step=Δx[1])
   gridy = range(gridmin[2], gridmax[2], step=Δx[2])
   gridz = range(gridmin[3], gridmax[3], step=Δx[3])

   gridx, gridy, gridz
end

function getinterp(A, gridx, gridy, gridz)

   @assert size(A,1) == 3 && ndims(A) == 4 "Only support 3D force field!"

   Ax = @view A[1,:,:,:]
   Ay = @view A[2,:,:,:]
   Az = @view A[3,:,:,:]

   itp = extrapolate(interpolate(Ax,
      BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpx = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Ay,
      BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpy = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Az,
      BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpz = scale(itp, gridx, gridy, gridz)

   interpx, interpy, interpz
end

"""
    prepare(grid, E, B; species=Proton)

Return a tuple consists of particle charge, mass for a prescribed `species` and
interpolated EM field functions. 
"""
function prepare(grid::CartesianGrid, E, B; species=Proton, q=1.0, m=1.0)

   q, m = getchargemass(species, q, m)

   gridx, gridy, gridz = makegrid(grid)
   
   interpEx, interpEy, interpEz = getinterp(E, gridx, gridy, gridz)
   interpBx, interpBy, interpBz = getinterp(B, gridx, gridy, gridz)

   q, m, (interpEx, interpEy, interpEz), (interpBx, interpBy, interpBz)
end

"""
    prepare(grid, E, B, F; species=Proton, q=1.0, m=1.0)

Return a tuple consists of particle charge, mass for a prescribed `species`,
analytic EM field functions, and external force `F`.
"""
function prepare(grid::CartesianGrid, E, B, F; species=Proton, q=1.0, m=1.0)

   q, m = getchargemass(species, q, m)

   gridx, gridy, gridz = makegrid(grid)

   interpEx, interpEy, interpEz = getinterp(E, gridx, gridy, gridz)
   interpBx, interpBy, interpBz = getinterp(B, gridx, gridy, gridz)
   interpFx, interpFy, interpFz = getinterp(F, gridx, gridy, gridz)
   
   q, m, (interpEx, interpEy, interpEz), (interpBx, interpBy, interpBz),
   (interpFx, interpFy, interpFz)
end

"""
    prepare(E, B; species=Proton, q=1.0, m=1.0)

Return a tuple consists of particle charge, mass for a prescribed `species` and
analytic EM field functions. Prescribed `species` are `Electron` and `Proton`;
other species can be manually specified with `species=Ion/User`, `q` and `m`.
"""
function prepare(E, B; species=Proton, q=1.0, m=1.0)

   q, m = getchargemass(species, q, m)

   q, m, E, B
end

"""
    prepare(E, B, F; species=Proton, q=1.0, m=1.0)

Return a tuple consists of particle charge, mass for a prescribed `species`,
analytic EM field functions, and external force `F`.
"""
function prepare(E, B, F; species=Proton, q=1.0, m=1.0)
   t = prepare(E, B; species, q, m)
   push!(t, F)
   t
end

"ODE equations for charged particle moving in static numerical EM field."
function trace_numeric!(dy, y, p, t)
   q, m, interpE, interpB = p
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*(getE(y, interpE) + y[4:6] × getB(y, interpB))
end

"ODE equations for charged particle moving in static numerical EM field and
external force field."
function trace_numeric_full!(dy, y, p, t)
   q, m, interpE, interpB, interpF = p
   dy[1:3] = y[4:6]
   dy[4:6] = (q*(getE(y, interpE) + y[4:6] × getB(y, interpB)) + 
      getF(t, interpF)) / m
end

"ODE equations for relativistic charged particle moving in static numerical EM
field."
function trace_numeric_relativistic!(dy, y, p, t)
   q, m, E, B = p
   if y[4]*y[4] + y[5]*y[5] + y[6]*y[6] ≥ c^2
      throw(ArgumentError("Particle faster than the speed of light!"))
   end
   γInv = √(1.0 - (y[4]*y[4] + y[5]*y[5] + y[6]*y[6])/c^2) 
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*γInv*(getE(y, interpE) + y[4:6] × getB(y, interpB))
end

"ODE equations for charged particle moving in static analytical EM field."
function trace_analytic!(dy, y, p, t)
   q, m, E, B = p
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*(E(y) + y[4:6] × (B(y[1:3])))
end

"ODE equations for charged particle moving in static analytical EM field and
external force field."
function trace_analytic_full!(dy, y, p, t)
   q, m, E, B, F = p
   dy[1:3] = y[4:6]
   dy[4:6] = (q*(E(y) + y[4:6] × (B(y[1:3]))) + F) / m
end

"ODE equations for relativistic charged particle moving in static analytical EM
field."
function trace_analytic_relativistic!(dy, y, p, t)
   q, m, E, B = p

   if y[4]*y[4] + y[5]*y[5] + y[6]*y[6] ≥ c^2
      throw(ArgumentError("Particle faster than the speed of light!"))
   end
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

# Return force at a given location.
function getF(xu, interpF)
   x = @view xu[1:3]
   u = @view xu[4:6]

   [interpF[1](x...), interpF[2](x...), interpF[3](x...)]
end

end
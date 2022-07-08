module TestParticle

using LinearAlgebra: norm, ×
using Meshes
using Interpolations
using StaticArrays

export prepare, trace!, trace_relativistic!
export trace, trace_relativistic
export Proton, Electron, Ion, User

include("utility/utility.jl")

"""
Type for the particles, `Proton`, `Electron`, `Ion`, or `User`.
"""
@enum Species Proton Electron Ion User

function getchargemass(species::Species, q, m)

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

   itp = extrapolate(interpolate(Ax, BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpx = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Ay, BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpy = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Az, BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpz = scale(itp, gridx, gridy, gridz)

   # Return field value at a given location.
   function get_field(xu)
      r = @view xu[1:3]

      return SA[interpx(r...), interpy(r...), interpz(r...)]
   end

   return get_field
end

"""
    prepare(grid, E, B; species=Proton) -> (q, m, E, B)

Return a tuple consists of particle charge, mass for a prescribed `species` and interpolated
EM field functions.
"""
function prepare(grid::CartesianGrid, E::TE, B::TB; species::Species=Proton, q=1.0, m=1.0) where {TE, TB}

   q, m = getchargemass(species, q, m)

   gridx, gridy, gridz = makegrid(grid)

   E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz) : E
   B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz) : B

   q, m, E, B
end

"""
    prepare(grid, E, B, F; species=Proton, q=1.0, m=1.0) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, analytic EM field functions, and external force `F`.
"""
function prepare(grid::CartesianGrid, E::TE, B::TB, F::TF; species::Species=Proton, q=1.0, m=1.0) where {TE, TB, TF}

   q, m = getchargemass(species, q, m)

   gridx, gridy, gridz = makegrid(grid)

   E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz) : E
   B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz) : B
   F = TF <: AbstractArray ? getinterp(F, gridx, gridy, gridz) : F

   q, m, E, B, F
end

"""
    prepare(E, B; species=Proton, q=1.0, m=1.0) -> (q, m, E, B)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m` and analytic EM field functions. Prescribed `species` are `Electron` and
`Proton`; other species can be manually specified with `species=Ion/User`, `q` and `m`.
"""
function prepare(E, B; species::Species=Proton, q=1.0, m=1.0)

   q, m = getchargemass(species, q, m)

   q, m, E, B
end

"""
    prepare(E, B, F; species=Proton, q=1.0, m=1.0) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, analytic EM field functions, and external force `F`.
"""
function prepare(E, B, F; species::Species=Proton, q=1.0, m=1.0)
   t = prepare(E, B; species, q, m)
   push!(t, F)
   t
end

"ODE equations for charged particle moving in static EM field."
function trace!(dy, y, p, t)
   q, m, E, B = p
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*(E(y) + y[4:6] × (B(y[1:3])))
end

function trace(y, p, t)
   q, m, E, B = p
   dx, dy, dz = y[4:6]
   dux, duy, duz = q/m*(E(y) + y[4:6] × (B(y[1:3])))
   SVector{6}(dx, dy, dz, dux, duy, duz)
end

"ODE equations for charged particle moving in static EM field and external force
field."
function trace_full!(dy, y, p, t)
   q, m, E, B, F = p
   dy[1:3] = y[4:6]
   dy[4:6] = (q*(E(y) + y[4:6] × (B(y[1:3]))) + F) / m
end

function trace_full(y, p, t)
   q, m, E, B, F = p
   dx, dy, dz = y[4:6]
   dux, duy, duz = (q*(E(y) + y[4:6] × (B(y[1:3]))) + F) / m
   SVector{6}(dx, dy, dz, dux, duy, duz)
end

"ODE equations for relativistic charged particle moving in static EM field."
function trace_relativistic!(dy, y, p, t)
   q, m, E, B = p

   if y[4]*y[4] + y[5]*y[5] + y[6]*y[6] ≥ c^2
      throw(ArgumentError("Particle faster than the speed of light!"))
   end
   γInv = √(1.0 - (y[4]*y[4] + y[5]*y[5] + y[6]*y[6])/c^2)
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*γInv*(E(y) + y[4:6] × (B(y[1:3])))
end

function trace_relativistic(y, p, t)
   q, m, E, B = p

   if y[4]*y[4] + y[5]*y[5] + y[6]*y[6] ≥ c^2
      throw(ArgumentError("Particle faster than the speed of light!"))
   end
   γInv = √(1.0 - (y[4]*y[4] + y[5]*y[5] + y[6]*y[6])/c^2)
   dx, dy, dz = y[4:6]
   dux, duy, duz = q/m*γInv*(E(y) + y[4:6] × (B(y[1:3])))
   SVector{6}(dx, dy, dz, dux, duy, duz)
end

end
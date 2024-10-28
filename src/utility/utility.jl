# Collection of utility functions and commonly used constants.

"Convert from spherical to Cartesian coordinates vector."
function sph2cart(r, ϕ, θ)
   r*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
end

include("constants.jl")
include("current_sheet.jl")
include("dipole.jl")
include("confinement.jl")

"""
    getchargemass(species::Species, q::AbstractFloat, m::AbstractFloat)

Return charge and mass for `species`. if `species = User`, input `q` and `m` are returned.
"""
function getchargemass(species::Species, q::AbstractFloat, m::AbstractFloat)
   if species == Proton
      q = qᵢ
      m = mᵢ
   elseif species == Electron
      q = qₑ
      m = mₑ
   end
   q, m
end

"Return uniform range from 2D/3D CartesianGrid."
function makegrid(grid::CartesianGrid)
   gridmin = coords(minimum(grid))
   gridmax = coords(maximum(grid))
   Δx = spacing(grid)
   dim = paramdim(grid)

   gridx = range(gridmin.x.val, gridmax.x.val, step=Δx[1].val)
   gridy = range(gridmin.y.val, gridmax.y.val, step=Δx[2].val)
   if dim == 3
      gridz = range(gridmin.z.val, gridmax.z.val, step=Δx[3].val)
      return gridx, gridy, gridz
   elseif dim == 2
      return gridx, gridy
   end
end

"""
    set_axes_equal(ax)

Set 3D plot axes to equal scale for Matplotlib.
Make axes of 3D plot have equal scale so that spheres appear as spheres and cubes as cubes.
Required since `ax.axis('equal')` and `ax.set_aspect('equal')` don't work on 3D.
"""
function set_axes_equal(ax)
   limits = zeros(2,3)
   limits[:,1] .= ax.get_xlim3d()
   limits[:,2] .= ax.get_ylim3d()
   limits[:,3] .= ax.get_zlim3d()
   origin = mean(limits, dims=1)
   radius = @. 0.5 * max(abs(limits[2,:] - limits[1,:]))
   x, y, z = origin
   ax.set_xlim3d([x - radius[1], x + radius[1]])
   ax.set_ylim3d([y - radius[2], y + radius[2]])
   ax.set_zlim3d([z - radius[3], z + radius[3]])
end

"""
    get_rotation_matrix(axis::AbstractVector, angle::Real) --> SMatrix{3,3}

Create a rotation matrix for rotating a 3D vector around a unit `axis` by an `angle` in
radians.
Reference: [Rotation matrix from axis and angle](https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)

# Example

```julia
using LinearAlgebra
v = [-0.5, 1.0, 1.0]
v̂ = normalize(v)
θ = deg2rad(-74)
R = get_rotation_matrix(v̂, θ)
```
"""
function get_rotation_matrix(v::AbstractVector{<:AbstractFloat}, θ::Real)
   sinθ, cosθ = sincos(eltype(v)(θ))
   tmp = 1 - cosθ
   m =  @SMatrix [
        cosθ+v[1]^2*tmp         v[1]*v[2]*tmp-v[3]*sinθ v[1]*v[3]*tmp+v[2]*sinθ;
        v[1]*v[2]*tmp+v[3]*sinθ cosθ+v[2]^2*tmp         v[2]*v[3]*tmp-v[1]*sinθ;
        v[1]*v[3]*tmp-v[2]*sinθ v[3]*v[2]*tmp+v[1]*sinθ cosθ+v[3]^2*tmp]
end

"Return the gyrofrequency."
function get_gyrofrequency(B::AbstractFloat=5e-9; q::AbstractFloat=qᵢ, m::AbstractFloat=mᵢ)
   ω = q * B / m
end

"Return the gyroradius."
function get_gyroradius(V::AbstractFloat, B::AbstractFloat; q::AbstractFloat=qᵢ, m::AbstractFloat=mᵢ)
   ω = get_gyrofrequency(B; q, m)
   r = V / ω
end

"Return the gyroperiod."
function get_gyroperiod(B::AbstractFloat=5e-9; q::AbstractFloat=qᵢ, m::AbstractFloat=mᵢ)
   ω = get_gyrofrequency(B; q, m)
   2π / ω
end

"Return velocity from relativistic γv in `sol`."
function get_velocity(sol)
   v = Array{eltype(sol.u[1]), 2}(undef, 3, length(sol))
   for is in axes(v, 2)
      γv = @view sol[4:6, is]
      γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
      v² = γ²v² / (1 + γ²v²/c^2)
      γ = 1 / √(1 - v²/c^2)
      for i in axes(v, 1)
         v[i,is] = γv[i] / γ
      end
   end

   v
end

"Return the energy [eV] from relativistic `sol`."
function get_energy(sol::AbstractODESolution; m=mᵢ, q=qᵢ)
   e = Vector{eltype(sol.u[1])}(undef, length(sol))
   for i in eachindex(e)
      γv = @view sol[4:6, i]
      γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
      v² = γ²v² / (1 + γ²v²/c^2)
      γ = 1 / √(1 - v²/c^2)
      e[i] = (γ-1)*m*c^2/abs(q)
   end

   e
end

"Calculate the energy [eV] of a relativistic particle from γv."
function get_energy(γv; m=mᵢ, q=qᵢ)
   γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
   v² = γ²v² / (1 + γ²v²/c^2)
   γ = 1 / √(1 - v²/c^2)

   (γ-1)*m*c^2/abs(q)
end

"Return velocity magnitude from energy in [eV]."
energy2velocity(Ek; m=mᵢ, q=qᵢ) = c*sqrt(1 - 1/(1+Ek*abs(q)/(m*c^2))^2)
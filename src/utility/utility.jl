# Collection of utility functions and commonly used constants.

using Statistics: mean

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
function makegrid(grid::CartesianGrid{3, T}) where T
   gridmin = coordinates(minimum(grid))
   gridmax = coordinates(maximum(grid))
   Δx = spacing(grid)

   gridx = range(gridmin[1], gridmax[1], step=Δx[1])
   gridy = range(gridmin[2], gridmax[2], step=Δx[2])
   gridz = range(gridmin[3], gridmax[3], step=Δx[3])

   gridx, gridy, gridz
end

function makegrid(grid::CartesianGrid{2, T}) where T
   gridmin = coordinates(minimum(grid))
   gridmax = coordinates(maximum(grid))
   Δx = spacing(grid)

   gridx = range(gridmin[1], gridmax[1], step=Δx[1])
   gridy = range(gridmin[2], gridmax[2], step=Δx[2])

   gridx, gridy
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
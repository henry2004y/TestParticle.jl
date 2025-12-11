# Collection of utility functions and commonly used constants.

"""
Convert from spherical to Cartesian coordinates vector.
"""
function sph2cart(r, θ, ϕ)
   sinθ, cosθ = sincos(θ)
   sinϕ, cosϕ = sincos(ϕ)
   SVector{3}(r * sinθ * cosϕ, r * sinθ * sinϕ, r * cosθ)
end

@inline @inbounds sph2cart(x) = sph2cart(x[1], x[2], x[3])

"""
Convert from Cartesian to spherical coordinates vector.
"""
function cart2sph(x, y, z)
   r = hypot(x, y, z)
   if r == 0
      return SVector{3, eltype(r)}(0, 0, 0)
   end
   θ = acos(z / r)
   ϕ = atan(y, x) |> mod2pi
   return SVector{3}(r, θ, ϕ)
end

@inline @inbounds cart2sph(x) = cart2sph(x[1], x[2], x[3])

"""
Convert a vector from spherical to Cartesian.
"""
function sph_to_cart_vector(vr, vθ, vϕ, θ, ϕ)
   sinθ, cosθ = sincos(θ)
   sinϕ, cosϕ = sincos(ϕ)
   vx = sinθ * cosϕ * vr + cosθ * cosϕ * vθ - sinϕ * vϕ
   vy = sinθ * sinϕ * vr + cosθ * sinϕ * vθ + cosϕ * vϕ
   vz = cosθ * vr - sinθ * vθ
   return SVector{3}(vx, vy, vz)
end

include("constants.jl")
include("loop.jl")
include("confinement.jl")

"""
     getchargemass(species::Species, q, m)

Return charge and mass for `species`.
For `species = Ion`, `q` and `m` are charge and mass numbers. For `species = User`, the input `q` and `m` are returned as is.
"""
function getchargemass(species::Species, q::Real, m::Real)
   if species == Proton
      q = qᵢ
      m = mᵢ
   elseif species == Electron
      q = qₑ
      m = mₑ
   elseif species == Ion
      q *= qᵢ
      m *= mᵢ
   end
   q, m
end

"""
Return uniform range from 2D/3D CartesianGrid.
"""
function makegrid(grid::CartesianGrid)
   gridmin = coords(minimum(grid))
   gridmax = coords(maximum(grid))
   Δx = spacing(grid)
   dim = paramdim(grid)

   gridx = range(gridmin.x.val, gridmax.x.val, step = Δx[1].val)
   gridy = range(gridmin.y.val, gridmax.y.val, step = Δx[2].val)
   if dim == 3
      gridz = range(gridmin.z.val, gridmax.z.val, step = Δx[3].val)
      return gridx, gridy, gridz
   elseif dim == 2
      return gridx, gridy
   elseif dim == 1
      return (gridx, )
   end
end

"""
Return ranges from 2D/3D RectilinearGrid.
"""
function makegrid(grid::RectilinearGrid)
   # Meshes.RectilinearGrid stores coordinates as vectors.
   # Access pattern: grid.x, grid.y, grid.z is typical if properties are exposed.
   # Based on Meshes.jl docs: RectilinearGrid(x, y, z).
   # We need to access these fields. Assuming standard field access works (xyz fields).
   # If not, we might need a specific accessor.
   # Checking current Meshes.jl implementation via GitHub source link in docs:
   # struct RectilinearGrid{M,C,T} <: Grid{M,C,T}
   #   xyz::Tuple
   # end

   if paramdim(grid) == 3
      return grid.xyz[1], grid.xyz[2], grid.xyz[3]
   elseif paramdim(grid) == 2
      return grid.xyz[1], grid.xyz[2]
   end
end

"""
     set_axes_equal(ax)

Set 3D plot axes to equal scale for Matplotlib.
Make axes of 3D plot have equal scale so that spheres appear as spheres and cubes as cubes.
Required since `ax.axis('equal')` and `ax.set_aspect('equal')` don't work on 3D.
"""
function set_axes_equal(ax)
   limits = zeros(2, 3)
   limits[:, 1] .= ax.get_xlim3d()
   limits[:, 2] .= ax.get_ylim3d()
   limits[:, 3] .= ax.get_zlim3d()
   origin = mean(limits, dims = 1)
   radius = @. 0.5 * max(abs(limits[2, :] - limits[1, :]))
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
   m = @SMatrix [cosθ+v[1]^2 * tmp v[1] * v[2] * tmp-v[3] * sinθ v[1] * v[3] * tmp+v[2] * sinθ;
                 v[1] * v[2] * tmp+v[3] * sinθ cosθ+v[2]^2 * tmp v[2] * v[3] * tmp-v[1] * sinθ;
                 v[1] * v[3] * tmp-v[2] * sinθ v[3] * v[2] * tmp+v[1] * sinθ cosθ+v[3]^2 * tmp]
end

"""
    get_gyrofrequency(B=5e-9; q=qᵢ, m=mᵢ)

Return the gyrofrequency [rad/s].

# Arguments

  - `B::AbstractFloat`: Magnetic field magnitude [T]. Default is 5 nT.
  - `q::AbstractFloat`: Charge [C]. Default is proton charge.
  - `m::AbstractFloat`: Mass [kg]. Default is proton mass.
"""
function get_gyrofrequency(
      B::AbstractFloat = 5e-9;
      q::AbstractFloat = qᵢ,
      m::AbstractFloat = mᵢ
)
   ω = q * B / m
end

"""
    get_gyroradius(V, B; q=qᵢ, m=mᵢ)

Return the gyroradius [m].

# Arguments

  - `V::AbstractFloat`: Velocity magnitude [m/s] (usually perpendicular to the magnetic field).
  - `B::AbstractFloat`: Magnetic field magnitude [T].
  - `q::AbstractFloat`: Charge [C]. Default is proton charge.
  - `m::AbstractFloat`: Mass [kg]. Default is proton mass.
"""
function get_gyroradius(
      V::AbstractFloat,
      B::AbstractFloat;
      q::AbstractFloat = qᵢ,
      m::AbstractFloat = mᵢ
)
   ω = get_gyrofrequency(B; q, m)
   r = V / ω
end

"""
    get_gyroperiod(B=5e-9; q=qᵢ, m=mᵢ)

Return the gyroperiod [s].

# Arguments

  - `B::AbstractFloat`: Magnetic field magnitude [T]. Default is 5 nT.
  - `q::AbstractFloat`: Charge [C]. Default is proton charge.
  - `m::AbstractFloat`: Mass [kg]. Default is proton mass.
"""
function get_gyroperiod(
      B::AbstractFloat = 5e-9;
      q::AbstractFloat = qᵢ,
      m::AbstractFloat = mᵢ
)
   ω = get_gyrofrequency(B; q, m)
   2π / ω
end

"""
Return velocity from relativistic γv in `sol`.
"""
function get_velocity(sol)
   v = Array{eltype(sol.u[1]), 2}(undef, 3, length(sol))
   for is in axes(v, 2)
      γv = @view sol[4:6, is]
      γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
      v² = γ²v² / (1 + γ²v² / c^2)
      γ = 1 / √(1 - v² / c^2)
      for i in axes(v, 1)
         v[i, is] = γv[i] / γ
      end
   end

   v
end

"""
Return the energy [eV] from relativistic `sol`.
"""
function get_energy(sol::AbstractODESolution; m = mᵢ, q = qᵢ)
   e = Vector{eltype(sol.u[1])}(undef, length(sol))
   for i in eachindex(e)
      γv = @view sol[4:6, i]
      γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
      v² = γ²v² / (1 + γ²v² / c^2)
      γ = 1 / √(1 - v² / c^2)
      e[i] = (γ - 1) * m * c^2 / abs(q)
   end

   e
end

"""
Calculate the energy [eV] of a relativistic particle from γv.
"""
function get_energy(γv; m = mᵢ, q = qᵢ)
   γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
   v² = γ²v² / (1 + γ²v² / c^2)
   γ = 1 / √(1 - v² / c^2)

   (γ - 1) * m * c^2 / abs(q)
end

"""
Return velocity magnitude from energy in [eV].
"""
energy2velocity(Ek; m = mᵢ, q = qᵢ) = c * sqrt(1 - 1 / (1 + Ek * abs(q) / (m * c^2))^2)
include("random.jl")

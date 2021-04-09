# Collection of utility functions and commonly used constants.

include("constants.jl")
include("current_sheet.jl")
include("dipole.jl")
include("confinement.jl")

"Convert from spherical to Cartesian coordinates vector."
function sph2cart(r, ϕ, θ)
   r*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
end
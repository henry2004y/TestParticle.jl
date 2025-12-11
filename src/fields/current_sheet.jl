# Current sheet model.

"""
     getB_CS_harris(B₀, L)

Return the magnetic field at location `r` near a current sheet with magnetic strength `B₀`
and sheet length `L`. The current sheet is assumed to lie in the z = 0 plane.
"""
function getB_CS_harris(r, B₀, L, Bn = 0.0)
   z = r[3]
   SA[B₀ * tanh(z / L), 0.0, Bn]
end

struct HarrisCurrentSheet{T} <: AbstractField{false}
   B₀::T
   L::T
   Bn::T
end

function (H::HarrisCurrentSheet)(xu)
   z = xu[3]
   SA[H.B₀ * tanh(z / H.L), 0.0, H.Bn]
end

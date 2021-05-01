module CurrentSheet
# Current sheet model.

"""
    getB_CS_harris(B₀, L)

Return the magnetic field at location `r` near a current sheet with magnetic strength `B₀`
and sheet length `L`.
"""
function getB_CS_harris(r, B₀, L)
   x = r[1]
   [0.0, 0.0, B₀ * tanh(x/L)]
end

end
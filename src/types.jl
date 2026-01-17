abstract type AbstractField{itd} <: Function end

"""
Abstract type for tracing solutions.
"""
abstract type AbstractTraceSolution{T, N, S} <: AbstractODESolution{T, N, S} end

"""
Abstract type for velocity distribution functions.
"""
abstract type VDF end

"""
Type for the particles: `Proton`, `Electron`.
"""
struct Species{M, Q}
    m::M
    q::Q
end

Ion(m, q = 1) = Species(m * mᵢ, q * qᵢ)
Ion(; m = 1, q = 1) = Species(m * mᵢ, q * qᵢ)

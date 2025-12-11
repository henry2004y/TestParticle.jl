abstract type AbstractField{itd} <: Function end

"""
Abstract type for tracing solutions.
"""
abstract type AbstractTraceSolution{T, N, S} <: AbstractODESolution{T, N, S} end

"""
Abstract type for velocity distribution functions.
"""
abstract type VDF end

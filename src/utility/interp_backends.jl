abstract type AbstractInterpolationBackend end

"""
    FastInterpolationsBackend

Interpolation backend using FastInterpolations.jl (default).
"""
struct FastInterpolationsBackend <: AbstractInterpolationBackend end

"""
    InterpolationsBackend

Interpolation backend using Interpolations.jl.
"""
struct InterpolationsBackend <: AbstractInterpolationBackend end

module TestParticle

using LinearAlgebra: norm, ×, ⋅
using Meshes: coordinates, spacing, embeddim, CartesianGrid
using Interpolations: interpolate, extrapolate, scale, BSpline, Linear, Quadratic, Cubic,
   Line, OnCell, Periodic
using SciMLBase: AbstractSciMLProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
   BasicEnsembleAlgorithm, EnsembleThreads, EnsembleSerial,
   DEFAULT_SPECIALIZATION, ODEFunction,
   LinearInterpolation
using Distributions: MvNormal
using StaticArrays
using ChunkSplitters
using PrecompileTools: @setup_workload, @compile_workload

export prepare, sample
export trace!, trace_relativistic!, trace_normalized!, trace, trace_relativistic,
   trace_relativistic_normalized!
export Proton, Electron, Ion, User
export Maxwellian, BiMaxwellian
export orbit, monitor
export TraceProblem

"""
Type for the particles, `Proton`, `Electron`, `Ion`, or `User`.
"""
@enum Species Proton Electron Ion User

"""
Abstract type for tracing solutions.
"""
abstract type AbstractTraceSolution{T, N, S} <: AbstractODESolution{T, N, S} end

include("utility/utility.jl")
include("utility/interpolation.jl")
include("sampler.jl")
include("prepare.jl")
include("equations.jl")
include("pusher.jl")

function orbit end

function monitor end

include("precompile.jl")

end
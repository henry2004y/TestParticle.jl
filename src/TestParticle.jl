module TestParticle

using LinearAlgebra: norm, ×, ⋅
using Statistics: mean, normalize
using Interpolations: interpolate, extrapolate, scale, BSpline, Linear, Quadratic, Cubic,
   Line, OnCell, Periodic, Flat
using SciMLBase: AbstractODEProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
   BasicEnsembleAlgorithm, EnsembleThreads, EnsembleSerial,
   DEFAULT_SPECIALIZATION, ODEFunction, 
   LinearInterpolation
using StaticArrays: SVector, @SMatrix, MVector, SA
using Meshes: coords, spacing, embeddim, CartesianGrid
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
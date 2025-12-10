module TestParticle

using LinearAlgebra: norm, ×, ⋅, diag
using Statistics: mean, normalize
using Interpolations: interpolate, interpolate!, extrapolate, scale, BSpline, Linear,
                      Quadratic, Cubic,
                      Line, OnCell, Periodic, Flat, Gridded
using SciMLBase: AbstractODEProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
                 BasicEnsembleAlgorithm, EnsembleThreads, EnsembleSerial,
                 DEFAULT_SPECIALIZATION, ODEFunction, ODEProblem,
                 LinearInterpolation
using StaticArrays: SVector, @SMatrix, MVector, SA, StaticArray
using Meshes: coords, spacing, paramdim, CartesianGrid
using ForwardDiff
using ChunkSplitters
using PrecompileTools: @setup_workload, @compile_workload
using Elliptic
import Tensors
import Base: +, *, /, setindex!, getindex
import LinearAlgebra: ×

export prepare, prepare_gc, sample, get_gc, get_gc_func
export trace!, trace_relativistic!, trace_normalized!, trace_relativistic_normalized!,
       trace, trace_relativistic, trace_relativistic_normalized, trace_gc!, trace_gc_1st!,
       trace_gc_drifts!, trace_gc_flr!, trace_gc_exb!, trace_fieldline!, trace_fieldline
export Proton, Electron, Ion, User
export Maxwellian, BiMaxwellian
export Kappa, BiKappa, SelfSimilar, BiSelfSimilar
export AdaptiveBoris
export get_gyrofrequency,
       get_gyroperiod, get_gyroradius, get_velocity, get_energy,
       energy2velocity
export orbit, monitor
export TraceProblem, Cartesian, CartesianNonUniform, Spherical, SphericalNonUniformR

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
include("gc.jl")
include("equations.jl")
include("pusher.jl")
include("multistep_boris.jl")
include("adaptive_boris.jl")
include("fieldline.jl")

function orbit end

function monitor end

include("precompile.jl")

end

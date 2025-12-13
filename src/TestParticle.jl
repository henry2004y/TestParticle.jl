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
using Meshes: coords, spacing, paramdim, CartesianGrid, RectilinearGrid, StructuredGrid
using ForwardDiff
using ChunkSplitters
using PrecompileTools: @setup_workload, @compile_workload
using Elliptic
import Tensors
import Base: +, *, /, setindex!, getindex
import LinearAlgebra: ×

export prepare, prepare_gc, get_gc, get_gc_func
export trace!, trace_relativistic!, trace_normalized!, trace_relativistic_normalized!,
       trace, trace_relativistic, trace_relativistic_normalized, trace_gc!, trace_gc_1st!,
       trace_gc_drifts!, trace_gc_flr!, trace_gc_exb!, trace_fieldline!, trace_fieldline
export Proton, Electron, Ion
export Maxwellian, BiMaxwellian
export Kappa, BiKappa
export AdaptiveBoris
export CurrentLoop, getB_loop
export get_gyrofrequency,
       get_gyroperiod, get_gyroradius, get_velocity, get_energy,
       energy2velocity,
       getB_zpinch, getB_bottle, getB_mirror, getB_tokamak_coil
export orbit, monitor
export TraceProblem, CartesianGrid, RectilinearGrid, StructuredGrid

"""
Type for the particles: `Proton`, `Electron`.
"""
struct Species{M, Q}
   m::M
   q::Q
end

Ion(m, q = 1) = Species(m * mᵢ, q * qᵢ)
Ion(; m = 1, q = 1) = Species(m * mᵢ, q * qᵢ)

include("types.jl")
include("fields/Fields.jl")
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

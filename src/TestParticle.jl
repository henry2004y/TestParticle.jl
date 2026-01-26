module TestParticle

using LinearAlgebra: norm, ×, ⋅, diag, normalize
using Statistics: mean
using Interpolations: interpolate, interpolate!, extrapolate, scale, BSpline, Linear,
    Quadratic, Cubic,
    Line, OnCell, Periodic, Flat, Gridded
using SciMLBase: AbstractODEProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
    BasicEnsembleAlgorithm, EnsembleThreads, EnsembleSerial,
    DEFAULT_SPECIALIZATION, ODEFunction, ODEProblem,
    LinearInterpolation, build_solution, ODESolution
using StaticArrays: SVector, @SMatrix, MVector, SA, StaticArray, SMatrix
using Meshes: coords, spacing, paramdim, CartesianGrid, RectilinearGrid, StructuredGrid
import ForwardDiff
using ChunkSplitters: index_chunks
using PrecompileTools: @setup_workload, @compile_workload
using MuladdMacro: @muladd

import Tensors
import Base: +, -, *, /, setindex!, getindex
import LinearAlgebra: ×

export prepare, prepare_gc, get_gc, get_gc_func
export trace!, trace_relativistic!, trace_normalized!, trace_relativistic_normalized!,
    trace, trace_relativistic, trace_relativistic_normalized,
    trace_gc!,
    trace_gc_drifts!, trace_gc_flr!, trace_gc_exb!, trace_fieldline!, trace_fieldline,
    get_gc_velocity, full_to_gc, gc_to_full
export Proton, Electron, Ion
export Maxwellian, BiMaxwellian
export Kappa, BiKappa
export AdaptiveBoris, AdaptiveHybrid
export CurrentLoop, getB_loop
export get_gyrofrequency,
    get_gyroperiod, get_gyroradius, get_velocity, get_energy,
    energy2velocity,
    get_curvature_radius, get_adiabaticity,
    sample_unit_sphere, get_number_density_flux,
    getB_zpinch, getB_bottle, getB_mirror, getB_tokamak_coil
export orbit, monitor
export TraceProblem, TraceGCProblem, TraceHybridProblem, CartesianGrid, RectilinearGrid, StructuredGrid

include("types.jl")
include("fields/Fields.jl")
include("utility/utility.jl")
include("utility/interpolation.jl")
include("sampler.jl")
include("prepare.jl")
include("gc.jl")
include("gc_solver.jl") # Added gc_solver.jl
include("equations.jl")
include("boris.jl")
include("adaptive_boris.jl")
include("hybrid.jl")
include("fieldline.jl")

function orbit end

function monitor end

include("precompile.jl")

end

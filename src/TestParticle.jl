module TestParticle

using LinearAlgebra: norm, ×, ⋅, diag, normalize
using Interpolations: interpolate, interpolate!, extrapolate, scale, BSpline, Linear,
    Quadratic, Cubic,
    Line, OnCell, Periodic, Flat, Gridded
using SciMLBase: AbstractODEProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
    BasicEnsembleAlgorithm, EnsembleThreads, EnsembleSerial, EnsembleDistributed,
    DEFAULT_SPECIALIZATION, ODEFunction, ODEProblem, remake,
    LinearInterpolation, build_solution, ODESolution
using Distributed: pmap
using StaticArrays: SVector, @SMatrix, MVector, SA, StaticArray, SMatrix
using Meshes: coords, spacing, paramdim, CartesianGrid, RectilinearGrid, StructuredGrid
import ForwardDiff
import DiffResults
using ChunkSplitters: index_chunks
using PrecompileTools: @setup_workload, @compile_workload
using MuladdMacro: @muladd
using KernelAbstractions: @kernel, @index, @Const, synchronize, Backend, CPU

import KernelAbstractions as KA
import Adapt
import Tensors
import Base: +, -, *, /, setindex!, getindex
import LinearAlgebra: ×

export prepare, prepare_gc, get_gc, get_gc_func
export trace!, trace_relativistic!, trace_normalized!, trace_relativistic_normalized!,
    trace, trace_relativistic, trace_normalized, trace_relativistic_normalized,
    trace_gc!,
    trace_gc_drifts!, trace_gc_flr!, trace_gc_exb!, trace_fieldline!, trace_fieldline,
    get_gc_velocity, full_to_gc, gc_to_full
export Proton, Electron, Ion
export Maxwellian, BiMaxwellian, Kappa, BiKappa
export AdaptiveBoris, AdaptiveHybrid
export get_gyrofrequency,
    get_gyroperiod, get_gyroradius, get_velocity, get_energy, get_mean_magnitude,
    energy2velocity, get_curvature_radius, get_adiabaticity,
    sample_unit_sphere, get_number_density_flux
export orbit, monitor
export get_fields, get_work
export LazyTimeInterpolator
export TraceProblem, TraceGCProblem, TraceHybridProblem, CartesianGrid, RectilinearGrid, StructuredGrid
export EnsembleSerial, EnsembleThreads, EnsembleDistributed, remake

include("types.jl")
include("utility/utility.jl")
include("utility/interpolation.jl")
include("sampler.jl")
include("prepare.jl")
include("gc.jl")
include("gc_solver.jl")
include("equations.jl")
include("boris.jl")
include("boris_kernel.jl")
include("adaptive_boris.jl")
include("hybrid.jl")
include("fieldline.jl")

function orbit end

function monitor end

include("precompile.jl")

end

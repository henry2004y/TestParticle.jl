module TestParticle

using LinearAlgebra: norm, ×, ⋅, diag, normalize
using FastInterpolations: constant_interp, linear_interp, cardinal_interp,
    interp, Extrap, PeriodicBC, ZeroCurvBC, gradient, deriv1,
    OnTheFly, PreCompute, FillExtrap, ClampExtrap, WrapExtrap, AbstractExtrap
using SciMLBase: AbstractODEProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
    BasicEnsembleAlgorithm,
    EnsembleThreads, EnsembleSerial, EnsembleDistributed, EnsembleSplitThreads,
    DEFAULT_SPECIALIZATION, ODEFunction, ODEProblem, remake,
    LinearInterpolation, build_solution, ODESolution, EnsembleSolution,
    DiscreteCallback, terminate!
using Distributed: pmap, nworkers
using StaticArrays: SVector, MVector, SA, StaticArray
using Meshes: coords, spacing, paramdim, CartesianGrid, RectilinearGrid, StructuredGrid,
    Plane, Disk, Point, normal, Sphere, area, Vec
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

export prepare, prepare_gc, get_gc, get_gc_func, ZeroField
export trace!, trace_relativistic!, trace_normalized!, trace_relativistic_normalized!,
    trace, trace_relativistic, trace_normalized, trace_relativistic_normalized,
    get_dx!, get_dv!,
    trace_gc!,
    trace_gc_drifts!, trace_gc_flr!, trace_gc_exb!,
    trace_fieldline!, trace_fieldline, TraceFieldlineProblem,
    get_gc_velocity, full_to_gc, gc_to_full
export Proton, Electron, Ion
export Maxwellian, BiMaxwellian, Kappa, BiKappa
export AdaptiveHybrid, AdaptiveBoris, Boris, MultistepBoris, MultistepBoris2, MultistepBoris4, MultistepBoris6
export get_gyrofrequency,
    get_gyroperiod, get_gyroradius, get_velocity, get_energy, get_mean_magnitude,
    energy2velocity, get_curvature_radius, get_adiabaticity,
    sample_unit_sphere, generate_sphere, sample_maxwellian,
    get_particle_flux, get_particle_fluxes, get_particle_crossings, get_first_crossing,
    sph2cart, cart2sph, sph2cartvec, cart2sphvec
export orbit, monitor
export get_fields, get_work
export LazyTimeInterpolator, build_interpolator
export TraceProblem, TraceGCProblem, TraceHybridProblem,
    CartesianGrid, RectilinearGrid, StructuredGrid
export EnsembleSerial, EnsembleThreads, EnsembleDistributed, EnsembleSplitThreads, remake,
    FillExtrap, ClampExtrap, WrapExtrap, OnTheFly, PreCompute,
    DiscreteCallback, TerminateOutside

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
include("hybrid.jl")
include("fieldline.jl")

function orbit end

function monitor end

include("precompile.jl")

end

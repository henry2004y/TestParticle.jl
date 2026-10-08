module TestParticle

using LinearAlgebra: norm, ×, ⋅, diag, normalize
using FastInterpolations: constant_interp, linear_interp, cardinal_interp,
    interp, Extrap, PeriodicBC, ZeroCurvBC, gradient, deriv1,
    OnTheFly, PreCompute, FillExtrap, ClampExtrap, WrapExtrap, AbstractExtrap,
    NoBC, NoExtrap
using SciMLBase: AbstractODEProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
    BasicEnsembleAlgorithm, EnsembleProblem,
    EnsembleThreads, EnsembleSerial, EnsembleDistributed, EnsembleSplitThreads,
    DEFAULT_SPECIALIZATION, ODEFunction, ODEProblem, remake,
    LinearInterpolation, build_solution, ODESolution, EnsembleSolution,
    DiscreteCallback, terminate!, EnsembleContext, isadaptive
import BorisPushers: Boris, AdaptiveBoris, MultistepBoris, AdaptiveMultistepBoris,
    MultistepBoris2, MultistepBoris4, MultistepBoris6, AbstractBoris,
    adapt_field_to_gpu, EnsembleKernel,
    get_q2m, get_EField, get_BField,
    update_velocity_boris, update_velocity_multistep, update_velocity,
    advance_boris, update_velocity_node, update_velocity_half, update_velocity_resync
import SciMLBase
import SciMLBase: solve
using Random: default_rng, AbstractRNG, Xoshiro
using StaticArrays: SVector, MVector, SA, StaticArray
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
    get_gc_velocity, full_to_gc, gc_to_full,
    get_B_parameters, get_E_parameters
export Proton, Electron, Ion
export Maxwellian, BiMaxwellian, Kappa, BiKappa
export AdaptiveBoris, AdaptiveMultistepBoris, AdaptiveHybrid,
    Boris, MultistepBoris,
    MultistepBoris2, MultistepBoris4, MultistepBoris6
export get_gyrofrequency,
    get_gyroperiod, get_gyroradius, get_velocity, get_energy, get_mean_magnitude,
    energy2velocity, get_curvature_radius, get_adiabaticity, adiabaticity_components,
    sample_unit_sphere, generate_sphere, sample_maxwellian,
    get_particle_flux, get_particle_fluxes, get_particle_crossings, get_first_crossing,
    sph2cart, cart2sph, sph2cartvec, cart2sphvec
export sample_velocity_ball, bin_centers, bin_velocity_space, project_vdf,
    analytic_projection, relative_l2, velocity_moments,
    vdf_grid_problem, vdf_backward, refine_vdf_window,
    vdf_backward_trace, vdf_forward_trace, embed_vdf
export orbit, monitor
export get_fields, get_work
export LazyTimeInterpolator, build_interpolator
export GPUGrid3D, GPUGrid2D, GPUGrid1D
export GPUSphericalGrid, GPUUniformAxis, GPUNonUniformAxis
export TraceProblem, TraceGCProblem, TraceHybridProblem, solve
export EnsembleProblem, EnsembleSerial, EnsembleThreads, EnsembleDistributed,
    EnsembleSplitThreads, EnsembleKernel, remake, FillExtrap, ClampExtrap, WrapExtrap,
    PeriodicBC, ZeroCurvBC, OnTheFly, PreCompute, DiscreteCallback, TerminateOutside

public ReturnCode, Field, Species, SpeciesDict, ZeroVector,
    CartesianGrid, RectilinearGrid, StructuredGrid,
    qₑ, mₑ, qᵢ, mᵢ, c, μ₀, ϵ₀, kB, Rₑ, BMoment_Earth, eV,
    get_q2m, get_EField, get_BField, jacobian, derivative_t,
    get_thermal_speed, get_cell_centers, get_magnetic_properties, get_work_rates_gc,
    makegrid

include("types.jl")
include("utility/utility.jl")
include("utility/interpolation.jl")
include("utility/gpu_grid.jl")
include("sampler.jl")
include("prepare.jl")
include("saveat.jl")
include("gc/gc.jl")
include("gc/gc_solver.jl")
include("gc/rk4_gc.jl")
include("gc/rk45_gc.jl")
include("equations.jl")
include("problem.jl")
include("boris/boris.jl")
include("boris/boris_solve.jl")
include("hybrid.jl")
include("fieldline.jl")

function orbit end

function monitor end

include("precompile.jl")

end

module TestParticle

using LinearAlgebra: norm, ×, ⋅, diag, normalize
using FastInterpolations: constant_interp, linear_interp, cardinal_interp,
    interp, Extrap, PeriodicBC, ZeroCurvBC, gradient, deriv1,
    OnTheFly, PreCompute, FillExtrap, ClampExtrap, WrapExtrap, AbstractExtrap,
    NoBC, NoExtrap, StorePolicy
using SciMLBase: AbstractODEProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
    BasicEnsembleAlgorithm, EnsembleProblem,
    EnsembleThreads, EnsembleSerial, EnsembleDistributed, EnsembleSplitThreads,
    DEFAULT_SPECIALIZATION, ODEFunction, ODEProblem, remake,
    LinearInterpolation, build_solution, ODESolution, EnsembleSolution,
    DiscreteCallback, terminate!, EnsembleContext, isadaptive
import BorisPushers: Boris, AdaptiveBoris, MultistepBoris, AdaptiveMultistepBoris,
    MultistepBoris2, MultistepBoris4, MultistepBoris6, AbstractBoris,
    adapt_field_to_gpu, EnsembleKernel,
    update_velocity_boris, update_velocity_multistep, update_velocity,
    advance_boris, update_velocity_node, update_velocity_half, update_velocity_resync
import SciMLBase
import SciMLBase: solve
using Random: default_rng, AbstractRNG, Xoshiro
using StaticArrays: SVector, MVector, SA, StaticArray, SMatrix, StaticVector
import ForwardDiff
import DiffResults
using ChunkSplitters: index_chunks
using PrecompileTools: @setup_workload, @compile_workload
using MuladdMacro: @muladd
using KernelAbstractions: @kernel, @index, @Const, synchronize, Backend, CPU

import KernelAbstractions as KA
import Adapt
import Base: +, -, *, /, setindex!, getindex
import LinearAlgebra: ×

export prepare, prepare_gc, prepare_canonical, ZeroField, PotentialField
export Proton, Electron, Ion
export TraceProblem, TraceCanonicalProblem, TraceGCProblem, TraceHybridProblem, solve
export trace!, trace_relativistic!, trace_normalized!, trace_relativistic_normalized!,
    trace, trace_relativistic, trace_normalized, trace_relativistic_normalized,
    trace_canonical!, trace_canonical,
    trace_canonical_relativistic!, trace_canonical_relativistic,
    trace_gc!, trace_gc_drifts!, trace_gc_flr!, trace_gc_exb!
export AdaptiveBoris, AdaptiveMultistepBoris, AdaptiveHybrid,
    Boris, MultistepBoris,
    MultistepBoris2, MultistepBoris4, MultistepBoris6
export vdf_backward_trace, vdf_forward_trace
export orbit, monitor
export get_fields, get_work
export build_interpolator
export EnsembleProblem, EnsembleSerial, EnsembleThreads, EnsembleDistributed,
    EnsembleSplitThreads, EnsembleKernel, remake, DiscreteCallback, TerminateOutside

public ReturnCode, Field, Species, SpeciesDict, ZeroVector,
    CartesianGrid, RectilinearGrid, StructuredGrid,
    qₑ, mₑ, qᵢ, mᵢ, c, μ₀, ϵ₀, kB, Rₑ, BMoment_Earth, eV,
    get_gc, get_gc_func, get_gc_velocity,
    full_to_gc, gc_to_full, velocity_to_canonical, canonical_to_velocity,
    canonical_hamiltonian,
    get_gyrofrequency, get_gyroperiod, get_gyroradius, get_velocity, get_energy,
    get_mean_magnitude, energy2velocity, get_curvature_radius, get_adiabaticity,
    adiabaticity_components,
    get_particle_flux, get_particle_fluxes, get_particle_crossings, get_first_crossing,
    Maxwellian, BiMaxwellian, Kappa, BiKappa,
    sample_unit_sphere, sample_maxwellian, sample_velocity_ball, generate_sphere,
    jacobian, derivative_t,
    get_thermal_speed, get_cell_centers, get_magnetic_properties,
    get_work_rates, get_work_rates_gc, get_work_channels,
    makegrid, get_dx!, get_dv!, get_B_parameters, get_E_parameters,
    LazyTimeInterpolator,
    GPUGrid1D, GPUGrid2D, GPUGrid3D, GPUSphericalGrid, GPUUniformAxis,
    GPUNonUniformAxis,
    sph2cart, cart2sph, sph2cartvec, cart2sphvec,
    bin_centers, bin_velocity_space, project_vdf, analytic_projection, relative_l2,
    velocity_moments, vdf_grid_problem, vdf_backward, refine_vdf_window,
    embed_vdf,
    FillExtrap, ClampExtrap, WrapExtrap, PeriodicBC, OnTheFly, PreCompute, StorePolicy

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

function orbit end

function monitor end

include("precompile.jl")

end

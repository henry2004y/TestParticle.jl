module TestParticle

using LinearAlgebra: norm, ×, ⋅, normalize
using FastInterpolations: constant_interp, linear_interp, cardinal_interp,
    Extrap, PeriodicBC, gradient, deriv1,
    OnTheFly, PreCompute, FillExtrap, ClampExtrap, WrapExtrap, AbstractExtrap,
    NoBC, NoExtrap, StorePolicy
using SciMLBase: AbstractODEProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
    BasicEnsembleAlgorithm, EnsembleProblem,
    EnsembleThreads, EnsembleSerial, EnsembleDistributed, EnsembleSplitThreads,
    DEFAULT_SPECIALIZATION, ODEFunction, remake,
    LinearInterpolation, build_solution, EnsembleSolution,
    DiscreteCallback, terminate!, EnsembleContext
import SciMLBase
import SciMLBase: solve
using Random: default_rng, AbstractRNG, Xoshiro
using StaticArrays: SVector, MVector, SA, StaticArray, SMatrix
import ForwardDiff
import DiffResults
using ChunkSplitters: index_chunks
using PrecompileTools: @setup_workload, @compile_workload
using MuladdMacro: @muladd

import Adapt
import Base: +, -, *, /
import LinearAlgebra: ×

export prepare, prepare_gc, prepare_canonical, ZeroField, PotentialField
export Proton, Electron, Ion
export TraceProblem, TraceCanonicalProblem, TraceGCProblem, TraceHybridProblem, solve
export trace!, trace_relativistic!, trace_normalized!, trace_relativistic_normalized!,
    trace, trace_relativistic, trace_normalized, trace_relativistic_normalized,
    trace_canonical!, trace_canonical,
    trace_canonical_relativistic!, trace_canonical_relativistic,
    trace_gc!, trace_gc_drifts!, trace_gc_flr!, trace_gc_exb!
export AdaptiveHybrid
export vdf_backward_trace, vdf_forward_trace
export orbit, monitor
export get_fields, get_work
export build_interpolator
export EnsembleProblem, EnsembleSerial, EnsembleThreads, EnsembleDistributed,
    EnsembleSplitThreads, remake, DiscreteCallback, TerminateOutside

public ReturnCode, Field, Species, SpeciesDict, ZeroVector,
    CartesianGrid, RectilinearGrid, StructuredGrid,
    qₑ, mₑ, qᵢ, mᵢ, c, μ₀, ϵ₀, kB, Rₑ, BMoment_Earth, eV,
    get_gc, get_gc_velocity,
    full_to_gc, gc_to_full, velocity_to_canonical, canonical_to_velocity,
    canonical_hamiltonian,
    get_gyrofrequency, get_gyroperiod, get_gyroradius, get_energy,
    energy2velocity, get_curvature_radius, get_adiabaticity,
    adiabaticity_components,
    get_particle_flux, get_particle_crossings, get_first_crossing,
    Maxwellian, BiMaxwellian, Kappa, BiKappa,
    sample_unit_sphere, sample_maxwellian, sample_velocity_ball, generate_sphere,
    jacobian, derivative_t,
    get_magnetic_properties, get_work_rates, get_work_channels,
    get_dx!, get_dv!,
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
include("utility/field_geometry.jl")
include("utility/work.jl")
include("utility/interpolation.jl")
include("utility/gpu_grid.jl")
include("sampler.jl")
include("prepare.jl")
include("saveat.jl")
include("gc/gc.jl")
include("gc/gc_solver.jl")
include("gc/rk4_gc.jl")
include("gc/rk45_gc.jl")
include("gc/equations.jl")
include("equations.jl")
include("canonical.jl")
include("problem.jl")
include("hybrid.jl")

function orbit end

function monitor end

include("precompile.jl")

end

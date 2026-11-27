module TestParticle

using LinearAlgebra: norm, ×, ⋅, normalize
using FastInterpolations: constant_interp, linear_interp, cardinal_interp,
    Extrap, PeriodicBC, gradient, deriv1,
    OnTheFly, PreCompute, FillExtrap, ClampExtrap, WrapExtrap, AbstractExtrap,
    NoBC, NoExtrap, StorePolicy
using SciMLBase: AbstractODEProblem, AbstractODEFunction, AbstractODESolution, ReturnCode,
    BasicEnsembleAlgorithm, EnsembleProblem,
    EnsembleThreads, EnsembleSerial, EnsembleDistributed, EnsembleSplitThreads,
    DEFAULT_SPECIALIZATION, FullSpecialize, ODEFunction, ODEProblem, remake,
    LinearInterpolation, build_solution,
    DiscreteCallback, terminate!, EnsembleContext, isadaptive
import SciMLBase
import SciMLBase: solve
using Random: default_rng, AbstractRNG, Xoshiro
using StaticArrays: SVector, SA, StaticArray, SMatrix
import ForwardDiff
import DiffResults
using ChunkSplitters: index_chunks
using FunctionWrappers: FunctionWrapper
using PrecompileTools: @setup_workload, @compile_workload
using MuladdMacro: @muladd
using PlasmaPopulations: Particle, Proton, Electron, Positron, Population, pdf,
    CartesianGridVDF,
    gyrofrequency, C_SI, E_SI, MU0_SI, EPS0_SI, KB_SI

import Adapt
import Base: +, -, *, /
import LinearAlgebra: ×

export prepare, prepare_canonical, prepare_erased, PotentialField
export Particle, Proton, Electron, Positron
export TraceProblem, TraceCanonicalProblem, TraceGCProblem, solve
export FullOrbit, RelativisticOrbit, Canonical, RelativisticCanonical,
    GuidingCenter, RelativisticGuidingCenter
export trace!, trace_relativistic!, trace, trace_relativistic,
    trace_canonical!, trace_canonical,
    trace_canonical_relativistic!, trace_canonical_relativistic,
    trace_gc!, trace_gc_relativistic!, trace_gc_relativistic,
    trace_gc_drifts!, trace_gc_flr!, trace_gc_exb!
export AdaptiveHybrid
export vdf_backward_trace, vdf_forward_trace
export orbit, monitor, save_trajectory
export get_fields, get_work_rates
export build_interpolator
export EnsembleProblem, EnsembleSerial, EnsembleThreads, EnsembleDistributed,
    EnsembleSplitThreads, remake, DiscreteCallback, TerminateOutside

public ReturnCode, Field, ErasedField, ErasedFunction,
    ZeroField, ZeroVector,
    CartesianGrid, RectilinearGrid, StructuredGrid,
    c, Rₑ, BMoment_Earth, eV,
    get_gc, get_gc_velocity,
    full_to_gc, gc_to_full, velocity_to_canonical, canonical_to_velocity,
    canonical_hamiltonian,
    get_gyroperiod, get_gyroradius, get_energy,
    energy2velocity, get_curvature_radius, get_adiabaticity,
    adiabaticity_components,
    get_particle_flux, get_particle_crossings, get_first_crossing,
    sample_unit_sphere, sample_velocity_ball, generate_sphere,
    jacobian, derivative_t,
    get_magnetic_properties, get_work_channels,
    get_dx!, get_dv!,
    LazyTimeInterpolator,
    GPUGrid, GPUUniformAxis, GPUNonUniformAxis,
    sph2cart, cart2sph, sph2cartvec, cart2sphvec,
    bin_velocity_space, project_vdf,
    FillExtrap, ClampExtrap, WrapExtrap, PeriodicBC, OnTheFly, PreCompute, StorePolicy

include("types.jl")
include("utility/utility.jl")
include("utility/field_geometry.jl")
include("utility/work.jl")
include("utility/interpolation.jl")
include("utility/gpu_grid.jl")
include("prepare.jl")
include("saveat.jl")
include("gc/gc.jl")
include("gc/equations.jl")
include("equations.jl")
include("canonical.jl")
include("problem.jl")
include("hybrid.jl")
include("deprecated.jl")

function orbit end

function monitor end

function save_trajectory end

function save_trajectory(args...; kwargs...)
    throw(
        ArgumentError(
            "save_trajectory requires JLD2.jl. Please run `using JLD2` to enable trajectory saving."
        )
    )
end

include("precompile.jl")

end

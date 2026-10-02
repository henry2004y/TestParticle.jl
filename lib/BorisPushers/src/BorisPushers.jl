module BorisPushers

using Reexport
@reexport using SciMLBase
import SciMLBase: AbstractODEProblem, EnsembleProblem, EnsembleSolution,
    EnsembleSerial, EnsembleThreads, BasicEnsembleAlgorithm,
    ReturnCode, build_solution, LinearInterpolation, solve
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveAlgorithm,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    AbstractController, AbstractControllerCache,
    alg_order, alg_cache, isfsal, initialize!, perform_step!,
    accept_step_controller, default_controller, setup_controller_cache,
    step_accept_controller!, step_reject_controller!, stepsize_controller!,
    _ode_addsteps!, default_linear_interpolation, get_fsalfirstlast,
    postamble!, _postamble!, ODEIntegrator, ODE_DEFAULT_ISOUTOFDOMAIN
using StaticArrays
using MuladdMacro
using LinearAlgebra
using KernelAbstractions: @kernel, @index, @Const, synchronize, Backend, CPU
import KernelAbstractions as KA
import Adapt
using ChunkSplitters: index_chunks

include("algorithms.jl")
include("alg_utils.jl")
include("parameters.jl")
include("boris_controller.jl")
include("boris_caches.jl")
include("dense_output.jl")
include("boris_step.jl")
include("boris_perform_step.jl")
include("boris_solve.jl")
include("saveat.jl")
include("boris_kernel.jl")

export Boris, AdaptiveBoris
export MultistepBoris, MultistepBoris2, MultistepBoris4, MultistepBoris6
export AdaptiveMultistepBoris
export AbstractBoris
export GPUBorisAlgorithm
export adapt_field_to_gpu
export SavingPlan, use_saveat
export get_q2m, get_EField, get_BField
export update_velocity_boris, update_velocity_multistep, update_velocity
export advance_boris, update_velocity_node, update_velocity_half, update_velocity_resync

end

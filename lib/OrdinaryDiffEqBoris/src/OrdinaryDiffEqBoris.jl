module OrdinaryDiffEqBoris

using Reexport
@reexport using SciMLBase
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveAlgorithm,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    AbstractController, AbstractControllerCache,
    alg_order, alg_cache, isfsal, initialize!, perform_step!,
    accept_step_controller, default_controller, setup_controller_cache,
    step_accept_controller!, step_reject_controller!, stepsize_controller!,
    _ode_addsteps!, default_linear_interpolation, get_fsalfirstlast
using StaticArrays
using MuladdMacro
using LinearAlgebra

include("algorithms.jl")
include("alg_utils.jl")
include("parameters.jl")
include("boris_controller.jl")
include("boris_caches.jl")
include("dense_output.jl")
include("boris_step.jl")
include("boris_perform_step.jl")

export Boris, AdaptiveBoris
export MultistepBoris, MultistepBoris2, MultistepBoris4, MultistepBoris6
export AdaptiveMultistepBoris
export AbstractBoris
export get_q2m, get_EField, get_BField
export update_velocity_boris, update_velocity_multistep, update_velocity
export boris_advance, boris_node_velocity, boris_half_velocity, boris_resync_velocity

end

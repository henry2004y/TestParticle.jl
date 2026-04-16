module OrdinaryDiffEqBoris

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveAlgorithm,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    alg_order, alg_cache, isfsal, initialize!, perform_step!


using SciMLBase
using RecursiveArrayTools
using StaticArrays
using MuladdMacro
using LinearAlgebra

include("algorithms.jl")
include("alg_utils.jl")
include("boris_caches.jl")
include("boris_perform_step.jl")

export Boris
export MultistepBoris
export AdaptiveBoris
export boris_velocity_update # Exposing strictly for debugging if needed
export update_velocity_multistep

end

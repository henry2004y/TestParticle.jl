module OrdinaryDiffEqBoris

using Reexport
@reexport using SciMLBase
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveAlgorithm,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    alg_order, alg_cache, isfsal, initialize!, perform_step!
using RecursiveArrayTools
using StaticArrays
using MuladdMacro
using LinearAlgebra

include("algorithms.jl")
include("alg_utils.jl")
include("parameters.jl")
include("boris_caches.jl")
include("boris_perform_step.jl")

export Boris
export MultistepBoris2, MultistepBoris4, MultistepBoris6
export get_q2m, get_EField, get_BField

end

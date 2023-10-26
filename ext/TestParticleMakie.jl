module TestParticleMakie

import TestParticle: orbit, monitor
using Makie
using SciMLBase: AbstractODESolution

include("typerecipes.jl")
include("interactive.jl")

# export orbit, monitor

end # module TestParticleMakie

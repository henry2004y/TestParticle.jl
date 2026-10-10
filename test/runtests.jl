using Test

# Shared fixtures. Each test file guards its own include of this, so any single
# test file can still be run on its own from the REPL.
include("test_common.jl")

@testset "TestParticle.jl" begin
    include("test_sampling.jl")
    include("test_numerical_field.jl")
    include("test_mixed_fields.jl")
    include("test_time_dependent_fields.jl")
    include("test_field_exceptions.jl")
    include("test_relativistic.jl")
    include("test_normalized_fields.jl")
    include("test_zerovector.jl")
    include("test_derivatives.jl")
end

include("test_boris.jl")
include("test_utility.jl")
include("test_phasespace.jl")
include("test_fieldline.jl")
include("test_gc.jl")

if "makie" in ARGS
    include("test_Makie.jl")
end

include("test_distributions.jl")
include("test_hybrid.jl")
include("test_adiabaticity.jl")
include("test_boris_kernel.jl")
include("test_boris_gpu.jl")
include("test_spherical_gpu.jl")
include("test_raw_output.jl")
include("test_boundary.jl")
include("test_adaptive_boris.jl")
include("test_adaptive_multistep_boris.jl")
include("test_reproducibility.jl")
include("test_seed_reproducibility.jl")
include("test_symplectic.jl")
include("test_traceproblem.jl")

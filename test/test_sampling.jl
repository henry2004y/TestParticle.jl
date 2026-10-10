module test_sampling

using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using Random
using LinearAlgebra: norm
# Loading the backend package is what attaches the sampling methods; both it and
# TestParticle export `Maxwellian`, so the binding has to be named explicitly.
using VelocityDistributionFunctions
import TestParticle: Maxwellian, BiMaxwellian

@testset "sampling" begin
    u0 = SA[0.0, 0.0, 0.0]
    p = 1.0e-9 # [Pa]
    n = 1.0e6 # [/m³]
    vdf = Maxwellian(u0, p, n)
    Random.seed!(1234)
    v = rand(vdf)
    @test sum(v) ≈ 238065.6009276599
    @test occursin("Maxwellian", repr(vdf))
    B = [1.0, 0.0, 0.0] # will be normalized internally
    vdf = BiMaxwellian(B, u0, p, p, n)
    Random.seed!(1234)
    v = rand(vdf)
    @test sum(v) ≈ -794362.2464141053
    @test occursin("BiMaxwellian", repr(vdf))

    # Maxwellian velocity sampling
    Tn = 3000.0
    m = 16 * TP.mᵢ
    v_sample = sample_maxwellian(Tn, m)
    @test length(v_sample) == 3

    v_off = sample_maxwellian(Tn, "O+"; offset = 2 * TP.eV)
    @test norm(v_off) > sqrt(2 * 2 * TP.qᵢ / (16 * TP.mᵢ))
end

end # module test_sampling

module test_zerovector

using Test
using TestParticle

@testset "ZeroVector" begin
    z = TestParticle.ZeroVector()
    @test length(z) == 3
    @test collect(z) == [0, 0, 0]
end

end # module test_zerovector

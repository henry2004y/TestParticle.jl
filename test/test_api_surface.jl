module test_api_surface

using Test
import TestParticle as TP

# Julia accepts `export`/`public` of names that are never defined; users only hit the
# `UndefVarError` on use, e.g. after a rename that misses the declaration.
@testset "Public names are defined" begin
    @test all(name -> isdefined(TP, name), names(TP))
end

end

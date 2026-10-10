module test_field_exceptions

using Test
using TestParticle: Field
using StaticArrays

"Time-dependent field that `Field` cannot decide on from its signature alone."
E_ambiguous(r, t) = SA[5.0e-11 * sin(2π * t), 0, 0]

"Field signature with three arguments, which `Field` does not support."
F_unsupported(r, v, t) = SA[r, v, t]

@testset "field exceptions" begin
    E = Field(E_ambiguous)
    @test_throws ArgumentError E([0, 0, 0])
    # An unsupported function form is accepted, but flagged as not time-dependent
    @test typeof(Field(F_unsupported)).parameters[1] == false
end

end # module test_field_exceptions

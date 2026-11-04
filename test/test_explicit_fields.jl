module test_explicit_fields
using Test, TestParticle, StaticArrays
using TestParticle: Field, PeriodicBC, FillExtrap, ClampExtrap, jacobian
struct StaticCallable end
(::StaticCallable)(x) = SVector(x[1]^2, x[2], 0.0)
struct BothCallable end
(::BothCallable)(x) = SVector(x[1], 0.0, 0.0)
(::BothCallable)(x, t) = SVector(x[1] + t, 0.0, 0.0)
@testset "Explicit fields" begin
    x = 0.0:1.0:3.0
    y = 0.0:1.0:2.0
    A = [SVector(a, b, 0.0) for a in x, b in y]
    f = build_interpolator(
        A, (x, y);
        bc = (PeriodicBC(; endpoint = :exclusive), FillExtrap(NaN))
    )
    @test f([3.5, 0.5, 0.0]) ≈ SVector(1.5, 0.5, 0.0)
    @test all(isnan, f([0.0, 3.0, 0.0]))
    @test_throws ArgumentError build_interpolator(A, (x, y); bc = (ClampExtrap(),))
    @test_throws ArgumentError build_interpolator(
        A, ([0.0, 1.0, 1.0, 3.0], y);
        bc = ClampExtrap()
    )
    @test_throws ArgumentError prepare(; B = A)
    r = SVector(2.0, 1.0, 0.0)
    p = prepare(; B = StaticCallable())
    @test jacobian(p.B, r, 0.0) ≈ [4 0 0; 0 1 0; 0 0 0]
    @test prepare_erased(; B = StaticCallable()).B(r, 0.0) == p.B(r, 0.0)
    @test prepare(; B = Field(BothCallable(); time_dependent = false)).B(r, 3.0)[1] == 2
    @test prepare(; B = Field(BothCallable(); time_dependent = true)).B(r, 3.0)[1] == 5
    for order in (0, 3)
        periodic = build_interpolator(
            A, (x, y); order,
            bc = (PeriodicBC(; endpoint = :exclusive), ClampExtrap())
        )
        @test periodic([4.2, 3.0, 0.0]) ≈ periodic([0.2, 2.0, 0.0])
    end
    r = range(1.0, 2.0; length = 4)
    theta = range(0.0, π; length = 4)
    phi = range(0.0, 2π; length = 9)[1:(end - 1)]
    data = fill(SVector(1.0, 0.0, 0.0), 4, 4, 8)
    policies = (ClampExtrap(), ClampExtrap(), PeriodicBC(; endpoint = :exclusive))
    spherical = build_interpolator(
        data, (r, theta, phi);
        bc = policies, gridtype = TestParticle.StructuredGrid
    )
    @test spherical([0.0, -1.5, 0.0]) ≈ SVector(0.0, -1.0, 0.0)
    @test_throws ArgumentError build_interpolator(
        data, (r, theta, phi ./ 2);
        bc = policies, gridtype = TestParticle.StructuredGrid
    )
    @test_throws MethodError prepare(; B = StaticCallable(), order = 3)
    # Simulation output need not match exactly across the longitude seam.
    seam = range(0.0, 2π; length = 9)
    mismatched = [SVector(1.0, 0.0, k / 100) for _ in r, _ in theta, k in eachindex(seam)]
    @test build_interpolator(
        TestParticle.StructuredGrid, mismatched, r, theta, seam
    ) isa Function
end
end

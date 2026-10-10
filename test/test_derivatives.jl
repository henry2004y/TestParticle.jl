module test_derivatives

using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using Meshes: CartesianGrid

"Linear ramp field, whose Jacobian is constant and known to ForwardDiff."
ramp_field(i, j, k = 0) = SA[Float64(i), Float64(j), Float64(k)]

@testset "Field Derivatives" begin
    @testset "3D Jacobian" begin
        x = range(0, 1, length = 5)
        y = range(0, 1, length = 5)
        z = range(0, 1, length = 5)
        A = [ramp_field(i, j, k) for i in 1:5, j in 1:5, k in 1:5]

        itp = build_interpolator(CartesianGrid, A, x, y, z, 1, FillExtrap(NaN))
        f = TP.Field(itp)

        pos = SA[0.5, 0.5, 0.5]
        @test TP.jacobian(f, pos, 0.0) ≈ TP.ForwardDiff.jacobian(r -> f(r, 0.0), pos)
    end

    @testset "2D Jacobian" begin
        x = range(0, 1, length = 5)
        y = range(0, 1, length = 5)
        A = [ramp_field(i, j) for i in 1:5, j in 1:5]

        itp = build_interpolator(CartesianGrid, A, x, y, 1, FillExtrap(NaN))
        f = TP.Field(itp)

        pos = SA[0.5, 0.5]
        # Bfunc is called as f(xu, t), where xu may be 2D or 3D; the 2D
        # interpolator only reads xu[1] and xu[2].
        @test TP.jacobian(f, pos, 0.0) ≈ TP.ForwardDiff.jacobian(r -> f(r, 0.0), pos)
    end

    @testset "1D Jacobian" begin
        x = range(0, 1, length = 5)
        A = [ramp_field(i, 0) for i in 1:5]

        itp = build_interpolator(CartesianGrid, A, x, 1, FillExtrap(NaN); dir = 1)
        f = TP.Field(itp)

        pos = SA[0.5]
        @test TP.jacobian(f, pos, 0.0) ≈ TP.ForwardDiff.jacobian(r -> f(r, 0.0), pos)
    end

    @testset "LazyTimeInterpolator derivatives" begin
        times = [0.0, 1.0, 2.0]
        loader(i) = pos -> SA[pos[1] * times[i], pos[2] * times[i], pos[3] * times[i]]
        itp = LazyTimeInterpolator(times, loader)
        f = TP.Field(itp)

        pos = SA[1.0, 2.0, 3.0]
        t = 0.5

        # Space Jacobian
        @test TP.jacobian(f, pos, t) ≈ TP.ForwardDiff.jacobian(r -> f(r, t), pos)
        # Time derivative
        @test TP.derivative_t(f, pos, t) ≈ TP.ForwardDiff.derivative(τ -> f(pos, τ), t)

        # Out-of-bounds Jacobian (clamp)
        @test TP.jacobian(f, pos, -1.0) == TP.jacobian(f, pos, 0.0)
        @test TP.jacobian(f, pos, 3.0) == TP.jacobian(f, pos, 2.0)

        # Out-of-bounds time derivative (zero)
        @test TP.derivative_t(f, pos, -1.0) == zero(f(pos, -1.0))
        @test TP.derivative_t(f, pos, 3.0) == zero(f(pos, 3.0))
    end
end

end # module test_derivatives

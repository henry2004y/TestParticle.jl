using Test
using TestParticle
import TestParticle as TP
using StaticArrays

@testset "Field Derivatives" begin
    let
        # 3D Jacobian Test
        x = range(0, 1, length = 5)
        y = range(0, 1, length = 5)
        z = range(0, 1, length = 5)
        A = [SA[Float64(i), Float64(j), Float64(k)] for i in 1:5, j in 1:5, k in 1:5]

        itp = TP.build_interpolator(CartesianGrid, A, x, y, z, 1, 1) # order=1, bc=1
        f = TP.Field(itp)

        pos = SA[0.5, 0.5, 0.5]
        J_analytical = TP.jacobian(f, pos, 0.0)
        J_forwarddiff = TP.ForwardDiff.jacobian(r -> f(r, 0.0), pos)

        @test J_analytical ≈ J_forwarddiff
    end

    let
        # 2D Jacobian Test
        x = range(0, 1, length = 5)
        y = range(0, 1, length = 5)
        A = [SA[Float64(i), Float64(j), 0.0] for i in 1:5, j in 1:5]

        itp = TP.build_interpolator(CartesianGrid, A, x, y, 1, 1)
        f = TP.Field(itp)

        pos = SA[0.5, 0.5]
        J_analytical = TP.jacobian(f, pos, 0.0)
        # Bfunc is called as f(xu, t), where xu might be 2D or 3D.
        # FieldInterpolator2D expects xu[1], xu[2].
        J_forwarddiff = TP.ForwardDiff.jacobian(r -> f(r, 0.0), pos)

        @test J_analytical ≈ J_forwarddiff
    end

    let
        # 1D Jacobian Test
        x = range(0, 1, length = 5)
        A = [SA[Float64(i), 0.0, 0.0] for i in 1:5]

        itp = TP.build_interpolator(CartesianGrid, A, x, 1, 1; dir = 1)
        f = TP.Field(itp)

        pos = SA[0.5]
        J_analytical = TP.jacobian(f, pos, 0.0)
        J_forwarddiff = TP.ForwardDiff.jacobian(r -> f(r, 0.0), pos)

        @test J_analytical ≈ J_forwarddiff
    end

    let
        # LazyTimeInterpolator Derivative Test
        times = [0.0, 1.0, 2.0]
        loader(i) = pos -> SA[pos[1] * times[i], pos[2] * times[i], pos[3] * times[i]]

        itp = LazyTimeInterpolator(times, loader)
        f = TP.Field(itp)

        pos = SA[1.0, 2.0, 3.0]
        t = 0.5

        # Space Jacobian
        J_analytical = TP.jacobian(f, pos, t)
        J_forwarddiff = TP.ForwardDiff.jacobian(r -> f(r, t), pos)
        @test J_analytical ≈ J_forwarddiff

        # Time Derivative
        dT_analytical = TP.derivative_t(f, pos, t)
        dT_forwarddiff = TP.ForwardDiff.derivative(τ -> f(pos, τ), t)
        @test dT_analytical ≈ dT_forwarddiff

        # Out-of-bounds Jacobian (clamp)
        @test TP.jacobian(f, pos, -1.0) == TP.jacobian(f, pos, 0.0)
        @test TP.jacobian(f, pos, 3.0) == TP.jacobian(f, pos, 2.0)

        # Out-of-bounds Time Derivative (zero)
        @test TP.derivative_t(f, pos, -1.0) == zero(f(pos, -1.0))
        @test TP.derivative_t(f, pos, 3.0) == zero(f(pos, 3.0))
    end
end

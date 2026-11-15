module test_derivatives

using Test
using TestParticle
using TestParticle: FillExtrap, LazyTimeInterpolator, PeriodicBC, WrapExtrap
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

        # Range grid and exact-node ForwardDiff derivatives
        let
            itp_r = LazyTimeInterpolator(0.0:1.0:2.0, loader)
            fr = TP.Field(itp_r)
            @test TP.jacobian(fr, pos, 0.0) ≈ TP.ForwardDiff.jacobian(r -> fr(r, 0.0), pos)
            @test TP.derivative_t(fr, pos, 0.0) ≈
                TP.ForwardDiff.derivative(τ -> fr(pos, τ), 0.0)
            @test TP.derivative_t(fr, pos, 1.0) ≈
                TP.ForwardDiff.derivative(τ -> fr(pos, τ), 1.0)
        end

        # FillExtrap derivatives
        let
            itp_nan = LazyTimeInterpolator(0.0:1.0:2.0, loader; bc = FillExtrap(NaN))
            fn = TP.Field(itp_nan)
            @test all(isnan, TP.jacobian(fn, pos, -1.0))
            @test all(isnan, TP.derivative_t(fn, pos, 3.0))

            itp_zero = LazyTimeInterpolator(0.0:1.0:2.0, loader; bc = FillExtrap(0.0))
            fz = TP.Field(itp_zero)
            @test TP.jacobian(fz, pos, -1.0) == zero(TP.jacobian(fz, pos, 0.5))
            @test TP.derivative_t(fz, pos, 3.0) == zero(pos)
        end

        # WrapExtrap and PeriodicBC derivatives
        let
            itp_wrap = LazyTimeInterpolator(0.0:1.0:2.0, loader; bc = WrapExtrap())
            fw = TP.Field(itp_wrap)
            for tq in (0.4, -0.4, 2.4)
                @test TP.jacobian(fw, pos, tq) ≈
                    TP.ForwardDiff.jacobian(r -> fw(r, tq), pos)
                @test TP.derivative_t(fw, pos, tq) ≈
                    TP.ForwardDiff.derivative(τ -> fw(pos, τ), tq)
            end

            itp_incl = LazyTimeInterpolator(0.0:1.0:2.0, loader; bc = PeriodicBC())
            fi = TP.Field(itp_incl)
            @test TP.derivative_t(fi, pos, -0.4) ≈
                TP.ForwardDiff.derivative(τ -> fi(pos, τ), -0.4)

            itp_excl = LazyTimeInterpolator(
                0.0:1.0:2.0, loader; bc = PeriodicBC(endpoint = :exclusive)
            )
            fe = TP.Field(itp_excl)
            for tq in (0.5, 2.5, -0.5, 3.5)
                @test TP.jacobian(fe, pos, tq) ≈
                    TP.ForwardDiff.jacobian(r -> fe(r, tq), pos)
                @test TP.derivative_t(fe, pos, tq) ≈
                    TP.ForwardDiff.derivative(τ -> fe(pos, τ), tq)
            end
        end

        # Shared multi-field component derivatives
        let
            multi = LazyTimeInterpolator(
                0.0:1.0:2.0,
                i -> (
                    r -> SA[r[1] * i, r[2], r[3]],
                    r -> SA[r[1] * times[i], r[2] * times[i], r[3] * times[i]],
                ),
            )
            _, B_itp = multi
            fB = TP.Field(B_itp)
            @test TP.jacobian(fB, pos, 0.5) ≈ TP.ForwardDiff.jacobian(r -> fB(r, 0.5), pos)
            @test TP.derivative_t(fB, pos, 0.5) ≈
                TP.ForwardDiff.derivative(τ -> fB(pos, τ), 0.5)
        end
    end
end

end # module test_derivatives

module test_interpolation
using Test
using TestParticle
using StaticArrays
import TestParticle as TP

@testset "Interpolation" begin
    @testset "Periodic Quadratic 1D" begin
        x = range(0.0, 1.0, length=11)
        f(xi) = sin(2π*xi)
        # Vector field for build_interpolator (3, nx) or SVector{3} array
        B = [SA[f(xi), 0.0, 0.0] for xi in x]
        
        # order=2, bc=2
        itp = TP.build_interpolator(B, x, 2, 2)
        
        # Test inside
        for xi in 0.1:0.2:0.9
            @test itp(SA[xi])[1] ≈ f(xi) atol=1e-2
        end
        
        # Test periodicity (bc=2 -> Extrap(:wrap))
        @test itp(SA[1.2])[1] ≈ f(1.2) atol=1e-2
        @test itp(SA[1.2])[1] ≈ f(0.2) atol=1e-2
        @test itp(SA[-0.8])[1] ≈ f(0.2) atol=1e-2
    end

    @testset "Periodic Quadratic 3D" begin
        x = range(0.0, 1.0, length=11)
        y = range(0.0, 1.0, length=11)
        z = range(0.0, 1.0, length=11)
        
        f(xi, yi, zi) = sin(2π*xi) + cos(2π*yi) + sin(2π*zi)
        
        # Create a 3D array of SVectors
        B = [SA[f(xi, yi, zi), 0.0, 0.0] for xi in x, yi in y, zi in z]
        
        # order=2, bc=2
        itp = TP.build_interpolator(B, x, y, z, 2, 2)
        
        # Test inside
        pt = SA[0.25, 0.5, 0.75]
        @test itp(pt)[1] ≈ f(pt[1], pt[2], pt[3]) atol=1e-2
        
        # Test periodicity
        pt_off = SA[1.25, -0.5, 2.75]
        pt_wrapped = SA[0.25, 0.5, 0.75]
        @test itp(pt_off)[1] ≈ itp(pt_wrapped)[1] atol=1e-8
    end

    @testset "Mixed Precision Error" begin
        x32 = range(0.0f0, 1.0f0, length=11)
        # Field in Float32 SVectors
        B32 = [SA[Float32(sin(2π*xi)), 0.0f0, 0.0f0] for xi in x32, yi in x32, zi in x32]
        # Grid in Float64
        x64 = range(0.0, 1.0, length=11)

        # Should work for order=1
        @test_nowarn TP.build_interpolator(B32, x64, x64, x64, 1, 1)

        # Should throw ArgumentError for order=2
        @test_throws ArgumentError TP.build_interpolator(B32, x64, x64, x64, 2, 1)
    end
end

end # module test_interpolation

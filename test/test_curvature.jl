using Test
using StaticArrays
using TestParticle
using LinearAlgebra

@testset "Curvature Radius" begin
   @testset "Zero Field" begin
      B_zero(x, t) = SA[0.0, 0.0, 0.0]
      x = SA[1.0, 0.0, 0.0]
      @test get_curvature_radius(x, 0.0, B_zero) == Inf
   end

   @testset "Uniform Field" begin
      B_uniform(x, t) = SA[0.0, 0.0, 1.0]
      x = SA[1.0, 0.0, 0.0]
      @test get_curvature_radius(x, 0.0, B_uniform) == Inf
   end

   @testset "Circular Field" begin
      # Field lines are circles in x-z plane: B = (-z, 0, x)
      # Radius of curvature at (R, 0, 0) is R.
      function B_circular(x, t)
         r = hypot(x[1], x[3])
         r == 0 && return SA[0.0, 0.0, 0.0]
         SA[-x[3], 0.0, x[1]]
      end

      R = 5.0
      @test get_curvature_radius(SA[R, 0.0, 0.0], 0.0, B_circular) ≈ R
      @test get_curvature_radius(SA[0.0, 0.0, R], 0.0, B_circular) ≈ R

      theta = π/4
      @test get_curvature_radius(SA[R*cos(theta), 0.0, R*sin(theta)], 0.0, B_circular) ≈ R
   end
end

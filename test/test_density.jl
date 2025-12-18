using TestParticle
using Meshes
using StaticArrays
using Test
using OrdinaryDiffEq

@testset "Number Density Calculation" begin
   # Mock TraceSolution
   struct MockSol
      t::StepRangeLen{
         Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
      u::Vector{Any} # Dummy
      interp::Function
   end

   (sol::MockSol)(t) = sol.interp(t)

   # 1. Test CartesianGrid (2D)
   @testset "CartesianGrid 2D" begin
      # Grid from (0,0) to (10,10) with 10x10 cells. Cell size 1x1. Volume = 1.
      grid = CartesianGrid((10, 10))

      # Particle 1: Stationary at (0.5, 0.5) -> Cell (1,1)
      p1 = MockSol(0.0:1.0:10.0, [], t -> SVector(0.5, 0.5, 0.0))
      # Particle 2: Moving along diagonal. At t=5, at (5.5, 5.5) -> Cell (6,6)
      p2 = MockSol(0.0:1.0:10.0, [], t -> SVector(t + 0.5, t + 0.5, 0.0))

      sols = [p1, p2]

      # Test Snapshot at t=0
      # p1 at (0.5, 0.5) -> (1,1)
      # p2 at (0.5, 0.5) -> (1,1)
      # Count at (1,1) should be 2. Density = 2 / 1 = 2.
      dens0 = TestParticle.get_number_density(sols, grid, 0.0)
      @test dens0[1, 1] ≈ 2.0
      @test sum(dens0) ≈ 2.0

      # Test Snapshot at t=5
      # p1 at (0.5, 0.5) -> (1,1)
      # p2 at (5.5, 5.5) -> (6,6)
      dens5 = TestParticle.get_number_density(sols, grid, 5.0)
      @test dens5[1, 1] ≈ 1.0
      @test dens5[6, 6] ≈ 1.0
      @test sum(dens5) ≈ 2.0

      # Test Time Average from t=0 to t=9 (10 steps)
      # p1 is always in (1,1). Adds 1 count for each of 10 steps. Total 10. Avg = 1.
      # p2 moves:
      # t=0: (0.5, 0.5) -> (1,1)
      # t=1: (1.5, 1.5) -> (2,2)
      # ...
      # t=9: (9.5, 9.5) -> (10,10)
      # Each diagonal cell gets 1 count. Avg = 1/10 = 0.1

      # (1,1) gets 10 from p1 and 1 from p2. Total 11. Avg = 1.1

      dens_avg = TestParticle.get_number_density(sols, grid, 0.0, 9.0, 1.0)
      @test dens_avg[1, 1] ≈ 1.1
      @test dens_avg[2, 2] ≈ 0.1
      @test dens_avg[10, 10] ≈ 0.1
      @test sum(dens_avg) ≈ 2.0 # Conservation of particles (averaged)
   end

   # 2. Test RectilinearGrid (2D) with non-uniform spacing
   @testset "RectilinearGrid 2D" begin
      # x: [0, 2, 4] -> 2 cells, widths 2 and 2
      # y: [0, 1, 5] -> 2 cells, heights 1 and 4
      # Cell (1,1): x[0,2], y[0,1]. Area = 2*1 = 2.
      # Cell (2,1): x[2,4], y[0,1]. Area = 2*1 = 2.
      # Cell (1,2): x[0,2], y[1,5]. Area = 2*4 = 8.
      # Cell (2,2): x[2,4], y[1,5]. Area = 2*4 = 8.

      x = [0.0, 2.0, 4.0]
      y = [0.0, 1.0, 5.0]
      grid = RectilinearGrid(x, y)

      # Particle in Cell (1,1). Pos (1.0, 0.5).
      p1 = MockSol(0.0:1.0:1.0, [], t -> SVector(1.0, 0.5, 0.0))

      sols = [p1]

      dens = TestParticle.get_number_density(sols, grid, 0.0)
      # Count is 1. Volume is 2. Density = 0.5.
      @test dens[1, 1] ≈ 0.5
      @test dens[1, 2] == 0.0

      # Particle in Cell (1,2). Pos (1.0, 3.0).
      p2 = MockSol(0.0:1.0:1.0, [], t -> SVector(1.0, 3.0, 0.0))
      sols = [p2]
      dens = TestParticle.get_number_density(sols, grid, 0.0)
      # Count is 1. Volume is 8. Density = 0.125.
      @test dens[1, 2] ≈ 0.125
   end

   # 3. Test CartesianGrid 3D
   @testset "CartesianGrid 3D" begin
      grid = CartesianGrid((2, 2, 2))
      # Volume = 1*1*1 = 1.

      p1 = MockSol(0.0:1.0:1.0, [], t -> SVector(0.5, 0.5, 0.5)) # (1,1,1)
      p2 = MockSol(0.0:1.0:1.0, [], t -> SVector(1.5, 1.5, 1.5)) # (2,2,2)

      sols = [p1, p2]
      dens = TestParticle.get_number_density(sols, grid, 0.0)

      @test dens[1, 1, 1] ≈ 1.0
      @test dens[2, 2, 2] ≈ 1.0
      @test dens[1, 2, 1] == 0.0
   end
end

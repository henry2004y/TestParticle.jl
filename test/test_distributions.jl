using TestParticle
using VelocityDistributionFunctions
using StaticArrays
using LinearAlgebra
using Statistics
using Test

@testset "Distributions" begin
   m = TestParticle.mᵢ

   # Maxwellian tests
   u0 = [10.0, 0.0, 0.0]
   vth = 1000.0
   n = 1e6
   maxwellian = Maxwellian(u0, vth^2 / 2 * m * n, n)
   @test maxwellian.vth ≈ vth
   v_m = rand(maxwellian)
   @test length(v_m) == 3

   # BiMaxwellian tests
   B = [1.0, 0.0, 0.0]
   vthpar = 1000.0
   vthperp = 500.0
   bimaxwellian = BiMaxwellian(B, u0, vthpar^2 / 2 * m * n, vthperp^2 / 2 * m * n, n)
   @test bimaxwellian.vth_para ≈ vthpar
   @test bimaxwellian.vth_perp ≈ vthperp
   v_bm = rand(bimaxwellian)
   @test length(v_bm) == 3

   # Kappa tests
   kappa = 4.0
   kdist = Kappa(u0, vth^2 / 2 * m * n, n, kappa)
   @test kdist.κ == kappa
   @test kdist.vth ≈ vth
   v_k = rand(kdist)
   @test length(v_k) == 3

   # BiKappa tests
   bikdist = BiKappa(B, u0, vthpar^2 / 2 * m * n, vthperp^2 / 2 * m * n, n, kappa)
   @test bikdist.κ == kappa
   @test bikdist.vth_para ≈ vthpar
   @test bikdist.vth_perp ≈ vthperp
   v_bk = rand(bikdist)
   @test length(v_bk) == 3

   # Statistical tests (variance check)
   N = 100000

   # Maxwellian Variance Check
   # Variance per component should be vth^2
   samples_m = [rand(maxwellian) - u0 for _ in 1:N]
   vars_m = [mean(2*v[i]^2 for v in samples_m) for i in 1:3]
   @test all(isapprox.(vars_m, vth^2, rtol = 0.05))

   # BiMaxwellian Variance Check
   # B is aligned with x
   samples_bm = [rand(bimaxwellian) - u0 for _ in 1:N]
   var_par_bm = mean(2*v[1]^2 for v in samples_bm) # parallel (x)
   vars_perp_bm = [mean(2*v[i]^2 for v in samples_bm) for i in 2:3] # perpendicular (y, z)
   @test isapprox(var_par_bm, vthpar^2, rtol = 0.05)
   @test all(isapprox.(vars_perp_bm, vthperp^2, rtol = 0.05))

   # Kappa Variance Check
   # Variance per component should be (2*kappa - 3) / (2*kappa) * vth^2 ??
   # Wait, the new package defines vth such that most probable speed v_p = vth * ...
   # Let's check the docstring or implementation.
   # _rand! uses scale = d.vth * sqrt(d.κ / ξ) where ξ ~ ChiSq(2κ - 1)
   # Variance of scale * Z (where Z ~ N(0,1)) is E[scale^2] * 1
   # E[scale^2] = vth^2 * κ * E[1/ξ]
   # If ξ ~ ChiSq(ν), E[1/ξ] = 1/(ν-2) for ν > 2
   # ν = 2κ - 1. So ν - 2 = 2κ - 3.
   # E[scale^2] = vth^2 * κ * (1 / (2κ - 3))
   # So variance = vth^2 * κ / (2κ - 3)

   theoretical_var_k = vth^2 * kappa / (2 * kappa - 3)

   samples_k = [rand(kdist) - u0 for _ in 1:N]
   var_k = mean(map(v -> v[1]^2, samples_k)) # x component
   @test isapprox(var_k, theoretical_var_k, rtol = 0.05)

   # BiKappa Variance Check
   # B is aligned with x
   # Similar to Kappa, we need to adjust expected variance based on how vth is defined.
   # Assuming vth in BiKappa follows same logic as Kappa but split.
   theoretical_var_par = vthpar^2 * kappa / (2 * kappa - 3)
   theoretical_var_perp = vthperp^2 * kappa / (2 * kappa - 3)

   samples_bk = [rand(bikdist) - u0 for _ in 1:N]
   var_par_bk = mean(map(v -> v[1]^2, samples_bk)) # parallel (x)
   var_perp_bk = mean(map(v -> v[2]^2, samples_bk)) # perpendicular (y)
   @test isapprox(var_par_bk, theoretical_var_par, rtol = 0.05)
   @test isapprox(var_perp_bk, theoretical_var_perp, rtol = 0.05)

   # Show methods check
   # Need to check expected output for new types
   # kdist is VelocityDistributionFunctions.Kappa
   @test occursin("Kappa", repr(kdist))
   @test startswith(repr(bikdist), "BiKappa")
end

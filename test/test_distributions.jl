using TestParticle
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
   maxwellian = Maxwellian(u0, vth * vth * m * n, n)
   @test maxwellian.vth ≈ vth

   # BiMaxwellian tests
   B = [1.0, 0.0, 0.0]
   vthpar = 1000.0
   vthperp = 500.0
   bimaxwellian = BiMaxwellian(B, u0, vthpar^2 * m * n, vthperp^2 * m * n, n)
   @test bimaxwellian.vthpar ≈ vthpar
   @test bimaxwellian.vthperp ≈ vthperp

   # Kappa tests
   kappa = 4.0
   kdist = Kappa(u0, vth * vth * m * n, n, kappa)
   @test kdist.kappa == kappa
   @test kdist.vth ≈ vth
   v_k = sample(kdist)
   @test length(v_k) == 3

   # BiKappa tests
   bikdist = BiKappa(B, u0, vthpar^2 * m * n, vthperp^2 * m * n, n, kappa)
   @test bikdist.kappa == kappa
   @test bikdist.vthpar ≈ vthpar
   @test bikdist.vthperp ≈ vthperp
   v_bk = sample(bikdist)
   @test length(v_bk) == 3

   # SelfSimilar tests
   s_exp = 3.0
   ssdist = SelfSimilar(u0, vth * vth * m * n, n, s_exp)
   @test ssdist.s == s_exp
   v_ss = sample(ssdist)
   @test length(v_ss) == 3

   # BiSelfSimilar tests
   p_exp = 3.0
   q_exp = 3.0
   bissdist = BiSelfSimilar(B, u0, vthpar^2 * m * n, vthperp^2 * m * n, n, p_exp, q_exp)
   @test bissdist.p_exp == p_exp
   @test bissdist.q_exp == q_exp
   v_bss = sample(bissdist)
   @test length(v_bss) == 3

   # Statistical tests (variance check)
   N = 200000

   # Kappa Variance Check
   # Variance per component should be vth^2
   samples_k = [sample(kdist) - u0 for _ in 1:N]
   var_k = mean(map(v -> v[1]^2, samples_k)) # x component
   @test isapprox(var_k, vth^2, rtol = 0.05)

   # BiKappa Variance Check
   # B is aligned with x
   samples_bk = [sample(bikdist) - u0 for _ in 1:N]
   var_par_bk = mean(map(v -> v[1]^2, samples_bk)) # parallel (x)
   var_perp_bk = mean(map(v -> v[2]^2, samples_bk)) # perpendicular (y)
   @test isapprox(var_par_bk, vthpar^2, rtol = 0.05)
   @test isapprox(var_perp_bk, vthperp^2, rtol = 0.05)

   # SelfSimilar Variance Check
   # Variance per component should be vth^2
   samples_ss = [sample(ssdist) - u0 for _ in 1:N]
   var_ss = mean(map(v -> v[1]^2, samples_ss)) # x component
   @test isapprox(var_ss, vth^2, rtol = 0.05)

   # BiSelfSimilar Variance Check
   # B aligned with x
   samples_bss = [sample(bissdist) - u0 for _ in 1:N]
   var_par_bss = mean(map(v -> v[1]^2, samples_bss)) # parallel
   var_perp_bss = mean(map(v -> v[2]^2, samples_bss)) # perpendicular
   @test isapprox(var_par_bss, vthpar^2, rtol = 0.05)
   @test isapprox(var_perp_bss, vthperp^2, rtol = 0.05) # Note: <v_perp^2> = 2 * vthperp^2, so component <v_y^2> = vthperp^2

   # Show methods check
   @test startswith(repr(kdist), "Isotropic Kappa distribution")
   @test startswith(repr(bikdist), "BiKappa distribution")
   @test startswith(repr(ssdist), "Isotropic Self-Similar distribution")
   @test startswith(repr(bissdist), "BiSelfSimilar distribution")
end

module test_distributions

using TestParticle
using StaticArrays
using LinearAlgebra
using Statistics
using Test
# Loading the backend package is what attaches the distribution constructors.
using VelocityDistributionFunctions
import TestParticle: Maxwellian, BiMaxwellian, Kappa, BiKappa

const m = TestParticle.mᵢ
const N = 100000

"Variance of the `i`-th component around `u0`, estimated from `samples`."
component_var(samples, i, u0) = mean((v[i] - u0[i])^2 for v in samples)

@testset "Distributions" begin
    u0 = [10.0, 0.0, 0.0]
    vth = 1000.0
    n = 1.0e6
    # B is aligned with x, so x is the parallel direction
    B = [1.0, 0.0, 0.0]
    vthpar = 1000.0
    vthperp = 500.0
    kappa = 4.0
    # `_rand!` draws scale = vth * sqrt(κ / ξ) with ξ ~ ChiSq(2κ - 1), and
    # E[1/ξ] = 1/(2κ - 3), so the per-component variance is vth² κ / (2κ - 3).
    κ_var = kappa / (2 * kappa - 3)

    # Maxwellian
    maxwellian = Maxwellian(u0, vth^2 / 2 * m * n, n)
    @test maxwellian.vth ≈ vth
    @test length(rand(maxwellian)) == 3

    # BiMaxwellian
    bimaxwellian = BiMaxwellian(B, u0, vthpar^2 / 2 * m * n, vthperp^2 / 2 * m * n, n)
    @test bimaxwellian.vth_para ≈ vthpar
    @test bimaxwellian.vth_perp ≈ vthperp
    @test length(rand(bimaxwellian)) == 3

    # Kappa
    kdist = Kappa(u0, vth^2 / 2 * m * n, n, kappa)
    @test kdist.κ == kappa
    @test kdist.vth ≈ vth
    @test length(rand(kdist)) == 3

    # BiKappa
    bikdist = BiKappa(B, u0, vthpar^2 / 2 * m * n, vthperp^2 / 2 * m * n, n, kappa)
    @test bikdist.κ == kappa
    @test bikdist.vth_para ≈ vthpar
    @test bikdist.vth_perp ≈ vthperp
    @test length(rand(bikdist)) == 3

    # Standard VDF constructors (forwarders)
    @test Maxwellian(u0, vth).vth ≈ vth
    bimaxwellian_std = BiMaxwellian(vthperp, vthpar, B; u0)
    @test bimaxwellian_std.vth_para ≈ vthpar
    @test bimaxwellian_std.vth_perp ≈ vthperp
    kappa_std = Kappa(vth, kappa; u0)
    @test kappa_std.κ == kappa
    @test kappa_std.vth ≈ vth
    bikappa_std = BiKappa(vthperp, vthpar, kappa, B; u0)
    @test bikappa_std.κ == kappa
    @test bikappa_std.vth_para ≈ vthpar
    @test bikappa_std.vth_perp ≈ vthperp

    # Maxwellian: the variance of every component is vth²/2
    samples_m = [rand(maxwellian) for _ in 1:N]
    for i in 1:3
        @test component_var(samples_m, i, u0) ≈ vth^2 / 2 rtol = 0.05
    end

    # BiMaxwellian: parallel along x, perpendicular along y and z
    samples_bm = [rand(bimaxwellian) for _ in 1:N]
    @test component_var(samples_bm, 1, u0) ≈ vthpar^2 / 2 rtol = 0.05
    for i in 2:3
        @test component_var(samples_bm, i, u0) ≈ vthperp^2 / 2 rtol = 0.05
    end

    samples_k = [rand(kdist) for _ in 1:N]
    @test component_var(samples_k, 1, u0) ≈ vth^2 * κ_var rtol = 0.05

    samples_bk = [rand(bikdist) for _ in 1:N]
    @test component_var(samples_bk, 1, u0) ≈ vthpar^2 * κ_var rtol = 0.05
    @test component_var(samples_bk, 2, u0) ≈ vthperp^2 * κ_var rtol = 0.05

    # Show methods
    @test occursin("Kappa", repr(kdist))
    @test occursin("BiKappa", repr(bikdist))
end

end # module test_distributions

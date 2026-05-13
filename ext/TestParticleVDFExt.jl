module TestParticleVDFExt

using TestParticle
import VelocityDistributionFunctions as VDF
using LinearAlgebra: norm
using StaticArrays: SA
using Random: AbstractRNG, default_rng
using PrecompileTools: @setup_workload, @compile_workload

function TestParticle.sample_maxwellian(rng::AbstractRNG, Tn, m; offset = 0.0, u0 = SA[0.0, 0.0, 0.0])
    vth = √(2 * TestParticle.kB * Tn / m)
    vdf = VDF.Maxwellian(vth; u0)
    v = rand(rng, vdf)
    if offset > 0
        v_mag = norm(v)
        E = 0.5 * m * v_mag^2 + offset
        v_mag_new = sqrt(2 * E / m)
        v *= (v_mag_new / v_mag)
    end
    return v
end

function TestParticle.sample_maxwellian(Tn::Real, m::Real; offset = 0.0, u0 = SA[0.0, 0.0, 0.0])
    return TestParticle.sample_maxwellian(default_rng(), Tn, m; offset, u0)
end

function TestParticle.sample_maxwellian(rng::AbstractRNG, Tn::Real, species::String; offset = 0.0, u0 = SA[0.0, 0.0, 0.0])
    haskey(TestParticle.SpeciesDict, species) || error("Unsupported species: $species")
    s = TestParticle.SpeciesDict[species]
    return TestParticle.sample_maxwellian(rng, Tn, s.m; offset, u0)
end

function TestParticle.sample_maxwellian(Tn::Real, species::String; offset = 0.0, u0 = SA[0.0, 0.0, 0.0])
    return TestParticle.sample_maxwellian(default_rng(), Tn, species; offset, u0)
end

function TestParticle.Maxwellian(args...; kwargs...)
    return VDF.Maxwellian(args...; kwargs...)
end

function TestParticle.Maxwellian(u0, p, n; m = TestParticle.mᵢ)
    vth = TestParticle.get_thermal_speed(p, n, m)

    return VDF.Maxwellian(vth; u0)
end

function TestParticle.BiMaxwellian(args...; kwargs...)
    return VDF.BiMaxwellian(args...; kwargs...)
end

function TestParticle.BiMaxwellian(B, u0, ppar, pperp, n; m = TestParticle.mᵢ)
    vpar = TestParticle.get_thermal_speed(ppar, n, m)
    vperp = TestParticle.get_thermal_speed(pperp, n, m)

    return VDF.BiMaxwellian(vperp, vpar, B; u0)
end

TestParticle.Kappa(args...; kwargs...) = VDF.Kappa(args...; kwargs...)

function TestParticle.Kappa(u0, p, n, kappa; m = TestParticle.mᵢ)
    vth = TestParticle.get_thermal_speed(p, n, m)

    return VDF.Kappa(vth, kappa; u0)
end

TestParticle.BiKappa(args...; kwargs...) = VDF.BiKappa(args...; kwargs...)

function TestParticle.BiKappa(B, u0, ppar, pperp, n, kappa; m = TestParticle.mᵢ)
    vpar = TestParticle.get_thermal_speed(ppar, n, m)
    vperp = TestParticle.get_thermal_speed(pperp, n, m)

    return VDF.BiKappa(vperp, vpar, kappa, B; u0)
end

@setup_workload begin
    @compile_workload begin
        u0 = [0.0, 0.0, 0.0]
        p = 1.0e-9
        n = 1.0e6
        vdf = Maxwellian(u0, p, n)
        v = rand(vdf)
        B = [1.0, 0.0, 0.0]
        vdf = BiMaxwellian(B, u0, p, p, n)
        v = rand(vdf)
        kappa = 4.0
        vdf = Kappa(u0, p, n, kappa)
        v = rand(vdf)
        vdf = BiKappa(B, u0, p, p, n, kappa)
        v = rand(vdf)
        TestParticle.sample_maxwellian(3000.0, "O+")
    end
end

end

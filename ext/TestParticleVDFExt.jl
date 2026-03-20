module TestParticleVDFExt

using TestParticle
import TestParticle as TP
import VelocityDistributionFunctions as VDF
import TestParticle: Maxwellian, BiMaxwellian, Kappa, BiKappa
using LinearAlgebra: norm
using StaticArrays: SA
using PrecompileTools: @setup_workload, @compile_workload

const species_name_to_type = Dict(
    "H+" => TP.Proton,
    "Proton" => TP.Proton,
    "e-" => TP.Electron,
    "Electron" => TP.Electron,
    "O+" => TP.Ion(16, 1),
    "He+" => TP.Ion(4, 1),
    "O2+" => TP.Ion(32, 1),
    "CO2+" => TP.Ion(44, 1),
)

function TP.sample_maxwellian(Tn, m; offset = 0.0, u0 = SA[0.0, 0.0, 0.0])
    vth = √(2 * TP.kB * Tn / m)
    vdf = VDF.Maxwellian(vth; u0)
    v = rand(vdf)
    if offset > 0
        v_mag = norm(v)
        E = 0.5 * m * v_mag^2 + offset
        v_mag_new = sqrt(2 * E / m)
        v *= (v_mag_new / v_mag)
    end
    return v
end

function TP.sample_maxwellian(Tn, species::String; offset = 0.0, u0 = SA[0.0, 0.0, 0.0])
    haskey(species_name_to_type, species) || error("Unsupported species: $species")
    s = species_name_to_type[species]
    return TP.sample_maxwellian(Tn, s.m; offset, u0)
end

function Maxwellian(args...; kwargs...)
    return VDF.Maxwellian(args...; kwargs...)
end

function Maxwellian(u0, p, n; m = TP.mᵢ)
    vth = TP.get_thermal_speed(p, n, m)

    return VDF.Maxwellian(vth; u0)
end

function BiMaxwellian(args...; kwargs...)
    return VDF.BiMaxwellian(args...; kwargs...)
end

function BiMaxwellian(B, u0, ppar, pperp, n; m = TP.mᵢ)
    vpar = TP.get_thermal_speed(ppar, n, m)
    vperp = TP.get_thermal_speed(pperp, n, m)

    return VDF.BiMaxwellian(vperp, vpar, B; u0)
end

Kappa(args...; kwargs...) = VDF.Kappa(args...; kwargs...)

function Kappa(u0, p, n, kappa; m = TP.mᵢ)
    vth = TP.get_thermal_speed(p, n, m)

    return VDF.Kappa(vth, kappa; u0)
end

BiKappa(args...; kwargs...) = VDF.BiKappa(args...; kwargs...)

function BiKappa(B, u0, ppar, pperp, n, kappa; m = TP.mᵢ)
    vpar = TP.get_thermal_speed(ppar, n, m)
    vperp = TP.get_thermal_speed(pperp, n, m)

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
        sample_maxwellian(3000.0, "O+")
    end
end

end

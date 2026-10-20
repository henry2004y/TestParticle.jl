# CODATA 2018 values
const qᵢ = 1.602176634e-19 # proton charge, [C]
const qₑ = -qᵢ # electron charge, [C]
const mₑ = 9.1093837015e-31 # electron mass, [kg]
const mᵢ = 1.67262192369e-27 # proton mass, [kg]
const c = 2.99792458e8 # speed of light, [m/s]
const μ₀ = 1.25663706212e-6 # Vacuum permeability, [H/m]
const ϵ₀ = 8.8541878128e-12 # Vacuum permittivity, [F/m]
const kB = 1.380649e-23 # Boltzmann constant, [J/K]
const Rₑ = 6371.0e3 # Earth radius, [m]
const BMoment_Earth = SVector(0.0, 0.0, -7.94e22) # [V*s/(A*m)]
const eV = qᵢ # electron volt, [J]

const Proton = Species(mᵢ, qᵢ)
const Electron = Species(mₑ, qₑ)

Ion(m, q = 1) = Species(m * mᵢ, q * qᵢ)
Ion(; m = 1, q = 1) = Species(m * mᵢ, q * qᵢ)

const SpeciesDict = Dict(
    "H+" => Proton,
    "Proton" => Proton,
    "e-" => Electron,
    "Electron" => Electron,
    "O+" => Ion(16, 1),
    "He+" => Ion(4, 1),
    "O2+" => Ion(32, 1),
    "CO2+" => Ion(44, 1),
)
